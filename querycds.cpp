// =============================================================================
// querycds.cpp  —  Stage 3: C++ query engine + benchmark timer
// =============================================================================
//
//  IT loads precomputed.bin, auto-detects which
//   construction mode was used (3 or 5), then evaluates a random Straight-Line
//   Program (SLP) from structure.txt using the compact data structure.
//   Timing is wall-clock, fixed repetitions, with an untimed warm-up phase.
//   Results are appended to benchmark_results_c++.csv.
//
// COMPILE
//   g++ -O3 -o querycds querycds.cpp
//
// RUN
//   ./querycds <group_name> <L> <reps> <seed> <warmup_reps>
//     group_name  : string label written to CSV (e.g. "SmallGroup(64,267)")
//     L           : SLP size (used for per-operation time calculation)
//     reps        : timed repetition count (default 50)
//     seed        : which SLP seed was used (for CSV tracking, not actual RNG)
//     warmup_reps : untimed warm-up repetitions (default 5)
//
// MODES
//   Mode 3 — Theorem 3 (Case 1)
//     11 array lookups per multiply, branch-free
//
//   Mode 5 — Theorem 4 nested inside Theorem 3 (Case 2)
//     29 array lookups per multiply (9 outer + 2x10 inner T4 calls)
//
//   Mode 4 — Standalone Theorem 4 (rarely used, kept for completeness)
// =============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
using namespace std;

// =============================================================================
// Binary block reader
// =============================================================================
struct Block {
    string tag;
    int rows = 0;
    int cols = 0;
    vector<int> data;
};

static string read_tag(ifstream &in) {
    int taglen;
    if (!in.read(reinterpret_cast<char*>(&taglen), sizeof(taglen))) return "";
    if (taglen <= 0) return "";
    string tag(taglen, '\0');
    in.read(&tag[0], taglen);
    return tag;
}

static Block read_block(ifstream &in) {
    Block b;
    b.tag = read_tag(in);
    if (b.tag.empty()) return b;
    in.read(reinterpret_cast<char*>(&b.rows), sizeof(b.rows));
    in.read(reinterpret_cast<char*>(&b.cols), sizeof(b.cols));
    // For 1-D arrays cols==0, so we read rows*1 ints.
    // For 2-D matrices we read rows*cols ints.
    int count = b.rows * max(1, b.cols);
    b.data.resize(count);
    if (count > 0)
        in.read(reinterpret_cast<char*>(b.data.data()), sizeof(int) * count);
    return b;
}

// =============================================================================
// Mode detection
// preparedsa writes "e_arr" as a block tag only in Case 2 (Mode 5).
// =============================================================================
int detect_mode(const unordered_map<string, Block>& B) {
    bool hasT4 = B.count("e_arr") && B.count("Flip_mat") && B.count("rede_arr");
    bool hasT3 = B.count("cL")    && B.count("cR")       && B.count("flipH_mat");
    if (hasT4 && hasT3) return 5;
    if (hasT4)          return 4;
    if (hasT3)          return 3;
    return -1;
}

int detect_orderG(const unordered_map<string, Block>& B) {
    if (B.count("cL"))    return B.at("cL").rows;
    if (B.count("e_arr")) return B.at("e_arr").rows;
    return -1;
}

// =============================================================================
// Theorem data structs containing required tables and arrays blocks
// =============================================================================
struct Theorem3Data {
    const Block *cL, *sL, *cR, *sR;         // decomposition arrays
    const Block *HH, *flipH, *flipR;         // flip step uses HH twice
    const Block *crossH, *crossR;            // cross step
    const Block *fuse, *inverse;             // reassembly and inverse
};

struct Theorem4Data {
    const Block *e_arr, *sR_arr, *sL_arr;   // exponent-residue decomposition
    const Block *NN, *Flip;                  // inner-N multiply and flip
    const Block *rede, *redN;                // exponent reduction tables
    const Block *Fuse, *inverse;             // reassembly and inverse
};

// =============================================================================
// multiply_theorem3
// =============================================================================
inline int multiply_theorem3(int g1, int g2, const Theorem3Data& T) {
    int cl1 = T.cL->data[g1],  sl1 = T.sL->data[g1];
    int sr1 = T.sR->data[g1],  r2  = T.cR->data[g2];
    // Step 2: HH[sL(g1)][sR(g2)]
    int h1  = T.HH->data[sl1 * T.HH->cols + T.sR->data[g2]];
    // Step 3: flip
    int h2  = T.flipH->data[cl1 * T.flipH->cols + h1];
    int r1  = T.flipR->data[cl1 * T.flipR->cols + h1];
    // Step 4: cross
    int h3  = T.crossH->data[r1 * T.crossH->cols + r2];
    int r3  = T.crossR->data[r1 * T.crossR->cols + r2];
    // Step 5: second HH
    int h4  = T.HH->data[h2 * T.HH->cols + h3];
    // Step 6: fuse -> global G-index of the result
    return T.fuse->data[h4 * T.fuse->cols + r3];
}

// =============================================================================
// multiply_theorem4
// =============================================================================
inline int multiply_theorem4(int g1, int g2, const Theorem4Data& T) {
    int alpha = T.e_arr->data[g1],  beta  = T.e_arr->data[g2];
    int sR_g1 = T.sR_arr->data[g1], sL_g2 = T.sL_arr->data[g2];
    // Step 2
    int n1    = T.NN->data[sR_g1 * T.NN->cols + sL_g2];
    // Step 3: flip n1 past g0^beta
    int n2    = T.Flip->data[n1 * T.Flip->cols + beta];
    // Step 4: reduce exponent
    int ell   = alpha + beta;
    int gamma = T.rede->data[ell],  redn = T.redN->data[ell];
    // Step 5
    int n3    = T.NN->data[redn * T.NN->cols + n2];
    // Step 6: reassemble
    return T.Fuse->data[gamma * T.Fuse->cols + n3];
}

// =============================================================================
// multiply_combined  —  Mode 5: T3 outer algorithm with T4 as inner algorithm
// =============================================================================
inline int multiply_combined(int g1, int g2,
                              const Theorem3Data& T3,
                              const Theorem4Data& T4,
                              const vector<int>& subgroupH) {
    int cl1 = T3.cL->data[g1], sl1 = T3.sL->data[g1];
    int sl2 = T3.sR->data[g2], r2  = T3.cR->data[g2];

    // Step 2 replacement: inner multiply sl1 * sl2 via T4.
    // sl1 and sl2 are local H-indices; T4 needs global G-indices.
    int h_global_1 = subgroupH[sl1];
    int h_global_2 = subgroupH[sl2];
    int h1_global  = multiply_theorem4(h_global_1, h_global_2, T4);
    // Convert T4 result (a global G-index) back to a local H-index.
    int h1 = T3.sR->data[h1_global];

    // Step 3: flip (T3's outer tables)
    int h2 = T3.flipH->data[cl1 * T3.flipH->cols + h1];
    int r1 = T3.flipR->data[cl1 * T3.flipR->cols + h1];
    // Step 4: cross (T3's outer tables)
    int h3 = T3.crossH->data[r1 * T3.crossH->cols + r2];
    int r3 = T3.crossR->data[r1 * T3.crossR->cols + r2];

    // Step 5 replacement: inner multiply h2 * h3 via T4.
    int hg2 = subgroupH[h2], hg3 = subgroupH[h3];
    int h4_global = multiply_theorem4(hg2, hg3, T4);
    int h4 = T3.sR->data[h4_global];

    // Step 6: fuse (T3's outer table)
    return T3.fuse->data[h4 * T3.fuse->cols + r3];
}

// =============================================================================
// SLP loader
//
// Reads structure.txt produced by structure_generator. Each line is one
// instruction in the Straight-Line Program:
//   C <id>      — load constant: result = group element with global index id
//   I <a>       — inversion:     result = inverse of instruction a's result
//   M <a> <b>   — multiply:      result = (instruction a's result) * (instruction b's result)
// =============================================================================
enum OpType { CONST_OP, INV_OP, MUL_OP };
struct Instruction { OpType type; int a, b; };

vector<Instruction> load_structure(const string& filename, int& maxLeaf) {
    ifstream in(filename);
    string line;
    vector<Instruction> ops;
    maxLeaf = -1;

    getline(in, line);  // skip "L <size>" header line

    while (getline(in, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        char type; ss >> type;

        if      (type == 'C') { int id;    ss >> id;    ops.push_back({CONST_OP, id, -1}); maxLeaf = max(maxLeaf, id); }
        else if (type == 'I') { int a;     ss >> a;     ops.push_back({INV_OP,   a,  -1}); }
        else if (type == 'M') { int a, b;  ss >> a >> b; ops.push_back({MUL_OP,  a,   b}); }
    }
    return ops;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[]) {

    // ---- Parse command-line arguments ---------------------------------------
    // All arguments are optional — defaults make it runnable without run_all.sh.
    std::string groupName = "UNKNOWN";
    int L           = -1;
    int fixed_reps  = 50;  // fixed rep count, matches FIXED_REPS in run_all.sh
    int slp_seed    = 0;   // only used for CSV tracking, not for seeding
    int warmup_reps = 5;   // untimed reps before the clock starts

    if (argc >= 2) groupName   = argv[1];
    if (argc >= 3) L           = std::stoi(argv[2]);
    if (argc >= 4) fixed_reps  = std::stoi(argv[3]);
    if (argc >= 5) slp_seed    = std::stoi(argv[4]);
    if (argc >= 6) warmup_reps = std::stoi(argv[5]);

    // structure.txt is always read from the current working directory.
    // run_all.sh regenerates it for each (group, L, seed) combination.
    const string STRUCTURE_FILE = "structure.txt";

    // ---- Load precomputed group data from binary file ----------------------
    // All blocks land in an unordered_map keyed by tag name.
    ifstream bin("precomputed.bin", ios::binary);
    unordered_map<string, Block> blocks;
    while (bin.peek() != EOF) {
        Block b = read_block(bin);
        if (b.tag.empty()) break;
        blocks.emplace(b.tag, move(b));
    }
    bin.close();

    int mode   = detect_mode(blocks);
    int orderG = detect_orderG(blocks);
    if (mode == -1 || orderG <= 0) {
        cerr << "Invalid or unrecognised group data in precomputed.bin\n";
        return 1;
    }

    // ---- Set up Theorem data struct pointers --------------------------------
    // These just point into the blocks map — no data is copied.
    Theorem3Data T3{}; Theorem4Data T4{};

    if (mode == 3 || mode == 5) {
        // Mode 3: HH_mat is present and read twice per multiply.
        // Mode 5: HH_mat is not written (skipped to keep O(n) space).
        //   multiply_combined never reads T3.HH — T4 handles those steps.
        T3 = { &blocks["cL"], &blocks["sL"], &blocks["cR"], &blocks["sR"],
               (mode == 3 ? &blocks["HH_mat"] : nullptr),
               &blocks["flipH_mat"], &blocks["flipR_mat"],
               &blocks["crossH_mat"], &blocks["crossR_mat"],
               &blocks["fuse_mat"], &blocks["inverse"] };
    }
    if (mode == 4 || mode == 5) {
        // Theorem 4 tables present in Mode 4 and Mode 5.
        T4 = { &blocks["e_arr"], &blocks["sR_arr"], &blocks["sL_arr"],
               &blocks["NN_mat"], &blocks["Flip_mat"],
               &blocks["rede_arr"], &blocks["redN_idx"],
               &blocks["Fuse_mat"], &blocks["inverse"] };
    }

    // In Mode 5 we need the global element list of the inner subgroup H = Gi
    // to convert between local H-indices and global G-indices inside
    // multiply_combined (Steps 2 and 5 of the outer Theorem 3 algorithm).
    vector<int> subgroupH;
    if (mode == 5)
        subgroupH = blocks["subgroupH"].data;

    // ---- Load SLP -----------------------------------------------------------
    int maxLeaf;
    vector<Instruction> ops = load_structure(STRUCTURE_FILE, maxLeaf);
    cout << "Mode          : " << mode       << "\n";
    cout << "Instructions  : " << ops.size() << "\n";
    cout << "Max leaf value: " << maxLeaf     << "\n";
    cout << "Group order   : " << orderG      << "\n";

    // ---- Warm-up (NOT timed) -----------------------------------------------
    // We run warmup_reps full SLP evaluations before starting the clock.
    
    {
        vector<int> wvals(ops.size());
        volatile int wsink = 0;
        for (int r = 0; r < warmup_reps; r++) {
            for (size_t i = 0; i < ops.size(); i++) {
                const Instruction& op = ops[i];
                if (op.type == CONST_OP) {
                    wvals[i] = op.a;
                } else if (op.type == INV_OP) {
                    int v = wvals[op.a];
                    wvals[i] = (mode == 3) ? T3.inverse->data[v] : T4.inverse->data[v];
                } else {
                    int v1 = wvals[op.a], v2 = wvals[op.b];
                    if      (mode == 3) wvals[i] = multiply_theorem3  (v1, v2, T3);
                    else if (mode == 4) wvals[i] = multiply_theorem4  (v1, v2, T4);
                    else                wvals[i] = multiply_combined  (v1, v2, T3, T4, subgroupH);
                }
            }
            wsink = wvals.back();
        }
        (void)wsink;
    }

    // ---- Timed SLP evaluation (fixed repetitions) --------------------------
    // We repeat the full SLP fixed_reps times and record total wall-clock time.
    vector<int> values(ops.size());
    int repetitions = fixed_reps;
    volatile int sink = 0;

    auto tstart = std::chrono::steady_clock::now();
    for (int r = 0; r < repetitions; r++) {
        for (size_t i = 0; i < ops.size(); i++) {
            const Instruction& op = ops[i];
            if (op.type == CONST_OP) {
                values[i] = op.a;
            } else if (op.type == INV_OP) {
                int v = values[op.a];
                // Both modes store inverse under the same tag name in the binary.
                values[i] = (mode == 3) ? T3.inverse->data[v] : T4.inverse->data[v];
            } else {
                int v1 = values[op.a], v2 = values[op.b];
                if      (mode == 3) values[i] = multiply_theorem3  (v1, v2, T3);
                else if (mode == 4) values[i] = multiply_theorem4  (v1, v2, T4);
                else                values[i] = multiply_combined  (v1, v2, T3, T4, subgroupH);
            }
        }
        sink = values.back();
    }
    auto tend = std::chrono::steady_clock::now();
    chrono::duration<double, milli> duration = tend - tstart;

    int rootValue = values.back();  // final result element index (for cross-checking with GAP)

    // ---- SLP statistics (sanity check) -------------------------------------
    int num_const = 0, num_inv = 0, num_mul = 0;
    for (const auto& op : ops) {
        if      (op.type == CONST_OP) num_const++;
        else if (op.type == INV_OP)   num_inv++;
        else if (op.type == MUL_OP)   num_mul++;
    }
    cout << "\n--- Tree Statistics ---\n";
    cout << "Total nodes      : " << ops.size() << "\n";
    cout << "Constants (C)    : " << num_const << " (" << (100.0 * num_const / ops.size()) << "%)\n";
    cout << "Inversions (I)   : " << num_inv   << " (" << (100.0 * num_inv   / ops.size()) << "%)\n";
    cout << "Multiplications  : " << num_mul   << " (" << (100.0 * num_mul   / ops.size()) << "%)\n";
    cout << "C - M            : " << (num_const - num_mul) << " (must be 1 for a valid binary tree)\n";

    cout << "\n--- Results ---\n";
    cout << "Root value       : " << rootValue << "\n";
    cout << "Total time (ms)  : " << duration.count() << "\n";
    cout << "Avg per rep (ms) : " << duration.count() / repetitions << "\n";

    // ---- Append to CSV -----------------------------------------------------
    // If the file doesn't exist yet, write the header first.
    // The CSV is appended to (not overwritten) so run_all.sh can process
    // all groups in one run with a single cumulative file.
    string filepath = "benchmark_results_c++.csv";
    bool fileExists = false;
    { ifstream chk(filepath); fileExists = chk.good(); }

    ofstream out(filepath, ios::app);
    if (!fileExists)
        out << "Engine,Group,Order,Mode,Instructions,Seed,Repetitions,Result,TotalTime_ms,TimePerRep_ms\n";

    out << "C++"
        << ",\"" << groupName << "\""
        << "," << orderG
        << "," << mode
        << "," << ops.size()
        << "," << slp_seed
        << "," << repetitions
        << "," << rootValue
        << "," << duration.count()
        << "," << duration.count() / repetitions
        << "\n";

    cout << "\nC++ results written to CSV.\n";
    return 0;
}