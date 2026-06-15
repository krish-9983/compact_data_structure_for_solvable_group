<<<<<<< HEAD
=======
// =============================================================================
// structure_generator.cpp  —  Random SLP (binary expression tree) generator
// =============================================================================
//
// PURPOSE
//   Generates a random binary expression tree — also called a Straight-Line
//   Program (SLP) in group theory — and writes it to structure.txt. Both
//   querycds and tree_benchmark.g read the same structure.txt, so both engines
//   execute the exact same sequence of group operations.
//
// COMPILE
//   g++ -O3 -o structure_generator structure_generator.cpp
//
// RUN
//   ./structure_generator <L> <group_order_minus_1> [seed]
//     L                   : total number of nodes (leaves + internal)
//     group_order_minus_1 : maximum group element index (= |G| - 1)
//     seed                : mt19937 seed for reproducibility
//


// =============================================================================

>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <set>
using namespace std;

enum OpType { CONST_OP, MUL_OP, INV_OP };

struct Instruction {
    OpType type;
<<<<<<< HEAD
    int a;
    int b;
};

 int main(int argc, char* argv[]) {

    int L = 1000;      // default number of nodes
    int order = 1151;  // default (group_order_minus_1)
    int seed  = 2;     // default seed — pass different seeds for different SLPs
=======
    int a;   // for CONST: element index; for INV/MUL: operand instruction index
    int b;   // only used for MUL (second operand); -1 otherwise
};

int main(int argc, char* argv[]) {

    // Default values — useful for quick ad-hoc testing without run_all.sh.
    int L     = 1000;   // number of nodes in the SLP
    int order = 1151;   // group_order - 1 (max valid element index)
    int seed  = 2;      // mt19937 seed
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033

    if (argc == 4) {
        L     = atoi(argv[1]);
        order = atoi(argv[2]);
        seed  = atoi(argv[3]);
<<<<<<< HEAD
    }
    else if (argc == 3) {
        L     = atoi(argv[1]);
        order = atoi(argv[2]);
    }
    else if (argc != 1) {
=======
    } else if (argc == 3) {
        L     = atoi(argv[1]);
        order = atoi(argv[2]);
    } else if (argc != 1) {
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
        cerr << "Usage: ./structure_generator <L> <group_order_minus_1> [seed]\n";
        cerr << "  seed controls which random SLP is generated.\n";
        cerr << "  Same seed + increasing L = prefix-locked growing SLPs.\n";
        return 1;
    }

<<<<<<< HEAD
    if (L < 2) {
        cerr << "Error: L must be >= 2\n";
        return 1;
    }
    
    // For a VALID BINARY TREE: C - M = 1
    int num_inv_target   = L / 4;
    int num_mul_target   = (L - 1 - num_inv_target); 
    num_mul_target   = (L - 1 - num_inv_target) / 2;
    int num_const_target = num_mul_target + 1;  // ALWAYS +1, not recalculated
num_inv_target = L - num_const_target - num_mul_target;

    cout << "Target distribution (VALID TREE):\n";
    cout << "  Constants      : " << num_const_target <<  "\n";
    cout << "  Inversions     : " << num_inv_target   << "\n";
    cout << "  Multiplications: " << num_mul_target   << "\n";
    cout << "  Tree constraint: C - M = " << (num_const_target - num_mul_target) 
         << "\n\n";

    mt19937 rng(seed); // seed passed as 3rd argument — same seed + increasing L gives prefix-locked SLPs
    uniform_int_distribution<> const_dist(0, order-1);
=======
    if (L < 2) { cerr << "Error: L must be >= 2\n"; return 1; }

    // ---- Compute target instruction counts ----------------------------------
    // For a valid binary tree: constants - multiplications = 1.
    // We target roughly 25% inversions and split the rest between constants
    // and multiplications (keeping constants : +1).
    
    int num_mul_target   = (L - 1 - num_inv_target) / 2;
    int num_const_target = num_mul_target + 1;  // +1 is the binary tree invariant
    num_inv_target = L - num_const_target - num_mul_target;  // fill remainder

    cout << "Target distribution (VALID TREE):\n";
    cout << "  Constants      : " << num_const_target << "\n";
    cout << "  Inversions     : " << num_inv_target   << "\n";
    cout << "  Multiplications: " << num_mul_target   << "\n";
    cout << "  Tree constraint: C - M = " << (num_const_target - num_mul_target) << "\n\n";

    // mt19937 with a fixed seed gives a reproducible sequence.
    // The SAME seed with a LARGER L extends the sequence.
    mt19937 rng(seed);
    uniform_int_distribution<> const_dist(0, order - 1);  // random group element index
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033

    vector<Instruction> program;
    program.reserve(L);

<<<<<<< HEAD
    vector<int> available;

    int num_const = 0;
    int num_inv = 0;
    int num_mul = 0;

    // Step 1: Create two initial constants
=======
    // available: instruction indices whose results are currently "available" as
    // operands. When an instruction consumes a result (for MUL or INV), we
    // remove it from available so it's not reused. When an instruction produces
    // a result, we add its index. This maintains the single-rooted tree structure.
    vector<int> available;

    int num_const = 0, num_inv = 0, num_mul = 0;

    // ---- Step 1: seed with two initial constants ----------------------------
    // We need at least 2 available operands before we can generate a MUL.
    // Starting with 2 constants ensures the first real instruction has choices.
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
    for (int i = 0; i < 2; i++) {
        int value = const_dist(rng);
        program.push_back({CONST_OP, value, -1});
        available.push_back(i);
        num_const++;
    }

<<<<<<< HEAD
    // Step 2: Build tree with target ratios
    for (int i = 2; i < L; i++) {

        int const_needed = num_const_target - num_const;
        int inv_needed   = num_inv_target - num_inv;
        int mul_needed   = num_mul_target - num_mul;

        // Build list of valid choices
        vector<pair<OpType, int>> choices;
        
        if (const_needed > 0) {
            choices.push_back({CONST_OP, const_needed});
        }
        
        if (inv_needed > 0 && available.size() >= 1) {
            choices.push_back({INV_OP, inv_needed});
        }
        
        if (mul_needed > 0 && available.size() >= 2) {
            choices.push_back({MUL_OP, mul_needed});
        }

        // Choose operation
        OpType chosen_op;
        
        if (choices.empty()) {
            chosen_op = CONST_OP;
        }
        else if (choices.size() == 1) {
            chosen_op = choices[0].first;
        }
        else {
            // Weighted random selection
            int total_weight = 0;
            for (auto& p : choices) total_weight += p.second;
            
            uniform_int_distribution<> choice_dist(0, total_weight - 1);
            int r = choice_dist(rng);
            
            int cumsum = 0;
            chosen_op = CONST_OP;
            for (auto& p : choices) {
                cumsum += p.second;
                if (r < cumsum) {
                    chosen_op = p.first;
                    break;
                }
            }
        }

        // Execute the chosen operation
=======
    // ---- Step 2: build the remaining L-2 instructions ----------------------
    for (int i = 2; i < L; i++) {
        // How many of each type do we still need?
        int const_needed = num_const_target - num_const;
        int inv_needed   = num_inv_target   - num_inv;
        int mul_needed   = num_mul_target   - num_mul;

        // Build a list of valid choices at this point.
        // INV needs >= 1 available operand; MUL needs >= 2.
        vector<pair<OpType, int>> choices;
        if (const_needed > 0)                          choices.push_back({CONST_OP, const_needed});
        if (inv_needed   > 0 && available.size() >= 1) choices.push_back({INV_OP,   inv_needed});
        if (mul_needed   > 0 && available.size() >= 2) choices.push_back({MUL_OP,   mul_needed});

        // Weighted random selection: instructions we need more of are
        // proportionally more likely to be chosen. This keeps the actual
        // counts close to the targets throughout the build rather than
        // overshooting and having to scramble at the end.
        OpType chosen_op = CONST_OP;
        if (choices.empty()) {
            chosen_op = CONST_OP;
        } else if (choices.size() == 1) {
            chosen_op = choices[0].first;
        } else {
            int total_weight = 0;
            for (auto& p : choices) total_weight += p.second;
            uniform_int_distribution<> choice_dist(0, total_weight - 1);
            int r = choice_dist(rng);
            int cumsum = 0;
            for (auto& p : choices) {
                cumsum += p.second;
                if (r < cumsum) { chosen_op = p.first; break; }
            }
        }

        // ---- Emit the chosen instruction ------------------------------------
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
        if (chosen_op == CONST_OP) {
            int value = const_dist(rng);
            program.push_back({CONST_OP, value, -1});
            available.push_back(i);
            num_const++;
<<<<<<< HEAD
        }
        else if (chosen_op == INV_OP) {
            uniform_int_distribution<> pick(0, available.size() - 1);
            int idx = pick(rng);
            int x = available[idx];
            
            // VALIDATE: x must be < i (must exist before this instruction)
=======

        } else if (chosen_op == INV_OP) {
            // Pick one random operand from available, use it,
            // and produce a new result (add i to available).
            uniform_int_distribution<> pick(0, available.size() - 1);
            int idx = pick(rng);
            int x = available[idx];

            // operand must come before use.
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
            if (x >= i) {
                cerr << "ERROR at instruction " << i << ": INV references future node " << x << "\n";
                return 1;
            }
<<<<<<< HEAD
            
=======

>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
            available.erase(available.begin() + idx);
            program.push_back({INV_OP, x, -1});
            available.push_back(i);
            num_inv++;
<<<<<<< HEAD
        }
        else {  // MUL_OP
            shuffle(available.begin(), available.end(), rng);
            
            int x = available.back();
            available.pop_back();
            
            int y = available.back();
            available.pop_back();
            
            // VALIDATE: both x and y must be < i
            if (x >= i || y >= i) {
                cerr << "ERROR at instruction " << i << ": MUL references future nodes " 
                     << x << ", " << y << "\n";
                return 1;
            }
            
=======

        } else {  // MUL_OP
            // Pick two operands (without replacement) by shuffling and popping
            // from the back. Both are used, one new result is produced.
            // Shuffle first so we don't always take the two most recent results.
            shuffle(available.begin(), available.end(), rng);
            int x = available.back(); available.pop_back();
            int y = available.back(); available.pop_back();

            if (x >= i || y >= i) {
                cerr << "ERROR at instruction " << i << ": MUL references future nodes "
                     << x << ", " << y << "\n";
                return 1;
            }

>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
            program.push_back({MUL_OP, x, y});
            available.push_back(i);
            num_mul++;
        }
    }

<<<<<<< HEAD
    // VALIDATION: Check all references are valid
=======
    // ---- Final validation ---------------------------------------------------
    // Check that every operand reference is strictly forward (a < instruction index).
    // If structure_generator has a bug this catches it before it corrupts the CSV.
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
    cout << "Validating structure...\n";
    bool valid = true;
    for (size_t i = 0; i < program.size(); i++) {
        if (program[i].type == INV_OP) {
            if (program[i].a < 0 || program[i].a >= (int)i) {
<<<<<<< HEAD
                cerr << "ERROR: Instruction " << i << " (INV) references invalid node " 
                     << program[i].a << "\n";
                valid = false;
            }
        }
        else if (program[i].type == MUL_OP) {
            if (program[i].a < 0 || program[i].a >= (int)i) {
                cerr << "ERROR: Instruction " << i << " (MUL) first operand " 
=======
                cerr << "ERROR: Instruction " << i << " (INV) references invalid node "
                     << program[i].a << "\n";
                valid = false;
            }
        } else if (program[i].type == MUL_OP) {
            if (program[i].a < 0 || program[i].a >= (int)i) {
                cerr << "ERROR: Instruction " << i << " (MUL) first operand "
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
                     << program[i].a << " invalid\n";
                valid = false;
            }
            if (program[i].b < 0 || program[i].b >= (int)i) {
<<<<<<< HEAD
                cerr << "ERROR: Instruction " << i << " (MUL) second operand " 
=======
                cerr << "ERROR: Instruction " << i << " (MUL) second operand "
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033
                     << program[i].b << " invalid\n";
                valid = false;
            }
        }
    }
<<<<<<< HEAD
    
    if (!valid) {
        cerr << "Structure validation FAILED!\n";
        return 1;
    }
    cout << "✓ All references are valid!\n\n";

    // Write structure.txt
    ofstream out("structure.txt");
    out << "L " << L << "\n";

    for (int i = 0; i < L; i++) {
        if (program[i].type == CONST_OP)
            out << "C " << program[i].a << "\n";
        else if (program[i].type == INV_OP)
            out << "I " << program[i].a << "\n";
        else
            out << "M " << program[i].a << " " << program[i].b << "\n";
    }

    out.close();

    // Final statistics
    cout << "Tree structure generated and saved to structure.txt\n\n";
    cout << "Actual distribution:\n";
    cout << "  Constants      : " << num_const << " (" << (100.0 * num_const / L) << "%)\n";
    cout << "  Inversions     : " << num_inv   << " (" << (100.0 * num_inv / L) << "%)\n";
    cout << "  Multiplications: " << num_mul   << " (" << (100.0 * num_mul / L) << "%)\n";
    cout << "  Total          : " << (num_const + num_inv + num_mul) << "\n\n";

    cout << "Tree validation:\n";
    cout << "  Formula: C - M should equal 1 for valid tree\n";
    cout << "  Result : " << num_const << " - " << num_mul << " = " 
         << (num_const - num_mul) << "\n";
    
    if (num_const - num_mul == 1) {
        cout << "  ✓ VALID BINARY TREE!\n";
    } else {
        cout << "  ✗ NOT a valid tree\n";
    }
    
    cout << "\n  Final available nodes: " << available.size() 
         << " (should be 1 for single rooted tree)\n";
=======
    if (!valid) { cerr << "Structure validation FAILED!\n"; return 1; }
    cout << "All references are valid!\n\n";

    // ---- Write structure.txt ------------------------------------------------
    // Always written to the current working directory (relative path).
    // run_all.sh calls structure_generator from the project root, so all three
    // programs (structure_generator, querycds, tree_benchmark.g) read/write
    // the same structure.txt in the same directory.
    ofstream out("structure.txt");
    out << "L " << L << "\n";
    for (int i = 0; i < L; i++) {
        if      (program[i].type == CONST_OP) out << "C " << program[i].a << "\n";
        else if (program[i].type == INV_OP)   out << "I " << program[i].a << "\n";
        else                                  out << "M " << program[i].a << " " << program[i].b << "\n";
    }
    out.close();

    cout << "Tree structure generated and saved to structure.txt\n\n";
    cout << "Actual distribution:\n";
    cout << "  Constants      : " << num_const << " (" << (100.0 * num_const / L) << "%)\n";
    cout << "  Inversions     : " << num_inv   << " (" << (100.0 * num_inv   / L) << "%)\n";
    cout << "  Multiplications: " << num_mul   << " (" << (100.0 * num_mul   / L) << "%)\n";
    cout << "  Total          : " << (num_const + num_inv + num_mul) << "\n\n";
    cout << "Tree validation:\n";
    cout << "  C - M = " << num_const << " - " << num_mul << " = " << (num_const - num_mul)
         << (num_const - num_mul == 1 ? "  (VALID)\n" : "  (INVALID)\n");
    cout << "  Final available nodes: " << available.size()
         << " (should be 1 for a single-rooted tree)\n";
>>>>>>> 6af4dcabb7a806ea2dd54a14a417d2455c9e6033

    return 0;
}
