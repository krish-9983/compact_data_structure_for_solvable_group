// =============================================================================
// preparedsa.cpp  —  Stage 2: C++ preprocessing
// =============================================================================
//
// PURPOSE
//  It reads cayley_parseable.txt (written by extractGroupInfo.g), decides which
//   construction case applies (Theorem 3 alone, or Theorem 4 nested inside
//   Theorem 3). Based on that it builds all the auxiliary lookup tables, and writes them to
//   precomputed.bin (binary) and precomputed.txt (human-readable text).
//

// COMPILE
//   g++ -O3 -o preparedsa preparedsa.cpp
//
// RUN
//   ./preparedsa

//It reads cayley_parseable.txt and group_order.txt from current working directory
//
// CASE SELECTION for -> checkSub3
//   Case 1 (Mode 3): some Gi satisfies sqrt(n)/2 <= |Gi| <= sqrt(n)
//                    -> Theorem 3 alone, 11 array lookups per multiply
//   Case 2 (Mode 5): no such Gi exists then
//                    -> Theorem 4 on (Gi, Gi+1) + Theorem 3 on (G, Gi)
//                    -> 29 array lookups per multiply
// =============================================================================

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<cctype>
#include<math.h>
#include<cmath>
#include<climits>
#include<unordered_map>
#include<algorithm>
#include<set>
#include<chrono>
using namespace std;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

// This holds one composition series member as loaded from cayley_parseable.txt.
// table[0] is the element list in global G-indices 
struct Subgroup {
    int order = 0;
    vector<vector<int>> table;  // table[i][j] = global G-index of H[i]*H[j]
};

// CosetReps stores one left representative per left coset and one right
// representative per right coset. It store in local indices. 
struct CosetReps {
    vector<int> left, right;
};

// QuotientData is only needed in Case 2 to verify that Gi/Gi+1 is cyclic
struct QuotientData {
    vector<vector<int>> table;   // Cayley table of Gi/Gi+1 (coset indices)
    vector<int> repsL;           // left coset representatives of Gi+1 in Gi
};


// ---------------------------------------------------------------------------
// Binary/text I/O helpers
// ---------------------------------------------------------------------------

// trim_local:It strips leading and trailing whitespace from a string.
static string trim_local(const string &s) {
    int a = -1, b = -1;
    for (int i = 0; i < (int)s.size(); ++i)
        if (!isspace((unsigned char)s[i])) { a = i; break; }
    for (int i = (int)s.size() - 1; i >= 0; --i)
        if (!isspace((unsigned char)s[i])) { b = i; break; }
    if (a == -1 || b == -1) return string();
    return s.substr(a, b - a + 1);
}

//write_bin_block:it writes a 1-D int vector as a tagged binary block.
// querycds reads these blocks one by one and stores them in an unordered_map
// keyed by tag name. cols=0 signals "this is a 1-D array, not a 2-D matrix".
void write_bin_block(ofstream &out, const string &tag, const vector<int> &v) {
    int tagLength = (int)tag.size();
    out.write((char*)&tagLength, sizeof(tagLength));
    out.write(tag.c_str(), tagLength);
    int rows = (int)v.size();
    int cols = 0;
    out.write((char*)&rows, sizeof(rows));
    out.write((char*)&cols, sizeof(cols));
    if (rows > 0)
        out.write((const char*)v.data(), rows * sizeof(int));
}

// write_bin_mat: writes a 2-D int matrix as a tagged binary block.
// querycds reconstructs the matrix as block.data[r * cols + c].
void write_bin_mat(ofstream &out, const string &tag, const vector<vector<int>> &M) {
    int taglen = tag.size();
    out.write((char*)&taglen, sizeof(taglen));
    out.write(tag.data(), taglen);
    int rows = M.size();
    int cols = rows ? M[0].size() : 0;
    out.write((char*)&rows, sizeof(rows));
    out.write((char*)&cols, sizeof(cols));
    for (int r = 0; r < rows; r++)
        if (cols > 0)
            out.write((char*)M[r].data(), cols * sizeof(int));
}

// Text versions of the above — written to precomputed.txt for debugging.
// Format: "<tag> <rows> <cols>\n<values...>\n\n"
void write_text_block(ostream &out, const string &tag, const vector<int> &v) {
    out << tag << " " << v.size() << " " << 0 << "\n";
    for (int i = 0; i < (int)v.size(); i++) {
        if (i) out << " ";
        out << v[i];
    }
    out << "\n\n";
}

void write_text_mat(ostream &out, const string &tag, const vector<vector<int>> &M) {
    int rows = M.size();
    int cols = rows ? M[0].size() : 0;
    out << tag << " " << rows << " " << cols << "\n";
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (c) out << " ";
            out << M[r][c];
        }
        out << "\n";
    }
    out << "\n";
}

// ---------------------------------------------------------------------------
// checkSub3: decide which construction case applies
// Returns:
//   {i, -1}  -> Case 1: use groups[i] as H for Theorem 3
//   {-1, i}  -> Case 2: groups[i] is the "largest above sqrt(n)" Gi
// ---------------------------------------------------------------------------
pair<int,int> checkSub3(vector<Subgroup>& groups) {
    double rootG     = sqrt((double)groups[0].order);
    int upperBound   = (int)floor(rootG);       // floor(sqrt(n))
    int lowerBound   = (int)ceil(rootG / 2.0);  // ceil(sqrt(n)/2)
    int largestAbove = -1;  // tracks the last index still above sqrt(n)

    for (int i = 1; i < (int)groups.size(); i++) {
        int ord = groups[i].order;

        // This subgroup lands exactly in the sweet spot — Case 1.
        if (ord <= upperBound && ord >= lowerBound)
            return {i, -1};

        if (ord > upperBound)
            largestAbove = i;
        else {
            // If largestAbove is still -1 here, the series started below
            // sqrt(n) immediately (tiny group) — return this index for Case 2.
            if (largestAbove == -1)
                return {-1, i};
        }
    }
    // Fell through the whole series — return the last above-sqrt(n) index.
    return {-1, largestAbove};
}

// ---------------------------------------------------------------------------
// coset_representative_cal: compute left and right coset representatives

// all indices here are LOCAL to G (i.e. row indices into G.table)
// ---------------------------------------------------------------------------
CosetReps coset_representative_cal(Subgroup& G, Subgroup& H, int n) {
    const vector<vector<int>>& tableG = G.table;
    vector<int> elementsH = H.table[0];  // global IDs of H's elements

    vector<int> seenL(G.order, 0), seenR(G.order, 0);
    CosetReps reps;

    // Build a map from global element ID -> local G-index (row in G.table).
    // G.table[0][i] = global ID of G's i-th element.
    vector<int> global_to_local(n, -1);
    for (int i = 0; i < G.order; i++)
        global_to_local[G.table[0][i]] = i;

    // Left coset reps: find an unseen element, add it as a rep,
    // then mark every g*h as seen (covers the whole left coset g*H).
    for (int g = 0; g < G.order; g++) {
        if (seenL[g]) continue;
        reps.left.push_back(g);
        for (int h : elementsH) {
            int h_local     = global_to_local[h];
            int prod_global = tableG[g][h_local];
            int prod_local  = global_to_local[prod_global];
            seenL[prod_local] = 1;
        }
    }

    // Right coset reps: same idea but multiply on the left (h*g covers H*g).
    for (int g = 0; g < G.order; g++) {
        if (seenR[g]) continue;
        reps.right.push_back(g);
        for (int h : elementsH) {
            int h_local     = global_to_local[h];
            int prod_global = tableG[h_local][g];
            int prod_local  = global_to_local[prod_global];
            seenR[prod_local] = 1;
        }
    }
    return reps;
}

// ---------------------------------------------------------------------------
// quotient_group_table: build the Cayley table of Gi/Gi+1
//
// Used only in Case 2 to check whether Gi/Gi+1 is cyclic before applying
// Theorem 4.
// ---------------------------------------------------------------------------
QuotientData quotient_group_table(vector<Subgroup>& series, int i) {
    Subgroup& Gi   = series[i-1];
    Subgroup& Gi_1 = series[i];
    int n = series[0].order;
    const vector<vector<int>>& tableG = Gi.table;

    CosetReps reps = coset_representative_cal(Gi, Gi_1, n);
    const vector<int>& R = reps.left;
    int k = R.size();

    // Map each element of Gi to its coset index.
    vector<int> elementToCoset(Gi.order, -1);
    vector<int> global_to_local_Gi(n, -1);
    for (int j = 0; j < Gi.order; j++)
        global_to_local_Gi[Gi.table[0][j]] = j;

    for (int idx = 0; idx < k; idx++) {
        int r_local = R[idx];
        for (int h_global : Gi_1.table[0]) {
            int h_local  = global_to_local_Gi[h_global];
            int x_global = tableG[r_local][h_local];
            int x_local  = global_to_local_Gi[x_global];
            elementToCoset[x_local] = idx;
        }
    }

    // Quotient table: Q[a][b] = coset index of R[a] * R[b].
    vector<vector<int>> Q(k, vector<int>(k));
    for (int a = 0; a < k; a++)
        for (int b = 0; b < k; b++) {
            int prod_global = tableG[R[a]][R[b]];
            int prod_local  = global_to_local_Gi[prod_global];
            Q[a][b] = elementToCoset[prod_local];
        }

    QuotientData QD;
    QD.repsL = reps.left;
    QD.table = Q;
    return QD;
}

// ---------------------------------------------------------------------------
// is_cyclic: check if a group given by its Cayley table is cyclic
// A group is cyclic iff some element generates the whole group. 
// ---------------------------------------------------------------------------
bool is_cyclic(const vector<vector<int>>& Q) {
    int k = Q.size();
    for (int g = 1; g < k; g++) {
        vector<bool> seen(k, false);
        int cur = 0;
        for (int step = 0; step < k; step++) {
            if (seen[cur]) break;
            seen[cur] = true;
            cur = Q[cur][g];
        }
        if (count(seen.begin(), seen.end(), true) == k)
            return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// build_right_decomposition / build_left_decomposition
//
// For every g in G, compute the unique decomposition:
//   right: g = sR[g] * cR[g]  where cR[g] is a right coset rep, sR[g] in H
//   left:  g = cL[g] * sL[g]  where cL[g] is a left coset rep, sL[g] in H
//
// cL/cR store the coset-rep INDEX, sL/sR store the H-element
// INDEX (not element). querycds uses these indices to look up subsequent
// tables during multiplication.
// ---------------------------------------------------------------------------
bool build_right_decomposition(
    int n,
    const vector<int>& H_elems,
    const vector<int>& repsR,
    const vector<vector<int>>& G_table,
    vector<int>& cR, vector<int>& sR
) {
    int m = H_elems.size();
    int k = repsR.size();
    for (int repIdx = 0; repIdx < k; repIdx++) {
        int r = repsR[repIdx];
        for (int hIdx = 0; hIdx < m; hIdx++) {
            int hglob = H_elems[hIdx];
            int g = G_table[hglob][r];  // g = h * r  (right decomp: g = sR * cR)
            if (cR[g] == -1) {
                cR[g] = repIdx;
                sR[g] = hIdx;
            }
        }
    }
    return true;
}

bool build_left_decomposition(
    int n,
    const vector<int>& H_elems,
    const vector<int>& repsL,
    const vector<vector<int>>& G_table,
    vector<int>& cL, vector<int>& sL
) {
    int m = H_elems.size();
    int k = repsL.size();
    for (int repIdx = 0; repIdx < k; repIdx++) {
        int l = repsL[repIdx];
        for (int hIdx = 0; hIdx < m; hIdx++) {
            int hglob = H_elems[hIdx];
            int g = G_table[l][hglob];  // g = l * h  (left decomp: g = cL * sL)
            if (cL[g] == -1) {
                cL[g] = repIdx;
                sL[g] = hIdx;
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// build_flip_RH
//
// Computes flipH and flipR tables.
// For each (l, h) in L x H we need to express l * h as h' * r where h' in H and r in R. 
// ---------------------------------------------------------------------------
void build_flip_RH(
    const vector<int>& repsL,
    const vector<int>& repsR,
    const vector<int>& H_elems,
    const vector<int>& cR, const vector<int>& sR,
    const vector<vector<int>>& G_table,
    vector<vector<int>>& flipH, vector<vector<int>>& flipR
) {
    int L = repsL.size();
    int m = H_elems.size();
    flipH.assign(L, vector<int>(m, -1));
    flipR.assign(L, vector<int>(m, -1));

    for (int lIdx = 0; lIdx < L; lIdx++) {
        int l = repsL[lIdx];
        for (int hIdx = 0; hIdx < m; hIdx++) {
            int h = H_elems[hIdx];
            int g = G_table[l][h];         // l * h
            flipH[lIdx][hIdx] = sR[g];     // H-component of the right decomp of l*h
            flipR[lIdx][hIdx] = cR[g];     // R-component of the right decomp of l*h
        }
    }
}

// ---------------------------------------------------------------------------
// build_cross_RR
//
// Computes crossH and crossR tables.
// For each (r1, r2) in R x R: r1 * r2 = CrossH[r1][r2] * CrossR[r1][r2]
// where CrossH is in H and CrossR is in R.
// ---------------------------------------------------------------------------
void build_cross_RR(
    const vector<int>& repsR,
    const vector<int>& H_elems,
    const vector<int>& cR, const vector<int>& sR,
    const vector<vector<int>>& G_table,
    vector<vector<int>>& CrossH, vector<vector<int>>& CrossR
) {
    int R = repsR.size();
    CrossH.assign(R, vector<int>(R, -1));
    CrossR.assign(R, vector<int>(R, -1));

    for (int i = 0; i < R; ++i) {
        int r1 = repsR[i];
        for (int j = 0; j < R; ++j) {
            int r2 = repsR[j];
            int g         = G_table[r1][r2];  // r1 * r2
            CrossH[i][j] = sR[g];             // H-component
            CrossR[i][j] = cR[g];             // R-component
        }
    }
}

// ---------------------------------------------------------------------------
// build_fuse
//
// Computes fuse table
// fuse[hIdx][rIdx] = global G-index of H_elems[hIdx] * repsR[rIdx].
// ---------------------------------------------------------------------------
void build_fuse(
    const vector<int>& H_elems,
    const vector<int>& repsR,
    const vector<vector<int>>& G_table,
    vector<vector<int>>& fuse
) {
    int m = H_elems.size();
    int R = repsR.size();
    fuse.assign(m, vector<int>(R, -1));
    for (int hIdx = 0; hIdx < m; hIdx++) {
        int h = H_elems[hIdx];
        for (int rIdx = 0; rIdx < R; rIdx++)
            fuse[hIdx][rIdx] = G_table[h][repsR[rIdx]];
    }
}

// ---------------------------------------------------------------------------
// build_HH
//
// Builds HH: the Cayley table of H in LOCAL H-indices.
// HH[i][j] = local H-index of H_elems[i] * H_elems[j].
//
// This is stored in local H-indices
//
// global_to_Hidx maps a global G element ID to its local H-index, pre-built
// by the caller from H_elems.
// ---------------------------------------------------------------------------
void build_HH(
    const vector<int>& H_elems,
    const vector<vector<int>>& G_table,
    const vector<int>& global_to_Hidx,
    vector<vector<int>>& HH
) {
    int m = H_elems.size();
    HH.assign(m, vector<int>(m, -1));
    for (int i = 0; i < m; i++) {
        int hgi = H_elems[i];
        for (int j = 0; j < m; j++) {
            int hgj  = H_elems[j];
            int prodg = G_table[hgi][hgj];
            HH[i][j] = global_to_Hidx[prodg];  // local H-index of the product
        }
    }
}

// ---------------------------------------------------------------------------
// find_identity / build_inverse_table
//
// find_identity: scan the table for the element e such that e*g = g*e = g
// for all g.
// build_inverse_table: for each element i, find j such that i*j = j*i = e.
// Stored as global G-indices.
// ---------------------------------------------------------------------------
int find_identity(const vector<vector<int>>& G_table) {
    int n = G_table.size();
    for (int i = 0; i < n; i++) {
        bool leftId = true, rightId = true;
        for (int j = 0; j < n; j++) {
            if (G_table[i][j] != j) leftId  = false;
            if (G_table[j][i] != j) rightId = false;
            if (!leftId && !rightId) break;
        }
        if (leftId && rightId) return i;
    }
    return -1;
}

vector<int> build_inverse_table(const vector<vector<int>>& G_table, int id) {
    int n = G_table.size();
    vector<int> inv(n, -1);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (G_table[i][j] == id && G_table[j][i] == id) {
                inv[i] = j;
                break;
            }
    return inv;
}

// ---------------------------------------------------------------------------
// Theorem 4 related functions (only used in Case 2)
// ---------------------------------------------------------------------------

// build_g0_powers: precompute g0^0, g0^1, ..., g0^maxExp using the full G table.
// g0 is the generator of the cyclic quotient Gi/Gi+1. We need powers up to
// 2k-2 (where k = [Gi:Gi+1]).
vector<int> build_g0_powers(const vector<vector<int>>& G_table, int g0, int id, int maxExp) {
    vector<int> temp(maxExp+1, -1);
    temp[0] = id;
    for (int i = 1; i <= maxExp; ++i)
        temp[i] = G_table[temp[i-1]][g0];
    return temp;
}

// build_e_sR_sL: for every element g in Gi, compute:
//   e[g]  = exponent i such that g = g0^i * n  for some n in N
//   sR[g] = N-index of n  (the "right N-component")
//   sL[g] = N-index of n' where g = n' * g0^i  (the "left N-component")
//
// We fill by iterating over all (g0^i, n) pairs.
bool build_e_sR_sL(
    const vector<vector<int>>& G_table,
    const vector<int>& N_elems, int k,
    const vector<int>& inv,
    vector<int>& e, vector<int>& sR, vector<int>& sL,
    vector<int>& gid_to_nidx, vector<int>& g0pow
) {
    int n = G_table.size();
    e.assign(n, -1); sR.assign(n, -1); sL.assign(n, -1);

    for (int i = 0; i < k; ++i) {
        int g0i = g0pow[i];
        for (int ni = 0; ni < (int)N_elems.size(); ni++) {
            int n_gid  = N_elems[ni];
            int g      = G_table[g0i][n_gid];   // g = g0^i * n
            e[g]       = i;
            sR[g]      = ni;
            // left N-component: n' = g * g0^{-i}
            int inv_g0i    = inv[g0i];
            int nprime_gid = G_table[g][inv_g0i];
            sL[g]          = gid_to_nidx[nprime_gid];
        }
    }
    return true;
}

// build_Flip: for each (n, i), compute n * g0^i gave in N * g0^j form.
// Specifically, Flip[n_idx][i] = N-index of the N-part when 
//   n * g0^i = g0^i * Flip(n, i)
vector<vector<int>> build_Flip(
    const vector<vector<int>>& G_table,
    const vector<int>& N_elems,
    const vector<int>& gid_to_nidx,
    const vector<int>& inv,
    const vector<int>& g0pow_k
) {
    int m = N_elems.size();
    int k = g0pow_k.size();
    vector<vector<int>> Flip(m, vector<int>(k, -1));

    for (int ni = 0; ni < m; ni++) {
        int n_gid = N_elems[ni];
        for (int i = 0; i < k; i++) {
            int g0i      = g0pow_k[i];
            int lhs      = G_table[n_gid][g0i];       // n * g0^i
            int flip_gid = G_table[inv[g0i]][lhs];    // g0^{-i} * (n * g0^i)
            Flip[ni][i]  = gid_to_nidx[flip_gid];
        }
    }
    return Flip;
}

// build_rede_redN: for each exponent sum ell = alpha+beta (ranging 0..2k-2),
// reduce it mod k to get rede[ell] (the canonical exponent in [0,k)), and
// compute the corresponding N-element: redN[ell] = N-index of g0^ell * g0^{-rede}.
void build_rede_redN(
    const vector<vector<int>>& G_table,
    const vector<int>& g0pow_full,
    const vector<int>& g0pow_k,
    const vector<int>& inv,
    const vector<int>& gid_to_nidx,
    vector<int>& rede, vector<int>& redN_idx
) {
    int maxEll = g0pow_full.size() - 1;
    int k      = g0pow_k.size();
    rede.assign(maxEll+1, -1);
    redN_idx.assign(maxEll+1, -1);

    for (int ell = 0; ell <= maxEll; ell++) {
        int rede_i     = ell % k;
        rede[ell]      = rede_i;
        int g_ell_gid  = g0pow_full[ell];
        int g_rede_gid = g0pow_k[rede_i];
        // redN[ell] is the N-element that "corrects" g0^ell back to g0^rede:
        //   g0^ell = g0^rede * redN[ell]
        int redN_gid   = G_table[inv[g_rede_gid]][g_ell_gid];
        redN_idx[ell]  = gid_to_nidx[redN_gid];
    }
}

// build_Fuse (Theorem 4): Fuse[gamma][ni] = global G-index of g0^gamma * N_elems[ni].
vector<vector<int>> build_Fuse(
    const vector<vector<int>>& G_table,
    const vector<int>& g0pow_k,
    const vector<int>& N_elems
) {
    int k = g0pow_k.size();
    int m = N_elems.size();
    vector<vector<int>> Fuse(k, vector<int>(m, -1));
    for (int i = 0; i < k; i++) {
        int g0i = g0pow_k[i];
        for (int ni = 0; ni < m; ni++)
            Fuse[i][ni] = G_table[g0i][N_elems[ni]];
    }
    return Fuse;
}

// ---------------------------------------------------------------------------
// write_theorem3_data / write_theorem4_data
//
// Bundle all the Theorem 3 or Theorem 4 arrays into precomputed.bin and
// precomputed.txt. The order matters for querycds: it reads blocks by tag
// name from an unordered_map.
//
// write_theorem4_data is always called BEFORE write_theorem3_data
// in Case 2. querycds detects Case 2 by checking for the presence of "e_arr"
// in the block map.
//
// writeHH controls whether HH_mat is written.
// Case 1: writeHH=true  — multiply_theorem3 reads HH_mat twice per query.
// Case 2: writeHH=false — HH_mat is skipped to keep O(n) space.
//   |Gi|^2 > n in Case 2, so writing it would break the space guarantee.
//   multiply_combined never reads HH_mat; T4 handles those two steps.
// ---------------------------------------------------------------------------
void write_theorem3_data(
    ofstream& bout, ofstream& tout,
    const vector<int>& subgroupH, const CosetReps& reps,
    const vector<int>& cL, const vector<int>& sL,
    const vector<int>& cR, const vector<int>& sR,
    const vector<vector<int>>& HH,
    const vector<vector<int>>& flipH, const vector<vector<int>>& flipR,
    const vector<vector<int>>& crossH, const vector<vector<int>>& crossR,
    const vector<vector<int>>& fuse,
    const vector<int>& inverse,
    bool writeHH = true
) {
    write_bin_block(bout, "subgroupH", subgroupH);
    write_bin_block(bout, "repsL",     reps.left);
    write_bin_block(bout, "repsR",     reps.right);
    write_bin_block(bout, "cL",        cL);
    write_bin_block(bout, "sL",        sL);
    write_bin_block(bout, "cR",        cR);
    write_bin_block(bout, "sR",        sR);
    if (writeHH)
        write_bin_mat(bout, "HH_mat",  HH);
    write_bin_mat  (bout, "flipH_mat", flipH);
    write_bin_mat  (bout, "flipR_mat", flipR);
    write_bin_mat  (bout, "crossH_mat",crossH);
    write_bin_mat  (bout, "crossR_mat",crossR);
    write_bin_mat  (bout, "fuse_mat",  fuse);
    write_bin_block(bout, "inverse",   inverse);

    tout << "# theorem3 data\n";
    write_text_block(tout, "subgroupH", subgroupH);
    write_text_block(tout, "repsL",     reps.left);
    write_text_block(tout, "repsR",     reps.right);
    write_text_block(tout, "cL",  cL);
    write_text_block(tout, "sL",  sL);
    write_text_block(tout, "cR",  cR);
    write_text_block(tout, "sR",  sR);
    if (writeHH)
        write_text_mat(tout, "HH_mat",  HH);
    write_text_mat  (tout, "flipH_mat", flipH);
    write_text_mat  (tout, "flipR_mat", flipR);
    write_text_mat  (tout, "crossH_mat",crossH);
    write_text_mat  (tout, "crossR_mat",crossR);
    write_text_mat  (tout, "fuse_mat",  fuse);
}

void write_theorem4_data(
    ofstream& bout, ofstream& tout,
    const vector<int>& N_elems,
    const vector<int>& g0pow_k, const vector<int>& g0pow_full,
    const vector<int>& e, const vector<int>& sR, const vector<int>& sL,
    const vector<vector<int>>& NN, const vector<vector<int>>& flip,
    const vector<int>& rede, const vector<int>& redN_idx,
    const vector<vector<int>>& fuse,
    const vector<int>& inverse
) {
    write_bin_block(bout, "N_elems",    N_elems);
    write_bin_block(bout, "g0pow_k",    g0pow_k);
    write_bin_block(bout, "g0pow_full", g0pow_full);
    write_bin_block(bout, "e_arr",      e);      // PRESENCE of "e_arr" signals Mode 5 to querycds
    write_bin_block(bout, "sR_arr",     sR);
    write_bin_block(bout, "sL_arr",     sL);
    write_bin_mat  (bout, "NN_mat",     NN);
    write_bin_mat  (bout, "Flip_mat",   flip);
    write_bin_block(bout, "rede_arr",   rede);
    write_bin_block(bout, "redN_idx",   redN_idx);
    write_bin_mat  (bout, "Fuse_mat",   fuse);
    write_bin_block(bout, "inverse",    inverse);

    tout << "# theorem4 data\n";
    write_text_block(tout, "N_elems",    N_elems);
    write_text_block(tout, "g0pow_k",    g0pow_k);
    write_text_block(tout, "g0pow_full", g0pow_full);
    write_text_block(tout, "e_arr",      e);
    write_text_block(tout, "sR_arr",     sR);
    write_text_block(tout, "sL_arr",     sL);
    write_text_mat  (tout, "NN_mat",     NN);
    write_text_mat  (tout, "Flip_mat",   flip);
    write_text_block(tout, "rede_arr",   rede);
    write_text_block(tout, "redN_idx",   redN_idx);
    write_text_mat  (tout, "Fuse_mat",   fuse);
}

// ===========================================================================
// main
// ===========================================================================
int main() {
    auto start = std::chrono::high_resolution_clock::now();

    // ---- Parse cayley_parseable.txt ----------------------------------------
    // groups[0] = full group G; groups[1..] = composition series subgroups
    string path = "cayley_parseable.txt";
    ifstream fin(path);
    if (!fin.is_open()) {
        cerr << "Error: cannot open file: " << path << "\n";
        return 1;
    }

    vector<Subgroup> groups;
    string raw;

    while (getline(fin, raw)) {
        string ln = trim_local(raw);
        if (ln.empty()) continue;

        istringstream iss(ln);
        int n;
        if (!(iss >> n) || n <= 0) continue;

        Subgroup g;
        g.order = n;
        g.table.assign(n, vector<int>(n));
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                iss >> g.table[r][c];
        groups.push_back(move(g));
    }
    fin.close();

    if (groups.empty()) {
        cerr << "Error: no groups parsed from " << path << "\n";
        return 1;
    }
    cout << "Parsed " << groups.size() << " subgroups.\n";
    cout << "Full group order: " << groups[0].order << "\n";

    // ---- Decide which case applies -----------------------------------------
    pair<int,int> isExist = checkSub3(groups);
    int orderG = groups[0].order;

    ofstream bout("precomputed.bin", ios::binary);
    if (!bout.is_open()) { cerr << "Cannot open precomputed.bin\n"; return 1; }
    ofstream tout("precomputed.txt");
    if (!tout.is_open()) { cerr << "Cannot open precomputed.txt\n"; return 1; }

    // Inverse table is needed by both cases (it goes into the binary for querycds).
    const Subgroup& G = groups[0];
    int id = find_identity(G.table);
    vector<int> inverse = build_inverse_table(G.table, id);

    // =========================================================================
    // CASE 1 — Theorem 3 alone (Mode 3)
    //
    // There exists a Gi with |Gi| in [sqrt(n)/2, sqrt(n)].
    // Theorem 3 applied to G with H = Gi gives an O(n) structure:
    // querycds will run multiply_theorem3 (11 array lookups per query).
    // =========================================================================
    if (isExist.first != -1) {
        int subIdx = isExist.first;
        cout << "Case 1: Theorem 3 with subgroup index " << subIdx
             << "  (|Gi| = " << groups[subIdx].order << ")\n";

        CosetReps reps = coset_representative_cal(groups[0], groups[subIdx], orderG);

        vector<int> subgroupH = groups[subIdx].table[0];  // global IDs of H's elements
        vector<int> cR(orderG, -1), sR(orderG, -1);
        vector<int> cL(orderG, -1), sL(orderG, -1);

        // Build a reverse map: global G-index -> local H-index
        vector<int> global_to_Hidx(orderG, -1);
        for (int i = 0; i < (int)subgroupH.size(); ++i)
            global_to_Hidx[subgroupH[i]] = i;

        vector<vector<int>> HH;
        build_HH(subgroupH, groups[0].table, global_to_Hidx, HH);
        build_right_decomposition(orderG, subgroupH, reps.right, groups[0].table, cR, sR);
        build_left_decomposition (orderG, subgroupH, reps.left,  groups[0].table, cL, sL);

        vector<vector<int>> flipH, flipR, crossH, crossR, fuse;
        build_flip_RH (reps.left, reps.right, subgroupH, cR, sR, groups[0].table, flipH, flipR);
        build_cross_RR(reps.right, subgroupH, cR, sR, groups[0].table, crossH, crossR);
        build_fuse    (subgroupH, reps.right, groups[0].table, fuse);

        write_theorem3_data(bout, tout,
            subgroupH, reps, cL, sL, cR, sR,
            HH, flipH, flipR, crossH, crossR, fuse, inverse);
            // writeHH defaults to true

    } else {
        // =====================================================================
        // CASE 2 — Theorem 4 nested inside Theorem 3 (Mode 5)
        //
        // No Gi lands in [sqrt(n)/2, sqrt(n)]. The composition series jumps:
        //   Gi  with |Gi| > sqrt(n)
        //   Gi+1 with |Gi+1| < sqrt(n)/2
        //
        // Step A: Apply Theorem 4 to the inner pair (Gi, Gi+1).
        // Step B: Apply Theorem 3 to the outer pair (G, Gi).
        // =====================================================================
        if (isExist.second == -1) {
            cout << "No suitable subgroup found.\n";
            return 0;
        }

        int sgroupInx = isExist.second;
        int nextInx   = sgroupInx + 1;

        int quotientIdx = (nextInx >= (int)groups.size()) ? sgroupInx : nextInx;
        QuotientData QD = quotient_group_table(groups, quotientIdx);

        if (is_cyclic(QD.table)) {
            const Subgroup& Nsub  = groups[nextInx];
            const Subgroup& Gisub = groups[sgroupInx];

            int orderGi = Gisub.order;
            int orderN  = Nsub.order;
            const vector<int>& N_elems = Nsub.table[0];  // global IDs of N's elements
            int k = orderGi / orderN;  // order of the cyclic quotient Gi/N

            // Build N membership and index maps
            vector<int> isInN(orderG, 0);
            vector<int> gid_to_nid(orderG, -1);
            for (int i = 0; i < (int)N_elems.size(); i++) {
                isInN[N_elems[i]]      = 1;
                gid_to_nid[N_elems[i]] = i;
            }

            // Find a generator g0 of the cyclic group Gi/N.
            int g0 = -1;
            for (int i = 0; i < orderGi; i++) {
                int cur = Gisub.table[0][i];
                if (isInN[cur]) continue;

                int pow = cur;
                bool earlyInN = false;
                for (int step = 1; step <= k - 2; step++) {
                    pow = G.table[pow][cur];
                    if (isInN[pow]) { earlyInN = true; break; }
                }
                if (earlyInN) continue;

                int kth = G.table[pow][cur];
                if (isInN[kth]) { g0 = cur; break; }
            }

            if (g0 == -1) {
                cerr << "Error: could not find generator g0 for Gi/N.\n";
                return 1;
            }

            // Precompute powers of g0 up to 2k-2 (needed for rede/redN reduction).
            int maxExp = 2*k - 2;
            vector<int> g0pow = build_g0_powers(G.table, g0, id, maxExp);
            vector<int> g0pow_k(g0pow.begin(), g0pow.begin() + k);

            // Step A: build all Theorem 4 tables for the inner pair (Gi, N).
            vector<int> e_arr, sR_arr, sL_arr;
            build_e_sR_sL(G.table, N_elems, k, inverse,
                          e_arr, sR_arr, sL_arr, gid_to_nid, g0pow_k);

            vector<vector<int>> NN;
            build_HH(N_elems, G.table, gid_to_nid, NN);

            vector<vector<int>> flip_mat =
                build_Flip(G.table, N_elems, gid_to_nid, inverse, g0pow_k);

            vector<int> rede, redN_idx;
            build_rede_redN(G.table, g0pow, g0pow_k, inverse, gid_to_nid, rede, redN_idx);

            vector<vector<int>> fuse_t4 = build_Fuse(G.table, g0pow_k, N_elems);

            // Write T4 data FIRST — querycds uses "e_arr" presence to detect Mode 5.
            write_theorem4_data(bout, tout,
                N_elems, g0pow_k, g0pow,
                e_arr, sR_arr, sL_arr,
                NN, flip_mat, rede, redN_idx, fuse_t4, inverse);

            cout << "Case 2, Step A: Theorem 4 on Gi (index " << sgroupInx
                 << ", |Gi|=" << orderGi << ") with N=Gi+1 (index "
                 << nextInx << ", |N|=" << orderN << ")\n";

            // Step B: build all Theorem 3 tables for the outer pair (G, Gi).
            // HH_mat is NOT built or written here: |Gi|^2 > n in Case 2,
            // which would exceed O(n) space. multiply_combined never reads
            // HH_mat in Mode 5 — T4 handles those two inner multiplications.
            CosetReps reps3 = coset_representative_cal(groups[0], groups[sgroupInx], orderG);

            vector<int> subgroupH3 = groups[sgroupInx].table[0];
            vector<int> cR3(orderG, -1), sR3(orderG, -1);
            vector<int> cL3(orderG, -1), sL3(orderG, -1);

            vector<int> global_to_Hidx3(orderG, -1);
            for (int j = 0; j < (int)subgroupH3.size(); ++j)
                global_to_Hidx3[subgroupH3[j]] = j;

            vector<vector<int>> HH3;  // left empty — not written (writeHH=false below)
            build_right_decomposition(orderG, subgroupH3, reps3.right, groups[0].table, cR3, sR3);
            build_left_decomposition (orderG, subgroupH3, reps3.left,  groups[0].table, cL3, sL3);

            vector<vector<int>> flipH3, flipR3, crossH3, crossR3, fuse3;
            build_flip_RH (reps3.left, reps3.right, subgroupH3, cR3, sR3, groups[0].table, flipH3, flipR3);
            build_cross_RR(reps3.right, subgroupH3, cR3, sR3, groups[0].table, crossH3, crossR3);
            build_fuse    (subgroupH3, reps3.right, groups[0].table, fuse3);

            write_theorem3_data(bout, tout,
                subgroupH3, reps3, cL3, sL3, cR3, sR3,
                HH3, flipH3, flipR3, crossH3, crossR3, fuse3, inverse,
                false);  // writeHH=false: skip HH_mat to maintain O(n) space

            cout << "Case 2, Step B: Theorem 3 on G with Gi (index " << sgroupInx << ")\n";

        } else {
            cerr << "Error: non-cyclic quotient for solvable group.\n";
            return 1;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    double ms = chrono::duration<double,milli>(end - start).count();
    cout << "preparedsa done in " << ms << " ms\n";
    return 0;
}