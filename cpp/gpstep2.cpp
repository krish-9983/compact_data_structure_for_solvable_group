#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include <algorithm>
#include <stdexcept>
#include<string>

struct Group {
    int n{};
    std::vector<std::vector<int>> mul; // mul[a][b] ∈ [0..n-1]
};

struct GroupChecker {
    static int find_identity(const Group& G) {
        for (int e = 0; e < G.n; ++e) {
            bool ok = true;
            for (int x = 0; x < G.n; ++x) {
                if (G.mul[e][x] != x || G.mul[x][e] != x) { ok = false; break; }
            }
            if (ok) return e;
        }
        return -1;
    }
    static bool has_inverses(const Group& G, int e) {
        for (int a = 0; a < G.n; ++a) {
            bool ok = false;
            for (int b = 0; b < G.n; ++b) {
                if (G.mul[a][b] == e && G.mul[b][a] == e) { ok = true; break; }
            }
            if (!ok) return false;
        }
        return true;
    }
    static bool is_group(const Group& G) {
        int e = find_identity(G);
        if (e < 0) return false;
        if (!has_inverses(G, e)) return false;
        // Lightweight associativity spot-check
        std::mt19937 rng(42);
        std::uniform_int_distribution<int> U(0, G.n - 1);
        for (int t = 0; t < std::min(1000, G.n * 10); ++t) {
            int a = U(rng), b = U(rng), c = U(rng);
            if (G.mul[G.mul[a][b]][c] != G.mul[a][G.mul[b][c]]) return false;
        }
        return true;
    }
};

struct Subgroup {
    std::vector<int> elems;        // elements of H
    std::vector<std::uint8_t> inH; // membership bitset size n
};

Subgroup make_subgroup(const Group& G, const std::vector<int>& elems) {
    Subgroup H;
    H.elems = elems;
    H.inH.assign(G.n, 0);
    for (int x : elems) {
        if (x < 0 || x >= G.n) throw std::runtime_error("Subgroup element out of range");
        H.inH[x] = 1;
    }
    return H;
}

// Left coset reps L: one per coset gH
std::vector<int> left_reps(const Group& G, const Subgroup& H) {
    std::vector<int> L;
    std::vector<std::uint8_t> seen(G.n, 0);
    for (int g = 0; g < G.n; ++g) {
        if (seen[g]) continue;
        L.push_back(g);
        for (int h : H.elems) seen[G.mul[g][h]] = 1;
    }
    return L;
}

// Right coset reps R: one per coset Hg
std::vector<int> right_reps(const Group& G, const Subgroup& H) {
    std::vector<int> R;
    std::vector<std::uint8_t> seen(G.n, 0);
    for (int g = 0; g < G.n; ++g) {
        if (seen[g]) continue;
        R.push_back(g);
        for (int h : H.elems) seen[G.mul[h][g]] = 1;
    }
    return R;
}

struct Splits {
    std::vector<int> cL, sL; // g = cL[g] * sL[g], cL ∈ L, sL ∈ H
    std::vector<int> cR, sR; // g = sR[g] * cR[g], sR ∈ H, cR ∈ R
};

// O(n·|L|·|H| + n·|R|·|H|)
Splits compute_splits(const Group& G, const Subgroup& H,
                      const std::vector<int>& L, const std::vector<int>& R) {
    Splits S;
    S.cL.assign(G.n, -1);
    S.sL.assign(G.n, -1);
    S.cR.assign(G.n, -1);
    S.sR.assign(G.n, -1);

    // Left split: find (l,h) with l∈L, h∈H, l*h = g
    for (int g = 0; g < G.n; ++g) {
        for (int l : L) {
            bool found = false;
            for (int h : H.elems) {
                if (G.mul[l][h] == g) {
                    S.cL[g] = l; S.sL[g] = h;
                    found = true; break;
                }
            }
            if (found) break;
        }
        if (S.cL[g] == -1) throw std::runtime_error("Left split failed for element " + std::to_string(g));
    }

    // Right split: find (h,r) with h∈H, r∈R, h*r = g
    for (int g = 0; g < G.n; ++g) {
        for (int r : R) {
            bool found = false;
            for (int h : H.elems) {
                if (G.mul[h][r] == g) {
                    S.cR[g] = r; S.sR[g] = h;
                    found = true; break;
                }
            }
            if (found) break;
        }
        if (S.cR[g] == -1) throw std::runtime_error("Right split failed for element " + std::to_string(g));
    }

    return S;
}

// Index maps: element -> position in L/R/H vectors (or -1 if not present)
struct IndexMaps {
    std::vector<int> idxL, idxR, idxH;
};

IndexMaps build_index_maps(const Group& G,
                           const std::vector<int>& L,
                           const std::vector<int>& R,
                           const Subgroup& H) {
    IndexMaps M;
    M.idxL.assign(G.n, -1);
    M.idxR.assign(G.n, -1);
    M.idxH.assign(G.n, -1);
    for (int i = 0; i < (int)L.size(); ++i) M.idxL[L[i]] = i;
    for (int i = 0; i < (int)R.size(); ++i) M.idxR[R[i]] = i;
    for (int i = 0; i < (int)H.elems.size(); ++i) M.idxH[H.elems[i]] = i;
    return M;
}

// Extension A tables
struct ExtATables {
    // Dimensions:
    // Flip*: |L| x |H|
    // Cross*: |R| x |R|
    // Fuse: |H| x |R|
    std::vector<std::vector<int>> FlipH, FlipR;
    std::vector<std::vector<int>> CrossH, CrossR;
    std::vector<std::vector<int>> Fuse;
};

// Build Flip/Cross/Fuse from splits
ExtATables build_extA_tables(const Group& G,
                             const Subgroup& H,
                             const std::vector<int>& L,
                             const std::vector<int>& R,
                             const Splits& S,
                             const IndexMaps& M) {
    ExtATables T;
    const int LH = (int)L.size(), RH = (int)R.size(), HH = (int)H.elems.size();

    T.FlipH.assign(LH, std::vector<int>(HH, -1));
    T.FlipR.assign(LH, std::vector<int>(HH, -1));
    T.CrossH.assign(RH, std::vector<int>(RH, -1));
    T.CrossR.assign(RH, std::vector<int>(RH, -1));
    T.Fuse.assign(HH, std::vector<int>(RH, -1));

    // Flip: for each l∈L, h∈H: x = l*h, then use right split of x
    for (int i = 0; i < LH; ++i) {
        int l = L[i];
        for (int j = 0; j < HH; ++j) {
            int h = H.elems[j];
            int x = G.mul[l][h];
            int hprime = S.sR[x]; // right split H-part
            int r = S.cR[x];      // right split R-rep
            T.FlipH[i][j] = hprime;
            T.FlipR[i][j] = r;
        }
    }

    // Cross: for each r1,r2∈R: x = r1*r2, then use right split of x
    for (int i = 0; i < RH; ++i) {
        int r1 = R[i];
        for (int j = 0; j < RH; ++j) {
            int r2 = R[j];
            int x = G.mul[r1][r2];
            T.CrossH[i][j] = S.sR[x]; // H-part
            T.CrossR[i][j] = S.cR[x]; // R-rep
        }
    }

    // Fuse: for each h∈H, r∈R: store h*r
    for (int i = 0; i < HH; ++i) {
        int h = H.elems[i];
        for (int j = 0; j < RH; ++j) {
            int r = R[j];
            T.Fuse[i][j] = G.mul[h][r];
        }
    }

    return T;
}

// Constant-time multiply via Extension A
int multiply_extA(const Group& G,
                  int g1, int g2,
                  const Splits& S,
                  const Subgroup& H,
                  const std::vector<int>& L,
                  const std::vector<int>& R,
                  const IndexMaps& M,
                  const ExtATables& T) {
    // Step 1: splits
    int l  = S.cL[g1];
    int s1 = S.sL[g1];
    int s2 = S.sR[g2];
    int r  = S.cR[g2];

    // Step 2: multiply in H
    int h1 = G.mul[s1][s2];

    // Step 3: Flip (l,h1) -> h2, r1
    int il = M.idxL[l];
    int ih1 = M.idxH[h1];
    if (il < 0 || ih1 < 0) throw std::runtime_error("Flip index error");
    int h2 = T.FlipH[il][ih1];
    int r1 = T.FlipR[il][ih1];

    // Step 4: Cross (r1,r) -> h3, r3
    int ir1 = M.idxR[r1];
    int ir  = M.idxR[r];
    if (ir1 < 0 || ir < 0) throw std::runtime_error("Cross index error");
    int h3 = T.CrossH[ir1][ir];
    int r3 = T.CrossR[ir1][ir];

    // Step 5: multiply in H
    int h4 = G.mul[h2][h3];

    // Step 6: Fuse (h4, r3)
    int ih4 = M.idxH[h4];
    int ir3 = M.idxR[r3];
    if (ih4 < 0 || ir3 < 0) throw std::runtime_error("Fuse index error");
    int g = T.Fuse[ih4][ir3];
    return g;
}

// Input format:
// n
// n x n Cayley table (0..n-1)
// h
// h integers: elements of H
// q
// q lines: g1 g2
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Group G;
    if (!(std::cin >> G.n)) return 0;
    G.mul.assign(G.n, std::vector<int>(G.n));
    for (int i = 0; i < G.n; ++i)
        for (int j = 0; j < G.n; ++j)
            std::cin >> G.mul[i][j];

    if (!GroupChecker::is_group(G)) {
        std::cout << "Invalid group Cayley table.\n";
        return 0;
    }

    int hsz; std::cin >> hsz;
    std::vector<int> he(hsz);
    for (int i = 0; i < hsz; ++i) std::cin >> he[i];
    Subgroup H = make_subgroup(G, he);

    // Build cosets and splits
    auto L = left_reps(G, H);
    auto R = right_reps(G, H);
    Splits S = compute_splits(G, H, L, R);
    IndexMaps M = build_index_maps(G, L, R, H);
    ExtATables T = build_extA_tables(G, H, L, R, S, M);

    std::cout << "OK. n=" << G.n
              << " |H|=" << hsz
              << " |L|=" << L.size()
              << " |R|=" << R.size()
              << "  (tables built)\n";

    int q; 
    if (!(std::cin >> q)) return 0;
    while (q--) {
        int a, b; std::cin >> a >> b;
        int got = multiply_extA(G, a, b, S, H, L, R, M, T);
        int ref = G.mul[a][b];
        std::cout << got;
        if (got != ref) std::cout << "  (mismatch! ref=" << ref << ")";
        std::cout << "\n";
    }
    return 0;
}
