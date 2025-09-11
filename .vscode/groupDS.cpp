#include <iostream>    
#include <vector>       
#include <cstdint>     
#include <random>      
#include <algorithm>
#include<string>
using namespace std;

struct Group {
    int n; 
    vector<vector<int>> mul; 
};

struct GroupChecker {
    static int find_identity(const Group &G) {
        for (int e = 0; e < G.n; ++e) {
            bool ok = true;
            for (int x = 0; x < G.n; ++x) {
                if (G.mul[e][x] != x || G.mul[x][e] != x) { ok = false; break; }
            }
            if (ok) return e;
        }
        return -1;
    }

    static bool has_inverses(const Group &G, int e) {
        for (int a = 0; a < G.n; ++a) {
            bool found = false;
            for (int b = 0; b < G.n; ++b) {
                if (G.mul[a][b] == e && G.mul[b][a] == e) {
                    found = true; break;
                }
            }
            if (!found) return false;
        }
        return true;
    }

    static bool is_group(const Group &G) {
        int e = find_identity(G);
        if (e == -1) return false;
        if (!has_inverses(G, e)) return false;
        return true; 
    }
};

struct Subgroup {                   // subgroup structure
    vector<int> elems;        // elements of H
    vector<uint8_t> inH; // membership bitset size n
};

Subgroup make_subgroup(const Group& G, const vector<int>& elems) {
    Subgroup H;
    H.elems = elems;
    H.inH.assign(G.n, 0);
    for (int x : elems) {
        if (x < 0 || x >= G.n) {
            throw std::runtime_error("Subgroup element out of range.");
        }
        H.inH[x] = 1;
    }
    return H;
}

// Build left coset reps L: one representative per coset gH
vector<int> left_reps(const Group& G, const Subgroup& H) {
    vector<int> L;
    vector<uint8_t> seen(G.n, 0);
    for (int g = 0; g < G.n; ++g) {
        if (seen[g]) continue;
        L.push_back(g);
        for (int h : H.elems) seen[G.mul[g][h]] = 1;
    }
    return L;
}

// Build right coset reps R: one representative per coset Hg
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

struct Splits {  // structure for cl,sl,cr,sr;
    std::vector<int> cL, sL; // g = cL[g] * sL[g],  cL ∈ L, sL ∈ H
    std::vector<int> cR, sR; // g = sR[g] * cR[g],  sR ∈ H, cR ∈ R
};

//below is the function to compute split arrays
Splits compute_splits(const Group& G, const Subgroup& H,
                      const vector<int>& L, const vector<int>& R) {
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
        if (S.cL[g] == -1) throw runtime_error("Left split failed for element " + to_string(g));
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
        if (S.cR[g] == -1) throw std::runtime_error("Right split failed for element " + to_string(g));
    }

    return S;
}



int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    Group G;
    G.n = n;
    G.mul.assign(n, vector<int>(n));

   
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> G.mul[i][j];
        }
    }

    if (!GroupChecker::is_group(G)) {
        cout << "Invalid group Cayley table!\n";
        return 0;
    }

    cout << "Group loaded. Order n = " << n << "\n";
    int hsz; std::cin >> hsz;
    std::vector<int> he(hsz);
    for (int i = 0; i < hsz; ++i) std::cin >> he[i];
    Subgroup H = make_subgroup(G, he);

    // Build coset reps and splits
    auto L = left_reps(G, H);
    auto R = right_reps(G, H);
    Splits S = compute_splits(G, H, L, R);

    std::cout << "OK. n=" << G.n
              << " |H|=" << hsz
              << " |L|=" << L.size()
              << " |R|=" << R.size() << "\n";

    // Dump results for verification
    std::cout << "Left reps L: ";
    for (int x : L) std::cout << x << " ";
    std::cout << "\nRight reps R: ";
    for (int x : R) std::cout << x << " ";
    std::cout << "\n";

    for (int g = 0; g < G.n; ++g) {
        std::cout << "g=" << g
                  << "   left: g = cL(" << S.cL[g] << ")*sL(" << S.sL[g] << ")"
                  << "   right: g = sR(" << S.sR[g] << ")*cR(" << S.cR[g] << ")\n";
    }
    return 0;
}
