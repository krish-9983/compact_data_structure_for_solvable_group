#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<cctype>
#include<math.h>
#include<unordered_set>
#include<climits>
using namespace std;

struct Subgroup {
    int order = 0;
    vector<vector<int>> table;
};
struct CosetReps 
{
    vector<int> left, right;
};

static string trim_local(const string &s) {
    int a = -1, b = -1;
    for (int i = 0; i < (int)s.size(); ++i) {
        if (!isspace((unsigned char)s[i])) { a = i; break; }
    }
    for (int i = (int)s.size() - 1; i >= 0; --i) {
        if (!isspace((unsigned char)s[i])) { b = i; break; }
    }
    if (a == -1 || b == -1) return string();
    return s.substr(a, b - a + 1);
}

// void printGroups(const vector<Subgroup>& groups) {
//     for (size_t k = 0; k < groups.size(); ++k) {
//         cout << "\nSubgroup " << (k+1)
//              << " (order " << groups[k].order << "):\n    ";
//         for (int j = 0; j < groups[k].order; ++j)
//             cout << (j+1) << (j+1==groups[k].order ? "\n" : "\t");
//         for (int i = 0; i < groups[k].order; ++i) {
//             cout << (i+1) << '\t';
//             for (int j = 0; j < groups[k].order; ++j)
//                 cout << groups[k].table[i][j]
//                      << (j+1==groups[k].order ? "\n" : "\t");
//         }
//     }
// }

pair<int, int> checkSub3(vector<Subgroup> & groups, int orderOfG) {
    double rootG = sqrt(orderOfG);

    int upperBound = (int)floor(rootG);
    int lowerBound = (int)ceil(rootG/2.0);
    int sml = INT_MIN;

    for(int i = 1; i < groups.size(); i++) {
        if(groups[i].order <= upperBound && groups[i].order >= lowerBound) {
            return {i,-1};
        }
        if(groups[i].order < lowerBound) {
            return {-1, i};            
        }
    }
    return {-1, -1};
}

// below functions are related to theorem 3:
//calculation of left and right represetative of subgroup in G
CosetReps coset_representative_cal(vector<Subgroup> & series, int desireSub){

const Subgroup& G = series[0];
int orderG = G.order;
const vector<vector<int>>& tableG = G.table;

const Subgroup& H = series[desireSub];
int orderH = H.order;
const vector<vector<int>>& tableH = H.table;

vector<int> elementsH = tableH[0];

vector<int> seenL(orderG,0), seenR(orderG,0);
CosetReps reps;

for(int g=0; g< orderG; g++) {
    if(seenL[g]) continue;
    reps.left.push_back(g);
    for(int h: elementsH) seenL[tableG[g][h]] = 1;
}

for(int g=0; g< orderG; g++) {
    if(seenR[g]) continue;
    reps.right.push_back(g);
    for(int h: elementsH) seenR[tableG[h][g]] = 1;
}

return reps;
}

//building right decomposition(cR and sR lists)
bool build_right_decomposition(
    int n,
    const vector<int>& H_elems,
    const vector<int>& repsR,
    const vector<vector<int>>& G_table,
    vector<int>& cR,
    vector<int>& sR
) {

int m = H_elems.size();
int k = repsR.size();

for (int repIdx = 0; repIdx < k; repIdx++) {
        int r = repsR[repIdx];
        for (int hIdx = 0; hIdx < m; hIdx++) {
            int hglob = H_elems[hIdx];
            int g = G_table[hglob][r];

            if (cR[g] == -1) {
                cR[g] = repIdx;
                sR[g] = hIdx;
            } 
        }
    }
return true;

}

//building left decomposition(sL and sL lists)
bool build_left_decomposition(
    int n,
    const vector<int>& H_elems,
    const vector<int>& repsL,
    const vector<vector<int>>& G_table,
    vector<int>& cL,
    vector<int>& sL
) {

int m = H_elems.size();
int k = repsL.size();

for (int repIdx = 0; repIdx < k; repIdx++) {
        int l = repsL[repIdx];
        for (int hIdx = 0; hIdx < m; hIdx++) {
            int hglob = H_elems[hIdx];
            int g = G_table[l][hglob];

            if (cL[g] == -1) {
                cL[g] = repIdx;
                sL[g] = hIdx;
            } 
        }
    }
return true;

}

//flip of l x h -> H' x r
void build_flip_RH(
    const vector<int>& repsL,
    const vector<int>& repsR,
    const vector<int>& H_elems,
    const vector<int>& cR, 
    const vector<int>& sR,
    const vector<vector<int>>& G_table,
    vector<vector<int>>& flipH,
    vector<vector<int>>& flipR 
) {
    int L = repsL.size();
    int m = H_elems.size();
    flipH.assign(L, vector<int>(m, -1));
    flipR.assign(L, vector<int>(m, -1));

    for (int lIdx = 0; lIdx < L; lIdx++) {
        int l = repsL[lIdx];
        for (int hIdx = 0; hIdx < m; hIdx++) {
            int h = H_elems[hIdx];
            int g = G_table[l][h];
            int rIdx = cR[g];
            int hPrimeIdx = sR[g];
            flipH[lIdx][hIdx] = hPrimeIdx;
            flipR[lIdx][hIdx] = rIdx;
        }
    }
}

//cross table RxR 
void build_cross_RR(
    const vector<int>& repsR,
    const vector<int>& H_elems, 
    const vector<int>& cR,          
    const vector<int>& sR,          
    const vector<vector<int>>& G_table,
    vector<vector<int>>& CrossH, 
    vector<vector<int>>& CrossR 
) {
    int R = repsR.size();
    int m = H_elems.size();
    CrossH.assign(R, vector<int>(R, -1));
    CrossR.assign(R, vector<int>(R, -1));

    for (int i = 0; i < R; ++i) {
        int r1 = repsR[i];
        for (int j = 0; j < R; ++j) {
            int r2 = repsR[j];
            int g = G_table[r1][r2];
            int hIdx = sR[g];
            int rPrimeIdx = cR[g];
            CrossH[i][j] = hIdx;
            CrossR[i][j] = rPrimeIdx;
        }
    }
}

// Build Fuse table H x R = g
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
        for (int rIdx = 0; rIdx < R; rIdx++) {
            int r = repsR[rIdx];
            fuse[hIdx][rIdx] = G_table[h][r];
        }
    }
}

//subgroup table for calculate multiply of 2 subgroup element which also give local indice of resulting element
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
            int hgj = H_elems[j];
            int prodg = G_table[hgi][hgj];
            int prod_hidx = global_to_Hidx[prodg];
            HH[i][j] = prod_hidx;
        }
    }
}


//muliply using compact data structure(theorem 3)
int multiply_using_cds_theorem3(
    int g1, int g2,
    const vector<int>& H_elems,
    const vector<int>& repsR,
    const vector<int>& cL,
    const vector<int>& sL,
    const vector<int>& cR,
    const vector<int>& sR,
    const vector<vector<int>>& flipH,
    const vector<vector<int>>& flipRH,
    const vector<vector<int>>& CrossH,
    const vector<vector<int>>& CrossR,
    const vector<vector<int>>& Fuse, 
    const vector<vector<int>>& HH 
) {
    int cl1 = cL[g1];
    int sl1 = sL[g1];  
    int sr1 = sR[g2]; 
    int r2Idx = cR[g2]; 

    int h1Idx = HH[sl1][sr1];

    int h2Idx = flipH[cl1][h1Idx];
    int r1Idx = flipRH[cl1][h1Idx];

    int h3Idx = CrossH[r1Idx][r2Idx];
    int r3Idx = CrossR[r1Idx][r2Idx];

    int h4Idx = HH[h2Idx][h3Idx];

    int final_global = Fuse[h4Idx][r3Idx];
    return final_global;
}


//below functions are regarding calculating theorem 4
//find identity;
int find_identity(const vector<vector<int>>& G_table) {
    int n = G_table.size();
    for (int i = 0; i < n; i++) {
        bool leftId = true, rightId = true;
        for (int j = 0; j < n; j++) {
            if (G_table[i][j] != j) leftId = false;
            if (G_table[j][i] != j) rightId = false;
            if (!leftId && !rightId) break;
        }
        if (leftId && rightId) return i;
    }
    return -1;
}

//building inverse table
vector<int> build_inverse_table(const vector<vector<int>>& G_table, int id) {
    int n = G_table.size();
    vector<int> inv(n, -1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (G_table[i][j] == id && G_table[j][i] == id) { inv[i] = j; break; }
        }
    }
    return inv;
}

//building g0_powers
vector<int> build_g0_powers(const vector<vector<int>>& G_table, int g0, int id, int maxExp) {
    vector<int> temp(maxExp+1, -1);
    temp[0] = id;
    for (int i = 1; i <= maxExp; ++i) temp[i] = G_table[temp[i-1]][g0];
    return temp;
}

//building decompostion of g into g0^ei x Sr and Sl x g0^ej and calculating lists: e, sR, sL
bool build_e_sR_sL(
    const vector<vector<int>>& G_table,
    const vector<int>& N_elems,
    int k,
    const vector<int>& inv,
    vector<int>& e,       
    vector<int>& sR, 
    vector<int>& sL,
    vector<int> & gid_to_nidx,
    vector<int>& g0pow
) {
    int n = G_table.size();
    e.assign(n, -1);
    sR.assign(n, -1);
    sL.assign(n, -1);

    for (int i = 0; i < k; ++i) {
        int g0i = g0pow[i];
        for (int ni = 0; ni < N_elems.size(); ni++) {
            int n_gid = N_elems[ni];
            int g = G_table[g0i][n_gid]; 
            e[g] = i;
            sR[g] = ni;
            int inv_g0i = inv[g0i];
            int nprime_gid = G_table[g][inv_g0i];
            sL[g] = gid_to_nidx[nprime_gid];
        }
    }
    return true;
}

//flip table which give n x g0^i -> g0^i x flip(n, i)
vector<vector<int>> build_Flip(
    const vector<vector<int>>& G_table,
    const vector<int>& N_elems,
    const vector<int>& gid_to_nidx,
    const vector<int>& inv, 
    const vector<int>& g0pow_k)
    {    
    int m = N_elems.size();
    int k = g0pow_k.size();
    vector<vector<int>> Flip(m, vector<int>(k, -1));

    for (int ni = 0; ni < m; ni++) {
        int n_gid = N_elems[ni];
        for (int i = 0; i < k; i++) {
            int g0i = g0pow_k[i]; 
            int lhs = G_table[n_gid][g0i];           
            int flip_gid = G_table[ inv[g0i] ][ lhs ];
            int flip_idx = gid_to_nidx[flip_gid];
            Flip[ni][i] = flip_idx;
        }
    }
    return Flip;
}


//to calculate lists rede and redN
void build_rede_redN(
    const vector<vector<int>>& G_table,
    const vector<int>& g0pow_full,
    const vector<int>& g0pow_k,
    const vector<int>& inv,
    const vector<int>& gid_to_nidx,
    vector<int>& rede,
    vector<int>& redN_idx    
) {
    int maxEll = g0pow_full.size() - 1;
    int k = g0pow_k.size();
    rede.assign(maxEll+1, -1);
    redN_idx.assign(maxEll+1, -1);

    for (int ell = 0; ell <= maxEll; ell++) {
        int rede_i = ell % k;          
        rede[ell] = rede_i;
        int g_ell_gid = g0pow_full[ell];
        int g_rede_gid = g0pow_k[rede_i];
        int redN_gid = G_table[inv[g_rede_gid]][g_ell_gid];
        int idx = gid_to_nidx[redN_gid];
        redN_idx[ell] = idx;
    }
}

//fuse table: 
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
        for (int ni = 0; ni < m; ni++) {
            Fuse[i][ni] = G_table[g0i][N_elems[ni]];
        }
    }
    return Fuse;
}

//multiply using cds (theorem 4)
int multiply_using_cds_theorem4(
    int g1, int g2,
    const vector<int>& e,
    const vector<int>& sR,
    const vector<int>& sL,
    const vector<vector<int>>& HH,
    const vector<vector<int>>& Flip, 
    const vector<int>& rede,
    const vector<int>& redN_idx,
    const vector<vector<int>>& Fuse
) {
    int alpha = e[g1];
    int beta  = e[g2];
    int n1_idx = HH[sR[g1]][sL[g2]];
    int n2_idx = Flip[n1_idx][beta];
    int ell = alpha + beta;
    int gamma = rede[ell];
    int redn_idx = redN_idx[ell]; 
    int n3_idx = HH[ redn_idx ][ n2_idx ];
    int result_gid = Fuse[gamma][n3_idx];
    return result_gid;
}



int main() {
    string path = "C:/Users/krish/Downloads/advanced_algorithm/goupCompactDS/cayley_parseable.txt";
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
    if (!(iss >> n)) continue;

    Subgroup g;
    g.order = n;
    g.table.assign(n, vector<int>(n));

    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            iss >> g.table[r][c];

    groups.push_back(move(g));

}
   fin.close();
pair<int, int>  isExist =  checkSub3(groups, groups[0].order);
int g1, g2;
cin >> g1 >> g2; // give 2 elements for more multiply of group g;
if (isExist.first != -1){
    //(theorem 3)
    CosetReps reps = coset_representative_cal(groups, isExist.first); 
    int orderG = groups[0].order;
    vector<int> cR(orderG, -1); // repIndex
    vector<int> sR(orderG, -1);
    vector<int> subgroupH = groups[isExist.first].table[0];
    vector<int> cL(orderG, -1); // repIndex
    vector<int> sL(orderG, -1);
    vector<vector<int>> flipH;
    vector<vector<int>> flipR;  
    vector<vector<int>> crossH, crossR;
    vector<vector<int>> fuse;
    vector<int> global_to_Hidx(orderG, -1);
    for (int i = 0; i < subgroupH.size(); ++i)
        global_to_Hidx[subgroupH[i]] = i;
    vector<vector<int>> HH;
    build_HH(subgroupH, groups[0].table, global_to_Hidx, HH);
    build_right_decomposition(orderG, subgroupH, reps.right, groups[0].table, cR, sR);
    build_left_decomposition(orderG, subgroupH, reps.left, groups[0].table, cL, sL);
    build_flip_RH(reps.left, reps.right, subgroupH, cR, sR, groups[0].table, flipH, flipR);
    build_cross_RR(reps.right, subgroupH, cR, sR, groups[0].table, crossH, crossR);
    build_fuse(subgroupH, reps.right, groups[0].table, fuse);
    int prod = multiply_using_cds_theorem3(g1, g2, subgroupH, reps.right,cL, sL, cR, sR, flipH, flipR, crossH, crossR, fuse, HH);
    cout << groups[0].table[g1][g2] << " " << prod << "     ";

    }
else {
    //thorem 4;
    int sgroupInx = isExist.second;
    const Subgroup& G = groups[0]; // group G
    const Subgroup& Nsub = groups[sgroupInx]; //normal subgroup N
    int orderG = G.order;
    int orderN = Nsub.order;
    const vector<int>& N_elems = Nsub.table[0];
    int k = orderG/orderN;  // quetient group size

    vector<int> invrse(orderG,-1);
     int id = find_identity(G.table);
    invrse = build_inverse_table(G.table, id);

    vector<int> isInN(orderG, 0); //for membership check in subgroup
    vector<int> gid_to_nid(orderG, -1);
        for (int i = 0; i < N_elems.size(); i++) {
            isInN[N_elems[i]] = 1;
            gid_to_nid[N_elems[i]] = i;
        }

    int g0 = -1; // g0(generator) for G/N calculation
    for (int i = 0; i < orderG; i++) {
            int cur = i;
            int ord = 1;
            while (ord <= k) {
                if (isInN[cur]) break;
                cur = G.table[cur][i];
                ord++;
            }
            if (ord == k && isInN[cur]) { g0 = i; break; }
        }
        int maxExp = 2*k - 2;
        vector<int> g0pow;
        vector<int> e, sR, sL;
        g0pow = build_g0_powers(G.table, g0, G.table[0][0], maxExp);
        vector<int> g0pow_k(g0pow.begin(), g0pow.begin() + k);
   
    bool ok = build_e_sR_sL(G.table, N_elems, k, invrse, e, sR, sL, gid_to_nid, g0pow_k);
vector<vector<int>> NN;
vector<vector<int>> flip;
vector<int> rede;
vector<int> redN_idx;
vector<vector<int>> fuse = build_Fuse(G.table, g0pow_k, N_elems);
build_HH(N_elems, G.table, gid_to_nid, NN);
flip = build_Flip(G.table, N_elems, gid_to_nid, invrse, g0pow);
build_rede_redN(G.table, g0pow, g0pow_k, invrse, gid_to_nid, rede, redN_idx);
int prod = multiply_using_cds_theorem4(g1,g2,e, sR, sL, NN, flip,rede, redN_idx, fuse);

cout << endl;

cout << G.table[g1][g2] << " " << prod;

    }
    return 0;

    
}
//give input 2 element of group G
//in output first first is output by actual caley table and second one is output of 0(n) dsa.
//if want to see output of theorem4. (first set if(theorem 3) to 0 so else work. then put desired index for choose subgroup from compostion series except group itself).
