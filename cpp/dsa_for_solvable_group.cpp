#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include<algorithm>
#include <regex>  // for splitting on the marker
using namespace std;

using Mask = vector<char>;

//group structure and inverse table;
struct Group {
    int n;
    vector<vector<int>> mul; // Cayley table: mul[a][b] = ab
    int e;                   // identity
    vector<int> inv;         // inverses

    Group(int n_, const vector<vector<int>>& m) : n(n_), mul(m) {
        compute_identity_and_inverses();
    }

    void compute_identity_and_inverses() {
        e = -1;
        for (int cand = 0; cand < n && e==-1; ++cand) {
            bool ok = true;
            for (int a=0;a<n && ok;++a) {
                if (mul[cand][a]!=a || mul[a][cand]!=a) ok=false;
            }
            if (ok) e = cand;
        }
        if (e==-1) throw runtime_error("No identity found");

        inv.assign(n,-1);
        for (int a=0;a<n;++a) {
            for (int b=0;b<n;++b) {
                if (mul[a][b]==e && mul[b][a]==e) { inv[a]=b; break; }
            }
            if (inv[a]==-1) throw runtime_error("No inverse for element");
        }
    }
};


// structur of layer info from layers.csv which we took from gap.
struct LayerInfo {
    int i;
    int sizeKi;
    int sizeKip;
    int k;        // order of quotient
    int g0;       // generator element (0-based)
    Mask Ki;
    Mask Kip;
    std::vector<int> reps; // optional, may be empty
};

//read text file getting from gap; 
vector<Mask> read_chain(const std::string& filename, int n) {
    ifstream in(filename);
    vector<Mask> chain;
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0]=='#') continue;
        istringstream ls(line);
        Mask m(n, 0);
        int id1;
        while (ls >> id1) {
            int id0 = id1 - 1; // convert to 0-based
            m[id0] = 1;
        }
        chain.push_back(move(m));
    }
    return chain;
}



// read layers.csv file
std::vector<LayerInfo> read_layers(const std::string& filename, int n) {
    std::ifstream in(filename);
    std::vector<LayerInfo> layers;
    std::string line;

    auto parse_ids = [&](const std::string& s) {
        std::istringstream ss(s);
        std::vector<int> ids;
        int id;
        while (ss >> id) ids.push_back(id - 1); // 0-based
        return ids;
    };

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        // split on ,|SEP|,
        std::vector<std::string> parts;
        std::regex sep(",\\|SEP\\|,");
        std::sregex_token_iterator it(line.begin(), line.end(), sep, -1), end;
        for (; it != end; ++it) parts.push_back(it->str());

        if (parts.size() < 3) {
            std::cerr << "Bad line: " << line << "\n";
            continue;
        }

        // --- meta part (before first ,|SEP|,)
        LayerInfo L;
        {
            std::string meta = parts[0];
            std::replace(meta.begin(), meta.end(), ',', ' ');
            std::istringstream ms(meta);
            ms >> L.i >> L.sizeKi >> L.sizeKip >> L.k >> L.g0;
            L.g0 -= 1; // 0-based
        }

        // --- Ki, Kip, reps
        auto KiIds  = parse_ids(parts[1]);
        auto KipIds = parse_ids(parts[2]);
        auto reps   = (parts.size() > 3 ? parse_ids(parts[3]) : std::vector<int>{});

        L.Ki  = Mask(n, 0);
        for (int v : KiIds) L.Ki[v] = 1;
        L.Kip = Mask(n, 0);
        for (int v : KipIds) L.Kip[v] = 1;
        L.reps = reps;

        layers.push_back(std::move(L));
    }
    return layers;
}



// pretty-print a mask as a list of IDs
void printMask(const Mask& m, bool oneBased = true) {
    std::cout << "{ ";
    for (int i = 0; i < (int)m.size(); ++i) {
        if (m[i]) {
            if (oneBased) std::cout << (i+1) << " ";
            else          std::cout << i << " ";
        }
    }
    std::cout << "}";
}


//below we are implementing theorem 4:

//power function useful for g0^t;
int power(const Group& G, int base, int exp) {
    if (exp==0) return G.e;
    if (exp<0) return power(G, G.inv[base], -exp);
    int res = G.e, b = base, e = exp;
    while (e>0) {
        if (e&1) res = G.mul[res][b];
        b = G.mul[b][b];
        e >>= 1;
    }
    return res;
}

//theorem 4 layer structure;
struct Layer {
    int k;    // quotient order
    int g0;   // generator
    Mask Ki, Kip;

    vector<int> e, sR, sL;  // splits
    vector<vector<int>> Flip;
    vector<int> rede, redN;
    vector<vector<int>> Fuse;

    Layer* lower = nullptr;

    Layer(const Group& G, const LayerInfo& info, Layer* low)
        : k(info.k), g0(info.g0), Ki(info.Ki), Kip(info.Kip),
          e(G.n,-1), sR(G.n,-1), sL(G.n,-1),
          Flip(G.n, vector<int>(info.k,-1)),
          rede(2*info.k-1,-1), redN(2*info.k-1,-1),
          Fuse(info.k, vector<int>(G.n,-1)),
          lower(low)
    {
        precompute(G);
    }

    void precompute(const Group& G) {
        // compute e, sR, sL
        for (int g=0; g<G.n; ++g) if (Ki[g]) {
            for (int t=0;t<k;++t) {
                int g0ti = power(G,g0,-t);
                int testR = G.mul[g0ti][g];
                if (Kip[testR]) {
                    e[g]=t; sR[g]=testR;
                    sL[g]=G.mul[g][g0ti];
                    break;
                }
            }
        }
        // Flip
        for (int n=0;n<G.n;++n) if (Kip[n]) {
            for (int i=0;i<k;++i) {
                int g0i=power(G,g0,i);
                int g0mi=power(G,g0,-i);
                Flip[n][i]=G.mul[g0mi][G.mul[n][g0i]];
            }
        }
        // rede, redN
        for (int l=0;l<=2*k-2;++l) {
            int g0l=power(G,g0,l);
            rede[l]=e[g0l];
            redN[l]=sR[g0l];
        }
        // Fuse
        for (int i=0;i<k;++i) {
            int g0i=power(G,g0,i);
            for (int n=0;n<G.n;++n) if (Kip[n]) {
                Fuse[i][n]=G.mul[g0i][n];
            }
        }
    }

    int multiply(const Group& G, int g1, int g2) {
        int alpha=e[g1], beta=e[g2];
        int n1 = lower? lower->multiply(G, sR[g1], sL[g2])
                       : G.mul[sR[g1]][sL[g2]];
        int n2 = Flip[n1][beta];
        int gamma=rede[alpha+beta];
        int c=redN[alpha+beta];
        int n3 = lower? lower->multiply(G,c,n2) : G.mul[c][n2];
        return Fuse[gamma][n3];
    }
};

//building layer stack
vector<Layer*> buildLayers(const Group& G, const vector<LayerInfo>& infos) {
    vector<Layer*> layers;
    Layer* low=nullptr;
    for (int i=infos.size()-1;i>=0;--i) {
        Layer* L=new Layer(G,infos[i],low);
        layers.push_back(L);
        low=L;
    }
    reverse(layers.begin(),layers.end());
    return layers;
}


int main() {
// Example Cayley table for D8 (8x8 matrix)
vector<vector<int>> mt  = {
    {0,1,2,3},
    {1,2,3,0},
    {2,3,0,1},
    {3,0,1,2}
    };
Group G(4,mt);

// load chain/layers from GAP
auto chain = read_chain("C:/Users/krish/Downloads/advanced algorithm/compact_DS_solvableGroup/gap/comp_series.txt", G.n);
auto infos = read_layers("C:/Users/krish/Downloads/advanced algorithm/compact_DS_solvableGroup/gap/layers.csv", G.n);


// build DS
auto layers=buildLayers(G,infos);

// test multiplication
bool ok=true;
for (int a=0;a<G.n;++a){
  for (int b=0;b<G.n;++b){
    int ds=layers[0]->multiply(G,a,b);
    if (ds!=G.mul[a][b]) {
        cout<<"Mismatch at "<<a<<","<<b<<"\n";
        ok=false;
    }
  }
}
cout << (ok?"Our created data structure is correct\n":"Error found, pplz check the program\n");
}

