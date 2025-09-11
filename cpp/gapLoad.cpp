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

//reading text file getting from gap; 
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


int main() {
    int n = 4; // group size (same n used in GAP)
    auto chain  = read_chain("C:/Users/krish/Downloads/advanced algorithm/thesis/gap/comp_series.txt", n);
    auto layers = read_layers("C:/Users/krish/Downloads/advanced algorithm/thesis/gap/layers.csv", n);

    std::cout << "Chain length: " << chain.size() << "\n";
    for (int i=0;i<chain.size();++i) {
        std::cout << "K" << i << " has " << std::count(chain[i].begin(), chain[i].end(), 1) << " elements\n";
    }

    std::cout << "\n--- Composition series ---\n";
for (int i = 0; i < chain.size(); ++i) {
    std::cout << "K" << i << " (size "
              << std::count(chain[i].begin(), chain[i].end(), 1)
              << "): ";
    printMask(chain[i]);   // print 1-based element IDs
    std::cout << "\n";
}


std::cout << "\n--- Layers ---\n";
for (auto& L : layers) {
    std::cout << "Layer " << L.i
              << " |Ki|=" << L.sizeKi
              << " |Kip|=" << L.sizeKip
              << " quotient=" << L.k
              << " g0=" << (L.g0+1) << "\n";  // +1 for human-friendly
    std::cout << "  Ki = "; printMask(L.Ki); std::cout << "\n";
    std::cout << "  Kip= "; printMask(L.Kip); std::cout << "\n";
    if (!L.reps.empty()) {
        std::cout << "  reps= "; printMask(Mask(n,0)); // TODO if reps exist
    }
}
}





