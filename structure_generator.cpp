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
    int a;
    int b;
};

 int main(int argc, char* argv[]) {

    int L = 1000;      // default number of nodes
    int order = 1151;  // default (group_order_minus_1)
    int seed  = 2;     // default seed — pass different seeds for different SLPs

    if (argc == 4) {
        L     = atoi(argv[1]);
        order = atoi(argv[2]);
        seed  = atoi(argv[3]);
    }
    else if (argc == 3) {
        L     = atoi(argv[1]);
        order = atoi(argv[2]);
    }
    else if (argc != 1) {
        cerr << "Usage: ./structure_generator <L> <group_order_minus_1> [seed]\n";
        cerr << "  seed controls which random SLP is generated.\n";
        cerr << "  Same seed + increasing L = prefix-locked growing SLPs.\n";
        return 1;
    }

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

    vector<Instruction> program;
    program.reserve(L);

    vector<int> available;

    int num_const = 0;
    int num_inv = 0;
    int num_mul = 0;

    // Step 1: Create two initial constants
    for (int i = 0; i < 2; i++) {
        int value = const_dist(rng);
        program.push_back({CONST_OP, value, -1});
        available.push_back(i);
        num_const++;
    }

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
        if (chosen_op == CONST_OP) {
            int value = const_dist(rng);
            program.push_back({CONST_OP, value, -1});
            available.push_back(i);
            num_const++;
        }
        else if (chosen_op == INV_OP) {
            uniform_int_distribution<> pick(0, available.size() - 1);
            int idx = pick(rng);
            int x = available[idx];
            
            // VALIDATE: x must be < i (must exist before this instruction)
            if (x >= i) {
                cerr << "ERROR at instruction " << i << ": INV references future node " << x << "\n";
                return 1;
            }
            
            available.erase(available.begin() + idx);
            program.push_back({INV_OP, x, -1});
            available.push_back(i);
            num_inv++;
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
            
            program.push_back({MUL_OP, x, y});
            available.push_back(i);
            num_mul++;
        }
    }

    // VALIDATION: Check all references are valid
    cout << "Validating structure...\n";
    bool valid = true;
    for (size_t i = 0; i < program.size(); i++) {
        if (program[i].type == INV_OP) {
            if (program[i].a < 0 || program[i].a >= (int)i) {
                cerr << "ERROR: Instruction " << i << " (INV) references invalid node " 
                     << program[i].a << "\n";
                valid = false;
            }
        }
        else if (program[i].type == MUL_OP) {
            if (program[i].a < 0 || program[i].a >= (int)i) {
                cerr << "ERROR: Instruction " << i << " (MUL) first operand " 
                     << program[i].a << " invalid\n";
                valid = false;
            }
            if (program[i].b < 0 || program[i].b >= (int)i) {
                cerr << "ERROR: Instruction " << i << " (MUL) second operand " 
                     << program[i].b << " invalid\n";
                valid = false;
            }
        }
    }
    
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

    return 0;
}
