// 9th August 2024 - Aryan Arora_22MA10077 - Basic Feasible Solution
#include <bits/stdc++.h>
using namespace std;

class LinearProgram {
private:
    int n, m; // Number of variables and constraints
    vector<vector<float>> a; // Coefficients matrix
    vector<float> b; // Right-hand side constants
    vector<float> c; // Coefficients of the objective function
    bool isMaximize; // Objective function type (maximize or minimize)

public:
    // Constructor to initialize the problem
    LinearProgram(int n, int m, const vector<vector<float>>& a, const vector<float>& b, const vector<float>& c, bool isMaximize)
        : n(n), m(m), a(a), b(b), c(c), isMaximize(isMaximize) {}

    // Function to find and print the optimal solution
    void solve() {
        vector<int> arr(n);
        iota(arr.begin(), arr.end(), 0);
        vector<vector<int>> combinations = getCombinations(arr, m);

        vector<float> optimalSolution(n, 0);
        float optimalValue = isMaximize ? -FLT_MAX : FLT_MAX;

        cout << fixed << setprecision(5);
        cout << "List of basic feasible solutions and their objective values:\n";

        for (auto& combination : combinations) {
            vector<vector<float>> subMatrix(m, vector<float>(m));
            vector<float> subRHS(m, 0); // Right-hand side for the sub-system
            vector<float> solution(n, 0);

            // Extract the sub-matrix and sub-RHS for current combination
            for (int i = 0; i < m; i++) {
                subRHS[i] = b[i];
                for (int j = 0; j < m; j++) {
                    subMatrix[i][j] = a[i][combination[j]];
                }
            }

            // Print the sub-matrix
            cout << "\nSub-matrix for combination ";
            for (int idx : combination) {
                cout << idx + 1 << " ";
            }
            cout << ":\n";
            printMatrix(subMatrix, subRHS);

            // Make the matrix diagonally dominant
            makeDominant(subMatrix, subRHS);

            // Solve using Gauss-Seidel method
            if (gaussSeidel(subMatrix, subRHS, solution)) {
                // Check if solution is basic feasible
                print(solution);
                bool isFeasible = all_of(solution.begin(), solution.end(), [](float val) { return val >= 0; });

                if (isFeasible) {
                    // Map the solution back to the original variable indices
                    vector<float> fullSolution(n, 0);
                    for (int i = 0; i < m; i++) {
                        fullSolution[combination[i]] = solution[i];
                    }

                    float objectiveValue = inner_product(fullSolution.begin(), fullSolution.end(), c.begin(), 0.0f);

                    // Print the solution and objective value
                    cout << "Basic feasible solution: ";
                    print(fullSolution);
                    cout << "Objective value: " << objectiveValue << endl;

                    // Update optimal solution and value
                    if ((isMaximize && objectiveValue > optimalValue) || (!isMaximize && objectiveValue < optimalValue)) {
                        optimalValue = objectiveValue;
                        optimalSolution = fullSolution;
                    }
                }
            }
        }

        // Print the optimal solution and value
        cout << "\nOptimal solution: ";
        print(optimalSolution);
        cout << "Optimal objective value: " << optimalValue << endl;
    }

private:
    vector<vector<int>> getCombinations(vector<int>& arr, int r) {
        vector<vector<int>> combinations;
        vector<int> current;
        combinationsUntil(arr, current, 0, r, combinations);
        return combinations;
    }

    void combinationsUntil(const vector<int>& arr, vector<int>& current, int start, int r, vector<vector<int>>& combinations) {
        if (current.size() == r) {
            combinations.push_back(current);
            return;
        }
        for (int i = start; i < (int)arr.size(); i++) {
            current.push_back(arr[i]);
            combinationsUntil(arr, current, i + 1, r, combinations);
            current.pop_back();
        }
    }

    bool diagonallyDominant(const vector<vector<float>>& a) {
        for (int i = 0; i < a.size(); i++) {
            float sum = accumulate(a[i].begin(), a[i].end(), 0.0f, [&](float acc, float val) { return acc + abs(val); });
            if (sum - abs(a[i][i]) >= abs(a[i][i])) {
                return false;
            }
        }
        return true;
    }

    float calculateError(const vector<float>& x0, const vector<float>& x1) {
        float err = 0;
        for (int i = 0; i < x0.size(); i++) {
            err += pow((x0[i] - x1[i]), 2);
        }
        return sqrt(err);
    }

    void makeDominant(vector<vector<float>>& a, vector<float>& b) {
        vector<int> indices(a.size());
        for (int i = 0; i < a.size(); i++) {
            int sum = accumulate(a[i].begin(), a[i].end(), 0.0f, [](float acc, float val) { return acc + abs(val); });
            for (int j = 0; j < a[i].size(); j++) {
                if (abs(a[i][j]) > sum - abs(a[i][j])) {
                    indices[i] = j;
                }
            }
        }
        set<int> st{indices.begin(), indices.end()};
        if (st.size() < indices.size()) {
            cout << "Not diagonally dominant in any case, Gauss-Seidel may not converge\n";
        } else {
            vector<vector<float>> temp(a.size(), vector<float>(a[0].size(), 0));
            vector<float> c(a.size(), 0);
            for (int i = 0; i < a.size(); i++) {
                temp[i] = a[indices[i]];
                c[i] = b[indices[i]];
            }
            a = temp;
            b = c;
        }
    }

    void print(const vector<float>& x) {
        for (int i = 0; i < x.size(); i++) {
            cout << "x_" << i + 1 << " = " << x[i] << " ";
        }
        cout << "\n";
    }

    void printMatrix(const vector<vector<float>>& matrix, vector<float>& subRHS) {
        for (int i = 0; i < matrix.size(); i++) {
            for (float val : matrix[i]) {
                cout << setw(10) << val << " ";
            }
            cout << " | " << subRHS[i] << endl;
        }
    }

    bool gaussSeidel(vector<vector<float>>& a, vector<float>& b, vector<float>& solution) {
        int n = a.size();
        vector<float> x(n, 0);
        cout << "Diagonally dominant (0 or 1): " << diagonallyDominant(a) << "\n";

        if (diagonallyDominant(a)) {
            vector<float> prev_x(n, 0);
            do {
                prev_x = x;
                for (int j = 0; j < n; j++) {
                    float sum = b[j];
                    for (int k = 0; k < n; k++) {
                        if (k != j) {
                            sum -= a[j][k] * x[k];
                        }
                    }
                    x[j] = sum / a[j][j];
                }
            } while (calculateError(x, prev_x) > 0.000001);

            solution = x;
            return true;
        }
        return false;
    }
};

int main() {
    int n, m;
    cout << "Enter the number of variables (n) and constraints (m): ";
    cin >> n >> m;

    vector<vector<float>> a(m, vector<float>(n));
    cout << "Enter the coefficients of the constraints matrix:\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cin >> a[i][j];
    }

    vector<float> b(m);
    cout << "Enter the constants on the right-hand side of the equations:\n";
    for (int i = 0; i < m; i++)
        cin >> b[i];

    vector<float> c(n);
    cout << "Enter the coefficients of the objective function:\n";
    for (int i = 0; i < n; i++)
        cin >> c[i];

    string opt;
    cout << "Do you want to maximize or minimize the objective function? (max/min): ";
    cin >> opt;
    bool isMaximize = (opt == "max");

    LinearProgram lp(n, m, a, b, c, isMaximize);
    lp.solve();

    return 0;
}
