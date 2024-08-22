// 9th August 2024 - Aryan Arora_22MA10077 - Basic Feasible Solution
#include <bits/stdc++.h>
using namespace std;

// Function prototypes
vector<vector<int>> getCombinations(vector<int> &arr, int r);
float error(const vector<float> &x0, const vector<float> &x1);
void print(const vector<float> &x);
bool diagonally_dominant(const vector<vector<float>> &a);
void make_dominant(vector<vector<float>> &a, vector<float> &b);
bool gauss_seidel(vector<vector<float>> &a, vector<float> &b, vector<float> &solution);
void printMatrix(const vector<vector<float>> &matrix, vector<float> &subRHS);
void combinationsUntil(const vector<int> &arr, vector<int> &current, int start, int r, vector<vector<int>> &combinations);

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

    vector<int> arr(n);
    for (int i = 0; i < n; i++)
        arr[i] = i;
    vector<vector<int>> combinations = getCombinations(arr, m);

    vector<float> optimalSolution(n, 0);
    float optimalValue = isMaximize ? -FLT_MAX : FLT_MAX;

    cout << fixed << setprecision(5);
    cout << "List of basic feasible solutions and their objective values:\n";

    for (auto combination : combinations) {
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
        make_dominant(subMatrix, subRHS);

        // Solve using Gauss-Seidel method
        if (gauss_seidel(subMatrix, subRHS, solution)) {
            // Check if solution is basic feasible
            print(solution);
            bool isFeasible = true;
            for (float val : solution) {
                if (val < 0) {
                    isFeasible = false;
                    break;
                }
            }

            if (isFeasible) {
                // Map the solution back to the original variable indices
                vector<float> fullSolution(n, 0);
                for (int i = 0; i < m; i++) {
                    fullSolution[combination[i]] = solution[i];
                }

                float objectiveValue = 0;
                for (int i = 0; i < n; i++) {
                    objectiveValue += fullSolution[i] * c[i];
                }

                // Print the solution and objective value
                cout << "Basic feasible solution: ";
                print(fullSolution);
                cout << "Objective value: " << objectiveValue << endl;

                // Update optimal solution and value
                if ((isMaximize && objectiveValue > optimalValue) ||
                    (!isMaximize && objectiveValue < optimalValue)) {
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

    return 0;
}

vector<vector<int>> getCombinations(vector<int> &arr, int r) {
    vector<vector<int>> combinations;
    vector<int> current;
    combinationsUntil(arr, current, 0, r, combinations);
    return combinations;
}

void combinationsUntil(const vector<int> &arr, vector<int> &current, int start, int r, vector<vector<int>> &combinations) {
    int curr_size = current.size();
    if (curr_size == r) {
        combinations.push_back(current);
        return;
    }
    for (int i = start; i < (int)arr.size(); i++) {
        current.push_back(arr[i]);
        combinationsUntil(arr, current, i + 1, r, combinations);
        current.pop_back();
    }
}

bool diagonally_dominant(const vector<vector<float>> &a) {
    int cols = a[0].size();
    int rows = a.size();
    for (int i = 0; i < rows; i++) {
        float sum = 0;
        for (int j = 0; j < cols; j++) {
            sum += abs(a[i][j]);
        }
        if (sum - abs(a[i][i]) >= abs(a[i][i]))
            return false;
    }
    return true;
}

float error(const vector<float> &x0, const vector<float> &x1) {
    float err = 0;
    int n = x0.size();
    for (int i = 0; i < n; i++) {
        err += pow((x0[i] - x1[i]), 2);
    }
    err = sqrt(err);
    return err;
}

void make_dominant(vector<vector<float>> &a, vector<float> &b) {
    int rows = a.size();
    int cols = a[0].size();
    vector<int> indices(rows);
    for (int i = 0; i < rows; i++) {
        int sum = 0;
        for (int j = 0; j < cols; j++)
            sum += abs(a[i][j]);
        for (int j = 0; j < cols; j++) {
            if (abs(a[i][j]) > sum - abs(a[i][j])) {
                indices[i] = j;
            }
        }
    }
    set<int> st{indices.begin(), indices.end()};
    if ((int)st.size() < (int)indices.size())
        cout << "Not diagonally dominant in any case, Gauss-Seidel may not converge\n";
    else {
        vector<vector<float>> temp(rows, vector<float>(cols, 0));
        vector<float> c(rows, 0);
        for (int i = 0; i < rows; i++) {
            temp[i] = a[indices[i]];
            c[i] = b[indices[i]];
        }
        a = temp;
        b = c;
    }
}

void print(const vector<float> &x) {
    int n = x.size();
    for (int i = 0; i < n; i++) {
        cout << "x_" << i + 1 << " = " << x[i] << " ";
    }
    cout << "\n";
    return;
}

void printMatrix(const vector<vector<float>> &matrix, vector<float> &subRHS) {
    int i{0};
    for (const auto &row : matrix) {
        for (float val : row) {
            cout << setw(10) << val << " ";
        }
        cout << " | " << subRHS[i++];
        cout << endl;
    }
}

bool gauss_seidel(vector<vector<float>> &a, vector<float> &b, vector<float> &solution) {
    int n = a.size();
    vector<float> x(n, 0);
    cout << "Diagonally dominant (0 or 1): " << diagonally_dominant(a) << "\n";

    if (diagonally_dominant(a)) {
        vector<float> prev_x(n, 0);  // Previous values of x
        do {
            prev_x = x;  // Update prev_x to the current values of x

            for (int j = 0; j < n; j++) {
                float sum = b[j];
                for (int k = 0; k < n; k++) {
                    if (k != j) {
                        sum -= a[j][k] * x[k];
                    }
                }
                x[j] = sum / a[j][j];
            }

        } while (abs(error(x, prev_x)) > 0.000001);

        solution = x;
        return true;
    } else {
        return false;
    }
}
