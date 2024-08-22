#include <bits/stdc++.h>

using namespace std;

class Simplex {
private:
    int rows, cols;
    vector<vector<float>> A; // Coefficient matrix
    vector<float> B; // RHS of constraints
    vector<float> C; // Coefficient of the objective function
    vector<int> basicVars; // Indices of basic variables
    vector<float> D; // Delta values
    int iterations;

public:
    Simplex(vector<vector<float>> matrix, vector<float> b, vector<float> c, vector<int> unrestrictedVars) {
        rows = matrix.size();
        cols = matrix[0].size();
        A = matrix;
        B = b;
        C = c;
        iterations = 0;

        handleUnrestrictedVariables(unrestrictedVars); // Handle unrestricted variables
        expandMatrix(); // Expand the matrix to include slack variables
        for (int i = 0; i < rows; i++) {
            basicVars.push_back(cols - rows + i); // Slack variables as initial basic variables
        }
    }

    void handleUnrestrictedVariables(const vector<int>& unrestrictedVars) {
        for (int index : unrestrictedVars) {
            // Add yj and zj
            for (int i = 0; i < rows; i++) {
                A[i].push_back(-A[i][index]); // Adding -zj
            }
            C.push_back(-C[index]); // Adding -czj
            cols++; // Increase the number of columns
        }
    }

    void expandMatrix() {
        for (auto &row : A) {
            row.resize(cols + rows, 0);
        }
        C.resize(cols + rows, 0);
        int slack_index = cols;
        for (int i = 0; i < rows; i++) 
            A[i][slack_index++] = 1; // Add slack variable
        cols += rows;
    }

    void calculateDelta() {
        D.clear();
        for (int j = 0; j < cols; j++) {
            float delta = 0;
            for (int i = 0; i < rows; i++) {
                delta += C[basicVars[i]] * A[i][j];
            }
            delta -= C[j];
            D.push_back(delta);
        }
    }

    void printMatrix(const vector<float> &matrix) {
        cout << setw(10) << " " << setw(10) << " ";
        for (float val : matrix) {
            cout << setw(10) << val; // Adjust width to align properly
        }
        cout << "\n";
    }


    void printTableau() {
        cout << "Simplex Tableau:\n";
        cout << "----------------------------------------------------------------------------------------------------------------\n";

        // Print the objective function coefficients
        printMatrix(C);

        cout << "----------------------------------------------------------------------------------------------------------------\n";

        cout << setw(10) << "Basic" << setw(10) << "CB";
        for (int j = 0; j < cols; j++) {
            cout << setw(10) << "x" + to_string(j + 1);
        }
        cout << setw(10) << "b\n";
        
        cout << "----------------------------------------------------------------------------------------------------------------\n";
        for (int i = 0; i < rows; i++) {
            cout << setw(10) << "x" + to_string(basicVars[i] + 1) << setw(10) << C[basicVars[i]];
            for (int j = 0; j < cols; j++) {
                cout << setw(10) << A[i][j];
            }
            cout << setw(10) << B[i] << "\n";
        }
        
        cout << "----------------------------------------------------------------------------------------------------------------\n";
        calculateDelta();
        cout << setw(10) << " " << setw(10) << " ";
        for (int j = 0; j < cols; j++) {
            cout << setw(10) << D[j];
        }
        cout << "\n----------------------------------------------------------------------------------------------------------------\n";
    }


    void iterate() {
        while (true) {
            iterations++;
            int enteringVar = findEnteringVariable();
            if (enteringVar == -1) {
                cout << "Optimal solution found!\n";
                break;
            }

            int leavingVar = findLeavingVariable(enteringVar);
            if (leavingVar == -1) {
                cout << "Unbounded solution!\n";
                break;
            }

            pivot(enteringVar, leavingVar);
            printTableau();
        }
        printFinalSolution();
    }

    int findEnteringVariable() {
        for (int j = 0; j < cols; j++) {
            if (D[j] < 0) return j;
        }
        return -1; // No negative delta means optimal solution is found
    }

    int findLeavingVariable(int enteringVar) {
        float minRatio = FLT_MAX;
        int leavingVar = -1;
        for (int i = 0; i < rows; i++) {
            if (A[i][enteringVar] > 0) {
                float ratio = B[i] / A[i][enteringVar];
                if (ratio < minRatio) {
                    minRatio = ratio;
                    leavingVar = i;
                }
            }
        }
        return leavingVar;
    }

    void pivot(int enteringVar, int leavingVar) {
        float pivotElement = A[leavingVar][enteringVar];
        for (int j = 0; j < cols; j++) {
            A[leavingVar][j] /= pivotElement;
        }
        B[leavingVar] /= pivotElement;

        for (int i = 0; i < rows; i++) {
            if (i != leavingVar) {
                float factor = A[i][enteringVar];
                for (int j = 0; j < cols; j++) {
                    A[i][j] -= factor * A[leavingVar][j];
                }
                B[i] -= factor * B[leavingVar];
            }
        }

        basicVars[leavingVar] = enteringVar;
    }

    void printFinalSolution() {
        cout << "\nFinal Optimal Solution:\n";
        for (int i = 0; i < cols - rows; i++) { // Only iterate over the original variables
            float value = 0;
            for (int j = 0; j < rows; j++) {
                if (basicVars[j] == i) {
                    value = B[j];
                    break;
                }
            }
            cout << "x" << i + 1 << " = " << value << "\n"; // Print only the original variable values
        }

        float optimalValue = 0;
        for (int i = 0; i < rows; i++) {
            optimalValue += C[basicVars[i]] * B[i];
        }
        cout << "Optimal Objective Value: " << optimalValue << "\n";
        cout << "Number of Iterations: " << iterations << "\n";
    }
};

int main() {
    int n, m;
    cout << "Enter number of variables (n), number of equations (m):\n";
    cin >> n >> m;

    vector<vector<float>> A(m, vector<float>(n));
    cout << "Enter the coefficients for each constraint (row-wise):\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the RHS values for constraints:\n";
    vector<float> B(m);
    for (int i = 0; i < m; i++) {
        cin >> B[i];
    }

    vector<float> C(n);
    cout << "Enter the coefficients of the objective function:\n";
    for (int i = 0; i < n; i++) {
        cin >> C[i];
    }

    string opt;
    cout << "Do you want to maximize or minimize the objective function? (max/min): ";
    cin >> opt;
    bool isMaximize = (opt == "max");

    if (!isMaximize) {
        for (int i = 0; i < n; i++) {
            C[i] *= -1;
        }
    }

    // Input unrestricted variables
    vector<int> unrestrictedVars;
    cout << "Enter the number of unrestricted variables: ";
    int numUnrestricted;
    cin >> numUnrestricted;
    cout << "Enter the indices of unrestricted variables (1-based): ";
    for (int i = 0; i < numUnrestricted; i++) {
        int index;
        cin >> index;
        unrestrictedVars.push_back(index - 1); // Convert to 0-based index
    }

    Simplex simplex(A, B, C, unrestrictedVars);
    simplex.printTableau();
    simplex.iterate();

    return 0;
}