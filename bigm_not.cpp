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
    vector<int> artificialVars; // Indices of artificial variables
    float minRatio = FLT_MAX;
    float pivotElement = FLT_MAX;
    bool simplex_possible = 1;
    float M = 1e6; // Large positive number for the Big M method

public:
    Simplex(vector<vector<float>> matrix, vector<float> b, vector<float> c, vector<int> unrestrictedVars, vector<string> eq) {
        rows = matrix.size();
        cols = matrix[0].size();
        A = matrix;
        B = b;
        C = c;
        make_rhs_less_than(matrix, b, eq);
        check_for_simplex(b);
        handleUnrestrictedVariables(unrestrictedVars); // Handle unrestricted variables
        addArtificialVariables(eq); // Add artificial variables for Big M method
        expandMatrix(); // Expand the matrix to include slack variables
        for (int i = 0; i < rows; i++) {
            basicVars.push_back(cols - rows + i); // Slack variables as initial basic variables
        }
    }

    void make_rhs_less_than(vector<vector<float>> &matrix, vector<float> &b, vector<string> &eq){
        for (int i = 0; i < (int) eq.size(); i++){
            if (eq[i] == ">="){
                for (int j = 0; j < (int) matrix[i].size(); j++)
                    matrix[i][j]*=-1;
                b[i]*=-1;
            }
        }
        return;
    }

    void check_for_simplex(const vector<float> &b) {
        for (int i = 0; i < (int)b.size(); i++) {
            if (b[i] < 0) {
                cout << "Infeasible solution: RHS should be non-negative for Simplex to work.\n";
                simplex_possible = false; // Stop the Simplex process
                return;
            }
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

    void addArtificialVariables(vector<string>& eq) {
        for (int i = 0; i < rows; i++) {
            if (eq[i] == ">=") {
                // Add artificial variable for >= constraints
                artificialVars.push_back(cols); 
                C.push_back(-M); // Subtract M * artificialVar in the objective function
                for (int j = 0; j < rows; j++) {
                    A[j].push_back((i == j) ? 1 : 0); // Add artificial variable to matrix
                }
                B[i] = B[i]; // Adjust RHS
                cols++;
            }
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

    void printbasicvars(){
        cout << "Basic Variables : \t";
        for (auto var : basicVars)
            cout << "x" << var+1 << " ";
        cout << "\n";
        return;
    }
    void printnonbasicvars(){
        cout << "Non Basic Variables : \t";
        for (int i  = 0; i < cols; i++){
            if (find(basicVars.begin(), basicVars.end(), i)==basicVars.end())
                cout << "x" << i+1 << " ";
        }
        cout << "\n";
        return;
    }
    void printDelta(){
        cout << "\nDelta (∆j) values:\n"; 
        for (int j = 0; j < cols; j++)
            cout << "∆" << j + 1 << " = " << D[j] << " ";
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
        int iterations = 0;
        while (simplex_possible) {
            iterations++;
            cout << "Iteration = " << iterations << "\n";
            calculateDelta();

            int enteringVar = findEnteringVariable();
            if (enteringVar == -1) {
                cout << "\n\nOptimal solution found!\n";
                cout << "Total iterations = " << iterations << "\n";
                cout << "\nFinal Optimal Solution:\n";
                break;
            }

            int leavingVar = findLeavingVariable(enteringVar);
            if (leavingVar == -1) {
                cout << "Unbounded solution detected!\n";
                simplex_possible = false; // Stop further iterations
                break;
            }

            cout << "Entering variable : x" << enteringVar + 1 << "\n";
            cout << "Leaving variable : x" << basicVars[leavingVar] + 1 << "\n";
            pivot(enteringVar, leavingVar);

            cout << "\nPivot element is " << pivotElement << " at (" << leavingVar + 1 << ", " << enteringVar + 1 << ")\n";
            cout << "MinRatio is : " << minRatio << "\n\n"; 

            printTableau();

            printFinalSolution(0);

            printbasicvars();

            printnonbasicvars();

            printDelta();
        }
        if (simplex_possible) {
            printFinalSolution(1); // Only print the optimal solution if it's feasible
        }
    }

    int findEnteringVariable() {
        float minDelta = *min_element(D.begin(), D.end());
        return (minDelta < 0) ? distance(D.begin(), find(D.begin(), D.end(), minDelta)) : -1;
    }

    int findLeavingVariable(int enteringVar) {
        minRatio = FLT_MAX;
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
        pivotElement = A[leavingVar][enteringVar];
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

    void printFinalSolution(int flag) {
        float finalSolution = 0;
        for (int i = 0; i < rows; i++) {
            if (basicVars[i] < cols - rows) {
                cout << "Variable x" << basicVars[i] + 1 << " = " << B[i] << "\n";
                finalSolution += C[basicVars[i]] * B[i];
            }
        }
        if (flag)
            cout << "Final Optimal Value = " << finalSolution << "\n";
    }
};

int main() {
    // Example Problem: Maximize Z = 3x1 + 2x2 
    vector<vector<float>> A = {{2, 1}, {1, 1}, {1, 3}};
    vector<float> B = {8, 6, 9};
    vector<float> C = {3, 2};
    vector<int> unrestrictedVars = {};
    vector<string> eq = {"<=", "<=", "<="};

    Simplex simplex(A, B, C, unrestrictedVars, eq);
    simplex.printTableau();
    simplex.iterate();

    return 0;
}
