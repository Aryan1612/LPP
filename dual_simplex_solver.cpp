/*
Aryan Arora
22MA10077
Dual Simplex
Date - 4th October 2024
*/

#include <bits/stdc++.h>

using namespace std;

class dual_simplex{
    private:
        int rows, cols;
        vector<vector<double>> A; // Coefficient matrix
        vector<double> B; // RHS of constraints
        vector<double> C; // Coefficient of the objective function
        vector<int> basicVars; // Indices of basic variables
        vector<double> D; // Delta values
        bool dual_simplex_possible = 1;
        bool isMaximize = 1;
        bool optimality = 0;

    public:
        dual_simplex(vector<vector<double>> matrix, vector<double> b, vector<double> c, vector<int> unrestrictedVars, vector<string> &comp, bool flag){
            rows = matrix.size();
            cols = matrix[0].size();
            A = matrix;
            B = b;
            C = c;
            isMaximize = flag;
            make_rhs_less_than(A, B, comp);
            handleUnrestrictedVariables(unrestrictedVars); // Handle unrestricted variables
            expandMatrix(); // Expand the matrix to include slack variables
            for (int i = 0; i < rows; i++) {
                basicVars.push_back(cols - rows + i); // Slack variables as initial basic variables
            }
            printTableau();
            iterate();
        }
        void make_rhs_less_than(vector<vector<double>> &matrix, vector<double> &b, vector<string> &eq){
            for (int i = 0; i < (int) eq.size(); i++){
                if (eq[i] == ">="){
                    for (int j = 0; j < (int) matrix[i].size(); j++)
                        matrix[i][j]*=-1;
                    b[i]*=-1;
                    eq[i] = "<=";
                }
            }
            return;
        }
    
    void handleUnrestrictedVariables(const vector<int>& unrestrictedVars) {
        for (int index : unrestrictedVars) {
            // Add yj and zj
            for (int i = 0; i < rows; i++) {
                //A[i].push_back(-A[i][index]); // Adding -zj
                A[i].insert(A[i].begin() + index + 1, -A[i][index]);
            }
            // C.push_back(-C[index]); // Adding -czj
            C.insert(C.begin() + index + 1, -C[index]);
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
            double delta = 0;
            for (int i = 0; i < rows; i++) {
                delta += C[basicVars[i]] * A[i][j];
            }
            delta -= C[j];
            if (abs(delta) < 10e-6)
                delta = 0;
            D.push_back(delta);
        }
    }
    void printMatrix(const vector<double> &matrix) {
        cout << setw(12) << " " << setw(12) << "cj";
        for (double val : matrix) {
            cout << setw(12) << val; // Adjust width to align properly
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
        cout << "Most negative Delta (∆j) value : " << "∆" << (distance(D.begin(), min_element(D.begin(), D.end())) + 1) << " = " << *min_element(D.begin(), D.end()) << "\n\n";
    }
    void printTableau() {
        if (!dual_simplex_possible)
            return;
        cout << "\nDual Tableau:\n";
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";

        // Print the objective function coefficients
        printMatrix(C);

        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";

        cout << setw(12) << "CB" << setw(12) << "xB";
        for (int j = 0; j < cols; j++) {
            cout << setw(12) << "x" + to_string(j + 1);
        }
        cout << setw(12) << "b\n";
        
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
        for (int i = 0; i < rows; i++) {
            cout << setw(12) << C[basicVars[i]]<< setw(12) << "x" + to_string(basicVars[i] + 1) ;
            for (int j = 0; j < cols; j++) {
                cout << setw(12) << A[i][j];
            }
            cout << setw(12) << B[i] << "\n";
        }
        
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
        calculateDelta();
        cout << setw(12) << " " << setw(12) << " ";
        for (int j = 0; j < cols; j++) {
            cout << setw(12) << D[j];
        }
        cout << "\n-----------------------------------------------------------------------------------------------------------------------------------------\n";
    }

    void iterate() {
        int iterations = 0;
        while (dual_simplex_possible) {
            iterations++;
            cout << "\n\nIteration = " << iterations << "\n";
            printbasicvars();
            printnonbasicvars();
            printSolution(optimality);
            calculateDelta();
            printDelta();
            
            int leavingVar = findLeavingVariable();
            if (leavingVar == -1) {
                cout << "\nOptimal solution found!\n";
                break;
            }

            int enteringVar = findEnteringVariable(leavingVar);
            if (enteringVar == -1) {
                cout << "No feasible solution.\n";
                break;
            }

            cout << "Entering variable : x" << enteringVar + 1 << "\n";
            cout << "Leaving variable : x" << basicVars[leavingVar] + 1 << "\n";
            pivot(enteringVar, leavingVar);
            printTableau();
        }
        if (dual_simplex_possible) {
            printSolution(1);
        }
    }

    int findEnteringVariable(int leavingVar) {
        double maxRatio = -FLT_MAX;
        int enteringVar = -1;

        for (int j = 0; j < cols; j++) {
            if (A[leavingVar][j] < 0) { // Only consider negative coefficients
                double ratio = D[j] / A[leavingVar][j];
                if (ratio > maxRatio) {
                    maxRatio = ratio;
                    enteringVar = j;
                }
            }
        }
        
        if (enteringVar == -1) {
            dual_simplex_possible = false;
            cout << "No feasible solution exists.\n";
        }
        return enteringVar;
    }

    int findLeavingVariable() {
        double minB = *min_element(B.begin(), B.end());
        if (minB >= 0) {
            // All variables are non-negative, optimal solution found
            dual_simplex_possible = false;
            return -1;
        }
        return distance(B.begin(), min_element(B.begin(), B.end())); // Index of most negative basic variable
    }

    void pivot(int enteringVar, int leavingVar) {
        // Get the pivot element from the matrix A
        double pivotElement = A[leavingVar][enteringVar];
        
        // Normalize the leaving row (leavingVar) by dividing by the pivot element
        for (int j = 0; j < cols; j++) {
            A[leavingVar][j] /= pivotElement;
        }
        B[leavingVar] /= pivotElement;

        // Perform row operations to zero out the entering variable in other rows
        for (int i = 0; i < rows; i++) {
            if (i != leavingVar) {
                double factor = A[i][enteringVar];
                for (int j = 0; j < cols; j++) {
                    A[i][j] -= factor * A[leavingVar][j];
                }
                B[i] -= factor * B[leavingVar];
            }
        }

        // Update the basic variables to reflect the pivot operation
        basicVars[leavingVar] = enteringVar;
    }

    void printSolution(bool optimality) {
        if (optimality)
            cout << "\nOptimal Solution : \n";
        else
            cout << "\nBasic Feasible Solution:\n";
        for (int i = 0; i < cols - rows; i++) { // Only iterate over the original variables
            double value = 0;
            for (int j = 0; j < rows; j++) {
                if (basicVars[j] == i) {
                    value = B[j];
                    break;
                }
            }
            cout << "x" << i + 1 << " = " << value << "\n"; // Print only the original variable values
        }

        double optimalValue = 0;
        for (int i = 0; i < rows; i++) {
            optimalValue += C[basicVars[i]] * B[i];
        }
        if (optimality){
            cout << "\nOptimal Objective Value: " << ((isMaximize ? optimalValue : -optimalValue)) << "\n\n";
        }   
        else    
            cout << "\nObjective Value: " << ((isMaximize ? optimalValue : -optimalValue)) << "\n";
    }
};

int main() {
    int n, m;
    cout << "Enter number of variables (n), number of equations (m):\n";
    cin >> n >> m;

    vector<vector<double>> A(m, vector<double>(n));
    vector<string> comp(m);
    vector<double> B(m);

    cout << "Enter the equations row-wise\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) 
            cin >> A[i][j];
        cin >> comp[i];
        cin >> B[i];
    }

    string opt;
    cout << "Do you want to maximize or minimize the objective function? (max/min): ";
    cin >> opt;

    vector<double> C(n);
    cout << "Enter the coefficients of the objective function:\n";
    for (int i = 0; i < n; i++) 
        cin >> C[i];
    
    bool isMaximize = (opt == "max");
    if (!isMaximize) {
        for (int i = 0; i < n; i++) 
            C[i] *= -1;
    }

    // Input unrestricted variables
    vector<int> unrestrictedVars;
    cout << "Enter the number of unrestricted variables: ";
    int numUnrestricted{0};
    cin >> numUnrestricted;
    cout << "Enter the indices of unrestricted variables (1-based): ";
    for (int i = 0; i < numUnrestricted; i++) {
        int index;
        cin >> index;
        unrestrictedVars.push_back(index-1+i); // Convert to 0-based index
    }
    cout << setprecision(3);

    dual_simplex dual_simplex(A, B, C, unrestrictedVars, comp, isMaximize);

    return 0;
}