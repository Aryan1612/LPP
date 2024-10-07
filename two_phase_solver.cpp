/*
Two_Phase_Implementation_1 (without checking for the decision vars as part of basis - 14th September 2024)
Aryan Arora
22MA10077
*/
#include <bits/stdc++.h>

using namespace std;

class two_phase {
private:
    int rows, cols;
    int total_slack{0}, total_surplus{0}, total_artificial{0};
    vector<vector<float>> A; // Coefficient matrix
    vector<float> B; // RHS of constraints
    vector<float> C; // Coefficient of the objective function
    vector<float> C_phase_I;
    vector<float> C_phase_II;
    vector<int> basicVars; // Indices of basic variables
    vector<float> D; // Delta values
    float minRatio = FLT_MAX;
    float pivotElement = FLT_MAX;
    bool two_phase_possible = 1;
    bool isMaximize = 1;
    bool is_phase2 = 0;

public:
    two_phase(vector<vector<float>> matrix, vector<float> b, vector<float> c, vector<int> unrestrictedVars, vector<string> &comp) {
        A = matrix;
        B = b;
        C = c;

        rows = matrix.size();
        cols = matrix[0].size();
        make_rhs_positive(comp);
        handleUnrestrictedVariables(unrestrictedVars); // Handle unrestricted variables
        count_variables(comp);
        expandMatrix(comp); // Expand the matrix to include slack variables, artificial vrariables and surplus variables

        rows = A.size();
        cols = A[0].size();
        basicVars = identify_basis();
        

        C_phase_I.resize(cols, 0);
        for (int l = cols-total_artificial; l < cols; l++)
            C_phase_I[l] = -1;
        run_phase_I();
    }
    void make_rhs_positive(vector<string> &comp) {
        int m = A.size();
        for (int i = 0; i < m; i++) {
            if (B[i] < 0) {
                B[i] *= -1;
                for (int j = 0; j < (int)A[i].size(); j++) {
                    A[i][j] *= -1;
                }
                // Reverse the inequality sign
                if (comp[i] == "<=") {
                    comp[i] = ">=";
                } else if (comp[i] == ">=") {
                    comp[i] = "<=";
                }
            }
        }
    }
    void count_variables(const vector<string> &comp) {
        for (const auto &c : comp) {
            if (c == "<=") {
                total_slack++;
            } else if (c == ">=") {
                total_surplus++;
                total_artificial++;
            } else if (c == "=") {
                total_artificial++;
            }
        }
    }

    void handleUnrestrictedVariables(const vector<int>& unrestrictedVars) {
        for (int index : unrestrictedVars) {
            // Add yj and zj
            for (int i = 0; i < (int) A.size(); i++) {
                A[i].push_back(-A[i][index]); // Adding -zj
            }
            C.push_back(-C[index]); // Adding -czj
            cols++; // Increase the number of columns
        }
    }

    void expandMatrix(vector<string> &comp) {
        int m = A.size();
        int n = A[0].size();

        // Resize each row to accommodate new variables
        for (auto &row : A) {
            row.resize(n + total_slack + total_surplus + total_artificial, 0);
        }

        int slack_index = n;
        int surplus_index = n + total_slack;
        int artificial_index = n + total_slack + total_surplus;

        for (int i = 0; i < m; i++) {
            if (comp[i] == "<=") {
                A[i][slack_index++] = 1;  // Add slack variable
            } else if (comp[i] == ">=") {
                A[i][surplus_index++] = -1;  // Add surplus variable
                A[i][artificial_index++] = 1;  // Add artificial variable
            } else if (comp[i] == "=") {
                A[i][artificial_index++] = 1;  // Add artificial variable
            }
        }
        C.resize(n+total_slack+total_artificial+total_surplus);
    }

    vector<int> identify_basis() {  
        vector<int> basis_indices;
        int m = A.size();
        for (int i = 0; i < m; i++){
            for (int j = cols-total_slack-total_surplus-total_artificial; j < cols; j++){
                int one_count = 0;
                if (A[i][j] == 1){
                    for (int k = 0; k < m; k++){
                        if(A[k][j] == 0)
                            continue;
                        else if (A[k][j] == 1)
                            one_count++;
                        else
                            one_count--;
                    }
                    if(one_count == 1)
                        basis_indices.push_back(j);
                }
            }
        }
        return basis_indices;
    }

    void calculateDelta(vector<float> &obj) {
        D.clear();
        for (int j = 0; j < cols; j++) {
            float delta = 0;
            for (int i = 0; i < rows; i++) {
                delta += obj[basicVars[i]] * A[i][j];
            }
            delta -= obj[j];
            if (abs(delta) < 10e-6)
                delta = 0;
            D.push_back(delta);
        }
    }

    void printMatrix(const vector<float> &matrix) {
        cout << setw(12) << " " << setw(12) << "cj";
        for (int i = 0; i < (int)matrix.size(); i++) {
            if (i >= cols - total_artificial && is_phase2)
                break;
            else
                cout << setw(12) << matrix[i];
        }
        cout << "\n";
    }

    void printbasicvars(){
        cout << "Basic Variables : \t";
        for (auto var : basicVars){
            if (var+1 <= cols-total_artificial)
                cout << "x" << var+1 << " ";
            else if(!is_phase2)
                cout << "A" << var+1-(cols-total_artificial) << " ";
        }
        cout << "\n";
        return;
    }
    void printnonbasicvars(){
        cout << "Non Basic Variables : \t";
        for (int i  = 0; i < cols; i++){
            if (find(basicVars.begin(), basicVars.end(), i)==basicVars.end()){
                if (i+1 <= cols-total_artificial)
                    cout << "x" << i+1 << " ";
                else if(!is_phase2)
                    cout << "A" << i+1-(cols-total_artificial) << " ";
            }
        }
        cout << "\n";
        return;
    }
    void printDelta(){
        cout << "\nDelta (∆j) values:\n"; 
        for (int j = 0; j < cols; j++){
            if(j < cols-total_artificial || !is_phase2)
                cout << "∆" << j + 1 << " = " << D[j] << " ";
        }
        cout << "\n";
        cout << "Most negative Delta (∆j) value : " << "∆" << (distance( D.begin(), min_element(D.begin(), D.end())) + 1) << " = " << *min_element(D.begin(), D.end()) << "\n\n";
    }
    void printTableau(vector<float> obj) {
        cout << "\n\ntwo_phase Tableau:\n";
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";

        // Print the objective function coefficients
        printMatrix(obj);
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";

        cout << setw(12) << "CB" << setw(12) << "xB";
        for (int j = 0; j < cols; j++) {
            if (j + 1 <= cols - total_artificial) {
                cout << setw(12) << "x" + to_string(j + 1);
            } else if (!is_phase2) {
                cout << setw(12) << "A" + to_string(j + 1 - (cols - total_artificial));
            }
        }
        cout << setw(12) << "b\n";
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
        for (int i = 0; i < rows; i++) {
            if (basicVars[i] + 1 <= cols - total_artificial) {
                cout << setw(12) << obj[basicVars[i]] << setw(12) << "x" + to_string(basicVars[i] + 1);
            } else if (!is_phase2) {
                cout << setw(12) << obj[basicVars[i]] << setw(12) << "A" + to_string(basicVars[i] + 1 - (cols - total_artificial));
            }
            for (int j = 0; j < cols; j++) {
                if (j >= cols - total_artificial && is_phase2) break;
                cout << setw(12) << A[i][j];
            }
            cout << setw(12) << B[i] << "\n";
        }
        cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
        calculateDelta(obj);
        cout << setw(12) << " " << setw(12) << " ";
        for (int j = 0; j < cols; j++) {
            if (j >= cols - total_artificial && is_phase2) break;
            cout << setw(12) << D[j];
        }
        cout << "\n-----------------------------------------------------------------------------------------------------------------------------------------\n";
    }


    void run_phase_I() {
        int iterations = 0;
        cout << "Phase I :";

        while (two_phase_possible) {
            iterations++;
            printTableau(C_phase_I);
            cout << "\n\nIteration = " << iterations << "\n";
            printbasicvars();
            printnonbasicvars();
            printFinalSolution(0);
            calculateDelta(C_phase_I);
            printDelta();

            int enteringVar = findEnteringVariable();
            if (enteringVar == -1) {
                cout << "Total iterations = " << iterations << "\n";
                break;
            }

            int leavingVar = findLeavingVariable(enteringVar);
            if (leavingVar == -1) {
                cout << "Unbounded solution detected!\n";
                two_phase_possible = false;
                break;
            }

            if (enteringVar > cols - total_artificial - 1)
                cout << "Entering variable : A" << enteringVar + 1 - (cols - total_artificial) << "\n";
            else
                cout << "Entering variable : x" << enteringVar + 1 << "\n";
            if (basicVars[leavingVar] > cols - total_artificial - 1)
                cout << "Leaving variable : A" << basicVars[leavingVar] + 1 - (cols - total_artificial) << "\n";
            else
                cout << "Leaving variable : x" << basicVars[leavingVar] + 1 << "\n";
            pivot(enteringVar, leavingVar);

            cout << "\nPivot element is " << pivotElement << " at (" << leavingVar + 1 << ", " << enteringVar + 1 << ")\n";
            cout << "MinRatio is : " << minRatio << "\n\n";
        }

        // Check for feasibility after Phase 1
        if (two_phase_possible) {
            for (int i = 0; i < (int)basicVars.size(); i++) {
                if (basicVars[i] > cols - total_artificial - 1) {
                    if (B[i] != 0) {
                        cout << "No feasible solution exists! Artificial variable still in basis at non-zero level!\n";
                        return;
                    }
                }
            }
            // If feasible, move to Phase 2
            is_phase2 = 1;
            cout << "\nENTERING PHASE-II\n\n";
            run_phase_II();
        }
    }

    void run_phase_II(){
        int iterations = 0;
        cout << "Phase II :\n";
        while (two_phase_possible) {
            iterations++;
            printTableau(C);
            cout << "\n\nIteration = " << iterations << "\n";
            printbasicvars();
            printnonbasicvars();
            printFinalSolution(0);
            calculateDelta(C);
            printDelta();

            int enteringVar = findEnteringVariable();
            if (enteringVar == -1) {
                cout << "Total iterations = " << iterations << "\n";
                cout << "\nFinal Solution:\n";
                break;
            }

            int leavingVar = findLeavingVariable(enteringVar);
            if (leavingVar == -1) {
                cout << "Unbounded solution detected!\n";
                two_phase_possible = false; // Stop further iterations
                break;
            }
            if (enteringVar > cols-total_artificial-1)
                cout << "Entering variable : A" << enteringVar + 1 - (cols-total_artificial) << "\n";
            else
                cout << "Entering variable : x" << enteringVar + 1 << "\n";
            if (basicVars[leavingVar] > cols-total_artificial-1)
                cout << "Leaving variable : A" << basicVars[leavingVar] + 1 - (cols-total_artificial) << "\n";
            else
                cout << "Leaving variable : x" << basicVars[leavingVar] + 1 << "\n";
            pivot(enteringVar, leavingVar);

            cout << "\nPivot element is " << pivotElement << " at (" << leavingVar + 1 << ", " << enteringVar + 1 << ")\n";
            cout << "MinRatio is : " << minRatio << "\n\n"; 
        }
        if (two_phase_possible) {
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
        vector<int> tieIndices;
        
        // First pass to find the minimum ratio
        for (int i = 0; i < rows; i++) {
            if (A[i][enteringVar] > 0) {
                float ratio = B[i] / A[i][enteringVar];
                if (ratio < minRatio) {
                    minRatio = ratio;
                    leavingVar = i;
                    tieIndices.clear(); // Clear previous ties
                    tieIndices.push_back(i);
                } else if (ratio == minRatio) {
                    tieIndices.push_back(i); // Add to tie list
                }
            }
        }
        // If there's a tie, give first preference to artificial variable then proceed further
        if (tieIndices.size() > 1) {
            cout << "\ndegeneracy detected!\n";
            for (int i = 0; i < (int) tieIndices.size(); i++){
                if(basicVars[tieIndices[i]] >= rows+total_slack+total_surplus-1)
                    return tieIndices[i];
            }
            leavingVar = *max_element(tieIndices.begin(), tieIndices.end());
        }
        return leavingVar;
    }
    void pivot(int enteringVar, int leavingVar) {
        pivotElement = A[leavingVar][enteringVar];
        for (int j = 0; j < cols; j++) 
            A[leavingVar][j] /= pivotElement;
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

    void check_alternative(){
        for (int i = 0; i < (int) D.size()-total_artificial; i++){
            if (D[i] == 0){
                if (find(basicVars.begin(), basicVars.end(), i) == basicVars.end()){
                    cout << "alternative solution found!\n";
                    int enteringVar = i;
                    int leavingVar = findLeavingVariable(enteringVar);
                    if (enteringVar > cols-total_artificial-1)
                        cout << "Entering variable : A" << enteringVar + 1 - (cols-total_artificial) << "\n";
                    else
                        cout << "Entering variable : x" << enteringVar + 1 << "\n";
                    if (basicVars[leavingVar] > cols-total_artificial-1)
                        cout << "Leaving variable : A" << basicVars[leavingVar] + 1 - (cols-total_artificial) << "\n";
                    else
                        cout << "Leaving variable : x" << basicVars[leavingVar] + 1 << "\n";
                    pivot(enteringVar, leavingVar);
                    cout << "\nPivot element is " << pivotElement << " at (" << leavingVar + 1 << ", " << enteringVar + 1 << ")\n";
                    cout << "MinRatio is : " << minRatio << "\n\n"; 
                    printTableau(C);
                    cout << "optimal solution : \n";
                    printFinalSolution(0);
                    printbasicvars();
                    printnonbasicvars();
                    printDelta();
                    break;
                }
            }
        }
    }

    void printFinalSolution(bool optimality) {
        if (optimality){
            for (int i = 0; i < (int)basicVars.size(); i++){
                if (basicVars[i] > cols-total_artificial-1){
                    if (B[i] != 0){
                        cout << "not a feasible solution!\n";
                        return;
                    }
                }
            }
            cout << "\nOptimal Solution : \n";
        }
        else
            cout << "\nBasic Feasible Solution:\n";
        for (int i = 0; i < cols - total_surplus - total_slack - total_artificial; i++) { // Only iterate over the original variables
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
        if (optimality){
            cout << "\nOptimal Objective Value: " << (isMaximize ? optimalValue : -optimalValue)  << "\n\n";
            check_alternative();
        }
        else    
            cout << "\nObjective Value: " << (isMaximize ? optimalValue : -optimalValue )<< "\n";
    }
};

int main() {
    int n, m;
    cout << "Enter number of variables (n), number of equations (m):\n";
    cin >> n >> m;

    vector<vector<float>> A(m, vector<float>(n));
    cout << "Enter the coefficients for each constraint (row-wise):\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) 
            cin >> A[i][j];
    }

    cout << "enter the equality signs for the equations:\n";
    vector<string> comp(m);
    for (int i = 0; i < m; i++)
        cin >> comp[i];

    cout << "Enter the RHS values for constraints:\n";
    vector<float> B(m);
    for (int i = 0; i < m; i++) 
        cin >> B[i];
        
    vector<float> C(n);
    cout << "Enter the coefficients of the objective function:\n";
    for (int i = 0; i < n; i++) 
        cin >> C[i];

    string opt;
    cout << "Do you want to maximize or minimize the objective function? (max/min): ";
    cin >> opt;
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
        unrestrictedVars.push_back(index - 1); // Convert to 0-based index
    }
    cout << setprecision(3);
    two_phase two_phase(A, B, C, unrestrictedVars, comp);
    return 0;
}
