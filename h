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
    float minRatio = FLT_MAX;
    float pivotElement = FLT_MAX;
    bool simplex_possible = true;

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
        expandMatrix(); // Expand the matrix to include slack variables
        initializeBasicVariables();
    }

    void make_rhs_less_than(vector<vector<float>> &matrix, vector<float> &b, vector<string> &eq) {
        for (int i = 0; i < eq.size(); i++) {
            if (eq[i] == ">=") {
                for (int j = 0; j < matrix[i].size(); j++)
                    matrix[i][j] *= -1;
                b[i] *= -1;
            }
        }
    }

    void check_for_simplex(const vector<float> &b) {
        for (int i = 0; i < b.size(); i++) {
            if (b[i] < 0) {
                cout << "Infeasible solution: RHS should be non-negative for Simplex to work.\n";
                simplex_possible = false;
                return;
            }
        }
    }

    void handleUnrestrictedVariables(const vector<int>& unrestrictedVars) {
        for (int index : unrestrictedVars) {
            for (int i = 0; i < rows; i++) {
                A[i].push_back(-A[i][index]);
            }
            C.push_back(-C[index]);
            cols++;
        }
    }

    void expandMatrix() {
        for (auto &row : A) {
            row.resize(cols + rows, 0);
        }
        C.resize(cols + rows, 0);
        int slack_index = cols;
        for (int i = 0; i < rows; i++) 
            A[i][slack_index++] = 1;
        cols += rows;
    }

    void initializeBasicVariables() {
        for (int i = 0; i < rows; i++) {
            basicVars.push_back(cols - rows + i);
        }
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
        cout << setw(10) << " ";
        for (float val : matrix) {
            cout << setw(10) << fixed << setprecision(2) << val;
        }
        cout << "\n";
    }

    void printbasicvars() {
        cout << "Basic Variables : \t";
        for (auto var : basicVars)
            cout << "x" << var + 1 << " ";
        cout << "\n";
    }

    void printnonbasicvars() {
        cout << "Non Basic Variables : \t";
        for (int i = 0; i < cols; i++) {
            if (find(basicVars.begin(), basicVars.end(), i) == basicVars.end())
                cout << "x" << i + 1 << " ";
        }
        cout << "\n";
    }

    void printDelta() {
        cout << "\nDelta (∆j) values:\n"; 
        for (int j = 0; j < cols; j++)
            cout << "∆" << j + 1 << " = " << fixed << setprecision(2) << D[j] << " ";
        cout << "\