//aryan arora - 22MA10077 - Simple Table done on 16th August 2024
#include <bits/stdc++.h>

using namespace std; 

void print_matrix(vector<float> &A){
    int m = A.size();
    cout << "\t\t\t\t\t\t  ";
    for (int i = 0; i < m; i++){
            cout << setw(10) << A[i] << setw(10);
    }
    cout << "\n";
    return;
}

void expand_matrix(vector<vector<float>> &A, vector<float> &C) {
    int m = A.size();
    int n = A[0].size();
    // Resize each row to accommodate new variables
    for (auto &row : A) {
        row.resize(n + m, 0);
    }
    C.resize(n+m, 0);
    int slack_index = n;
    for (int i = 0; i < m; i++) 
        A[i][slack_index++] = 1;  // Add slack variable
}

class Simplex {

    private: 
        int rows, cols; 
        vector<vector<float>> A; // Coefficient matrix 
        vector<float> B; // RHS of constraints 
        vector<float> C; // Coefficient of the objective function 
        vector<int> basicVars; // To store indices of basic variables
        vector<float> D; 
        float maximum; 

    public: // Constructor 
        Simplex(vector<vector<float>> matrix, vector<float> b, vector<float> c) { 
            rows = matrix.size(); // Number of constraints 
            cols = matrix[0].size(); // Number of variables (including slack) 
            A = matrix; B = b; C = c; // Initialize basic variables (last m variables initially) 
            for (int i = 0; i < rows; i++) 
                basicVars.push_back(cols - rows + i); // Slack variables as initial basic variables 
        }
        
        // Function to print the initial Simplex Tableau 
        void printTableau() { 
            cout << "Initial Simplex Tableau:\n";
            cout << "----------------------------------------------------------------------------------------------------------------\n"; 
            print_matrix(C); 
            cout << "----------------------------------------------------------------------------------------------------------------\n"; 

            cout << setw(10) << "Basic" << setw(10) << "CB" << setw(10); 
            for (int j = 0; j < cols; j++) 
                cout << setw(10) << "x" + to_string(j + 1); 
            cout << setw(10) << "b\n"; 
            cout << "----------------------------------------------------------------------------------------------------------------\n"; 
            for (int i = 0; i < rows; i++) { 
                cout << setw(10) << "x" + to_string(basicVars[i] + 1) << setw(10) << C[basicVars[i]] << setw(10);
                for (int j = 0; j < cols; j++) { 
                    cout << setw(10) << A[i][j]; 
                } 
                cout << setw(10) << B[i] << "\n"; 
            } 
            cout << "----------------------------------------------------------------------------------------------------------------\n"; 
            cout << setw(10) << " " << setw(10) << " " << setw(10);
            calculateDelta(); 
            for (int j = 0; j < cols; j++) { 
                cout << setw(10) << D[j]; 
            }
            cout << "\n----------------------------------------------------------------------------------------------------------------\n";  
        } 
        // Function to calculate Delta (∆j) values 
        void calculateDelta() { 
            for (int j = 0; j < cols; j++) { 
                float delta = 0; 
                for (int i = 0; i < rows; i++) { 
                    delta += C[basicVars[i]] * A[i][j];
                } 
                delta -= C[j];
                D.push_back(delta); 
                //cout << "∆" << j + 1 << " = " << delta << "\n";
            } 
        } // Function to set up the initial tableau and calculate initial values 
        void setupInitialTableau() { 
            printTableau(); calculateDelta(); 
        }
        void printDelta(){
            cout << "\nDelta (∆j) values:\n"; 
            cout << "----------------------------------------------------------------------------------------------------------------\n"; 
            for (int j = 0; j < cols; j++)
                cout << "∆" << j + 1 << " = " << D[j] << "\n";
        }
}; 
                
int main() { 
    int n, m;  
    cout << "Enter no of variables (n), no of equations (m):\n";
    cin >> n >> m;

    vector<vector<float>> A(m, vector<float>(n));  // coefficients
    cout << "Enter the coefficients for each constraint (row-wise):\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cin >> A[i][j];
    }

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

    // Convert minimization to maximization if needed
    if (!isMaximize) {
        for (int i = 0; i < n; i++) {
            C[i] *= -1;
        }
    }

    expand_matrix(A, C);

    // vector<vector<float>> A = { {2, 1, 1, 1, 0}, 
    //                             {1, 3, 2, 0, 1}, 
    //                             {2, 1, 2, 0, 0}}; 
    // // Example RHS of constraints 
    // vector<float> B = {180, 300, 240}; 
    // // Example Objective function coefficients (maximize) 
    // vector<float> C = {-6, -5, -4, 0, 0}; 
    // Maximization problem, so coefficients are negative 
    Simplex simplex(A, B, C); 
    simplex.setupInitialTableau(); 
    simplex.printDelta();
    return 0; 
}
