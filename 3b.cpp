//aryan arora - 22MA10077 - surplus, slack and artifical variables
#include <bits/stdc++.h>

using namespace std;

// Function to make RHS positive and adjust the inequality signs
void make_rhs_positive(vector<vector<float>> &A, vector<float> &b, vector<string> &comp) {
    int m = A.size();
    for (int i = 0; i < m; i++) {
        if (b[i] < 0) {
            b[i] *= -1;
            for (int j = 0; j < A[i].size(); j++) {
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

// Function to determine the total number of slack, surplus, and artificial variables
void count_variables(const vector<string> &comp, int &total_slack, int &total_surplus, int &total_artificial) {
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

// Function to expand the coefficient matrix to include slack, surplus, and artificial variables
void expand_matrix(vector<vector<float>> &A, const vector<string> &comp, int total_slack, int total_surplus, int total_artificial) {
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
}

// Function to identify the basis by finding the identity matrix columns
vector<int> identify_basis(const vector<vector<float>> &A) {
    vector<int> basis_indices;
    int m = A.size();
    int n = A[0].size();

    for (int j = 0; j < n; j++) {
        bool is_basis = true;
        int one_count = 0;
        for (int i = 0; i < m; i++) {
            if (A[i][j] == 1) one_count++;
            else if (A[i][j] != 0) {
                is_basis = false;
                break;
            }
        }
        if (is_basis && one_count == 1) {
            basis_indices.push_back(j);
        }
    }

    return basis_indices;
}

int main() {
    int n, m;  // m constraints, n variables;
    cout << "Enter no of variables (n), no of equations (m):\n";
    cin >> n >> m;

    vector<vector<float>> A(m, vector<float>(n));  // coefficients
    cout << "Enter the coefficients for each constraint (row-wise):\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cin >> A[i][j];
    }

    cout << "Enter '<=' | '>=' | '=' for all the m constraints:\n";
    vector<string> comp(m);
    for (int i = 0; i < m; i++)
        cin >> comp[i];

    cout << "Enter the RHS values for constraints:\n";
    vector<float> b(m);  // value on RHS of constraints
    for (int i = 0; i < m; i++)
        cin >> b[i];

    cout << "Enter 1 for x_i to be positive or 0 for unrestricted (for all values of i from 1 to n):\n";
    vector<int> sign(n);  // 1: positive, else unrestricted
    for (int i = 0; i < n; i++)
        cin >> sign[i];

    vector<float> c(n);
    cout << "Enter the coefficients of the objective function:\n";
    for (int i = 0; i < n; i++)
        cin >> c[i];

    string opt;
    cout << "Do you want to maximize or minimize the objective function? (max/min): ";
    cin >> opt;
    bool isMaximize = (opt == "max");

    // Convert minimization to maximization if needed
    if (!isMaximize) {
        for (int i = 0; i < n; i++) {
            c[i] *= -1;
        }
    }

    // Make RHS positive and adjust inequality signs
    make_rhs_positive(A, b, comp);

    // Handling unrestricted variables
    for (int i = 0; i < n; i++) {
        if (sign[i] == 0) {  // If the variable is unrestricted
            for (int j = 0; j < m; j++) {
                A[j].insert(A[j].begin() + i + 1, -A[j][i]);  // Add x2_j at the same position (negative of x_j)
            }
            c.insert(c.begin() + i + 1, -c[i]);  // Add corresponding x2_j coefficient in the objective function
            n++;
            i++;  // Skip the next index since we've added a new variable
        }
    }

    int total_slack = 0, total_surplus = 0, total_artificial = 0;
    count_variables(comp, total_slack, total_surplus, total_artificial);

    // Expand the coefficient matrix to include all slack, surplus, and artificial variables
    expand_matrix(A, comp, total_slack, total_surplus, total_artificial);

    // Identify the basis
    vector<int> basis_indices = identify_basis(A);

    // Printing the standard form
    cout << "\nStandard Form:\n";
    cout << "Maximize Z = ";
    for (int i = 0; i < c.size(); i++) {
        if (i > 0) cout << " + ";
        cout << c[i] << "x" << (i + 1);
    }
    cout << "\nSubject to:\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < A[i].size(); j++) {
            if (j > 0) cout << " + ";
            cout << A[i][j] << "x" << (j + 1);
        }
        cout << " = " << b[i] << "\n";
    }

    // Print out the variable categories
    cout << "Slack variables: ";
    for (int i = n; i < n + total_slack; i++) cout << "x" << (i + 1) << " ";
    cout << "\nSurplus variables: ";
    for (int i = n + total_slack; i < n + total_slack + total_surplus; i++) cout << "x" << (i + 1) << " ";
    cout << "\nArtificial variables: ";
    for (int i = n + total_slack + total_surplus; i < n + total_slack + total_surplus + total_artificial; i++) cout << "x" << (i + 1) << " ";
    cout << endl;

    // Print the basis
    cout << "\nBasis (columns forming the identity matrix):\n";
    if (basis_indices.empty()) {
        cout << "No basis columns found.\n";
    } else {
        for (const auto &index : basis_indices) {
            cout << "x" << (index + 1) << " ";
        }
        cout << endl;
    }

    return 0;
}
