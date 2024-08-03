//2nd August 2024 - Aryan Arora_22MA10077 - Basic Feasible Soultion
#include <bits/stdc++.h>
using namespace std;

vector <vector <int> > getCombinations(vector<int> &arr, int r);
float error(const vector<float> &x0, const vector<float> &x1);
void print(const vector<float> &x);
bool diagonally_dominant(const vector<vector<float>> &a);
void make_dominant(vector<vector<float>> &a, vector<float> &b);
void gauss_seidel(vector<vector<float>> &a, vector<float> &b);
void combinationsUntil(const vector<int> &arr, vector<int> &current, int start, int r, vector<vector<int>> &combinations);


int main(){
    int n, m;
    cin >> n >> m;
    vector<vector<float>> a(m, vector<float> (n));
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++)
            cin >> a[i][j];
    }
    vector<float> b(m);
    for (int i = 0; i < m; i++)
        cin >> b[i];
    
    vector<int> arr(n);
    for (int i = 0; i < n; i++)
        arr[i] = i;
    vector<vector<int>> combinations = getCombinations(arr, m);

    for (auto combination : combinations){
        vector <float> c{b};
        vector<vector<float> > arr (m, vector<float> (0));
        for (auto it:combination){
            for (int i = 0; i < m; i++){
                arr[i].push_back(a[i][it]);
            }
        }
        cout << "\nfor basic variables : ";
        for (auto it : combination){
            cout << "x_"<<it << " ";
        }
        cout << "\n";
        make_dominant(arr, c);
        for (int i = 0; i < m; i++){
            for (int j = 0; j < m; j++){
                cout << arr[i][j] << " ";
            }
            cout << c[i];
            cout << "\n";
        }
        cout << "\n";
        gauss_seidel(arr, c);
    }
}


vector <vector <int> > getCombinations(vector<int> &arr, int r){
    vector<vector<int>> combinations;
    vector<int> current;
    combinationsUntil(arr, current, 0, r, combinations);
    return combinations;
}
void combinationsUntil(const vector<int> &arr, vector<int> &current, int start, int r, vector<vector<int>> &combinations){
    int curr_size = current.size();
    if (curr_size == r){
        combinations.push_back(current);
        return;
    }
    for (int i = start; i < (int) arr.size(); i++){
        current.push_back(arr[i]);
        combinationsUntil(arr, current, i+1, r, combinations);
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
void make_dominant(vector<vector<float>> &a, vector<float> &b){
    int rows = a.size();
    int cols = a[0].size();
    vector<int> indices(rows);
    for (int i = 0; i < rows; i++){
        int sum = 0;
        for (int j = 0; j < cols; j++)
            sum+= abs(a[i][j]);
        for (int j = 0; j < cols; j++){
            if (abs(a[i][j]) > sum - abs(a[i][j])){
                indices[i] = j;
            }
        }
    }
    set<int> st{indices.begin(), indices.end()};
    if ((int)st.size() < (int)indices.size())
        cout << "not diagonlly dominant in any case, hence no solution for this case\n";
    else{
        vector< vector<float> > temp{a};
        vector<float> c{b};
        for (int i = 0; i < (int)indices.size(); i++){
            temp[indices[i]] = a[i];
            c[indices[i]] = b[i];
        }
        a = temp;
        b = c;
    }
}
void print(const vector<float> &x){
    int n = x.size();
    for (int i = 0; i < n; i++){
        cout << x[i] << " ";
     }
    cout << "\n";
    return;
}

void gauss_seidel(vector<vector<float>> &a, vector<float> &b){
    int n = a.size();
    vector<float> x(n, 0);
    cout << "diagonally dominant (0 or 1) " << diagonally_dominant(a) << "\n\n";

    if (diagonally_dominant(a)){
        // Input vector b
        for (int i = 0; i < n; i++){
            cin >> b[i];
        }

        vector<float> prev_x(n, 0);  // Previous values of x
        do {
            prev_x = x;  // Update prev_x to the current values of x

            for (int j = 0; j < n; j++){
                float sum = b[j];
                for (int k = 0; k < n; k++){
                    if (k != j){
                        sum -= a[j][k] * x[k];
                    }
                }
                x[j] = sum / a[j][j];
            }

            // cout << "Iteration " << iteration++ << ": ";
            // for (int i = 0; i < n; i++){
            //     cout << x[i] << " ";
            // }
            //cout << "\n";

        } while (error(x, prev_x) > 0.000001);

        cout << "Final result: ";
        for (int i = 0; i < n; i++){
            cout << x[i] << " ";
        }
        cout << "\n";
        int flag1{1}, flag2{1};
        for (int i = 0; i < n; i++){
            if (x[i] < 0){
                cout << "basic solution\t";
                flag1 = -1;
                break;
            }
        }
        if (flag1 == 1)
            cout << "basic feasible solution\t";
        for (int i = 0; i < n; i++){
            if (x[i] < 0.00001 && x[i] >= 0 ){
                cout << "degenerate solution\n";
                flag2 = -1;
                break;
            }
        }
        if (flag2 == 1)
            cout << "non-degenerate solution\n";
    }
}
