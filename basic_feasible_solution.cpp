#include <bits/stdc++.h>
using namespace std;

class BasicFeasibleSolution {
private:
    int n, m;
    vector<vector<float>> a;
    vector<float> b;

public:
    BasicFeasibleSolution(int n, int m, vector<vector<float>> &a, vector<float> &b)
        : n(n), m(m), a(a), b(b) {}

    vector<vector<int>> getCombinations(vector<int> &arr, int r) {
        vector<vector<int>> combinations;
        vector<int> current;
        combinationsUntil(arr, current, 0, r, combinations);
        return combinations;
    }

    void solve() {
        vector<int> arr(n);
        for (int i = 0; i < n; i++)
            arr[i] = i;
        vector<vector<int>> combinations = getCombinations(arr, m);

        for (auto combination : combinations) {
            vector<float> c{b};
            vector<vector<float>> subMatrix(m, vector<float>(0));
            for (auto it : combination) {
                for (int i = 0; i < m; i++) {
                    subMatrix[i].push_back(a[i][it]);
                }
            }
            cout << fixed << setprecision(5) << "for basic variables: ";
            for (auto it : combination) {
                cout << "x_" << it << " ";
            }
            cout << "\n";
            make_dominant(subMatrix, c);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    cout << subMatrix[i][j] << " ";
                }
                cout << c[i];
                cout << "\n";
            }
            cout << "\n";
            gauss_seidel(subMatrix, c);
        }
    }

private:
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
        return sqrt(err);
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
        if ((int)st.size() < (int)indices.size()) {
            cout << "Not diagonally dominant in any case, Gauss-Seidel not applicable\n";
        } else {
            vector<vector<float>> temp{a};
            vector<float> c{b};
            for (int i = 0; i < (int)indices.size(); i++) {
                temp[indices[i]] = a[i];
                c[indices[i]] = b[i];
            }
            a = temp;
            b = c;
        }
    }

    void print(const vector<float> &x) {
        for (float val : x) {
            cout << val << " ";
        }
        cout << "\n";
    }

    void gauss_seidel(vector<vector<float>> &a, vector<float> &b) {
        int n = a.size();
        vector<float> x(n, 0);
        cout << "Diagonally dominant (0 or 1): " << diagonally_dominant(a) << "\n\n";

        if (diagonally_dominant(a)) {
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

            } while (error(x, prev_x) > 0.000001);

            cout << "Final result: ";
            print(x);
            bool isBasicFeasible = true;
            bool isDegenerate = false;

            for (float xi : x) {
                if (xi < 0) {
                    cout << "Basic solution\t";
                    isBasicFeasible = false;
                    break;
                }
            }
            if (isBasicFeasible)
                cout << "Basic feasible solution\t";

            for (float xi : x) {
                if (xi < 0.00001 && xi >= 0) {
                    cout << "Degenerate solution\n\n";
                    isDegenerate = true;
                    break;
                }
            }
            if (!isDegenerate)
                cout << "Non-degenerate solution\n\n";
        }
    }
};

int main() {
    int n, m;
    cin >> n >> m;

    vector<vector<float>> a(m, vector<float>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cin >> a[i][j];
    }

    vector<float> b(m);
    for (int i = 0; i < m; i++)
        cin >> b[i];

    BasicFeasibleSolution bfs(n, m, a, b);
    bfs.solve();

    return 0;
}
