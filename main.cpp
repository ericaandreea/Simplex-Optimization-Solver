#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <string>

using namespace std;

const double EPSILON = 1e-8;
const double BIG_M = 1000000.0;

class UniversalSimplex {
private:
    vector<vector<double>> tableau;
    vector<int> basic_variables;
    int num_rows, num_cols;
    int num_orig_vars;
    int first_artificial_idx;
    int num_artificial;
    bool is_max;

    void printTableau(int iteration, int pivot_row = -1, int pivot_col = -1) {
        cout << "\nTABEL SIMPLEX - Iteratia " << iteration << endl;
        cout << string(10, ' ') << "|";
        for (int j = 0; j < num_cols - 1; ++j) {
            string var_name = "x" + to_string(j + 1);
            cout << setw(10) << var_name << "|";
        }
        cout << setw(10) << "RHS" << "|" << endl;
        cout << string(11 + num_cols * 11, '-') << endl;

        for (int i = 0; i < num_rows; ++i) {
            if (i < num_rows - 1) {
                cout << "x" << setw(2) << basic_variables[i] + 1 << " (E" << i + 1 << ") |";
            } else {
                cout << setw(10) << "f(x) / Z |";
            }

            for (int j = 0; j < num_cols; ++j) {
                string prefix = (i == pivot_row && j == pivot_col) ? "->" : "  ";
                cout << prefix;
                if (abs(tableau[i][j]) < EPSILON) cout << setw(8) << "0.00";
                else cout << setw(8) << fixed << setprecision(2) << tableau[i][j];
                cout << "|";
            }
            cout << endl;
        }
        cout << string(11 + num_cols * 11, '-') << endl;
    }

public:
    UniversalSimplex(vector<vector<double>> A, vector<double> b, vector<double> c, vector<int> signs, bool is_maximize) {
        is_max = is_maximize;
        int m = A.size();
        int n = c.size();
        num_orig_vars = n;

        for (int i = 0; i < m; ++i) {
            if (b[i] < 0) {
                b[i] = -b[i];
                for (int j = 0; j < n; ++j) A[i][j] = -A[i][j];
                signs[i] = -signs[i];
            }
        }

        int num_slack_surplus = 0;
        num_artificial = 0;
        for (int s : signs) {
            if (s != 0) num_slack_surplus++;
            if (s >= 0) num_artificial++;
        }

        first_artificial_idx = n + num_slack_surplus;
        num_rows = m + 1;
        num_cols = n + num_slack_surplus + num_artificial + 1;
        tableau.assign(num_rows, vector<double>(num_cols, 0.0));
        basic_variables.resize(m);

        int current_slack = n;
        int current_artificial = first_artificial_idx;

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) tableau[i][j] = A[i][j];
            tableau[i][num_cols - 1] = b[i];

            if (signs[i] == -1) { // <=
                tableau[i][current_slack] = 1.0;
                basic_variables[i] = current_slack;
                current_slack++;
            } else if (signs[i] == 1) { // >=
                tableau[i][current_slack] = -1.0;
                current_slack++;
                tableau[i][current_artificial] = 1.0;
                basic_variables[i] = current_artificial;
                current_artificial++;
            } else { // =
                tableau[i][current_artificial] = 1.0;
                basic_variables[i] = current_artificial;
                current_artificial++;
            }
        }

        for (int j = 0; j < n; ++j) tableau[m][j] = is_max ? -c[j] : c[j];

        int temp_artif = first_artificial_idx;
        for (int i = 0; i < m; ++i) {
            if (signs[i] >= 0) {
                tableau[m][temp_artif] = BIG_M;
                for (int j = 0; j < num_cols; ++j) {
                    tableau[m][j] -= BIG_M * tableau[i][j];
                }
                temp_artif++;
            }
        }
    }

    void solve(vector<double>& solution, double& optimal_value) {
        int iteration = 0;
        while (true) {
            int pivot_col = -1;
            double min_val = -EPSILON;
            for (int j = 0; j < num_cols - 1; ++j) {
                if (tableau[num_rows - 1][j] < min_val) {
                    min_val = tableau[num_rows - 1][j];
                    pivot_col = j;
                }
            }

            printTableau(iteration, -1, pivot_col);

            if (pivot_col == -1) break;

            int pivot_row = -1;
            double min_ratio = 1e30;
            for (int i = 0; i < num_rows - 1; ++i) {
                if (tableau[i][pivot_col] > EPSILON) {
                    double ratio = tableau[i][num_cols - 1] / tableau[i][pivot_col];
                    if (ratio < min_ratio) {
                        min_ratio = ratio;
                        pivot_row = i;
                    }
                }
            }

            if (pivot_row == -1) throw runtime_error("Problema are solutie nemarginita!");

            cout << ">>> Pivot: R" << pivot_row + 1 << ", C" << pivot_col + 1 << " | Intra x" << pivot_col + 1 << ", Iese x" << basic_variables[pivot_row] + 1 << endl;

            double pivot_val = tableau[pivot_row][pivot_col];
            for (int j = 0; j < num_cols; ++j) tableau[pivot_row][j] /= pivot_val;

            for (int i = 0; i < num_rows; ++i) {
                if (i != pivot_row) {
                    double factor = tableau[i][pivot_col];
                    for (int j = 0; j < num_cols; ++j) tableau[i][j] -= factor * tableau[pivot_row][j];
                }
            }
            basic_variables[pivot_row] = pivot_col;
            iteration++;
        }

        for (int i = 0; i < num_rows - 1; ++i) {
            if (basic_variables[i] >= first_artificial_idx && tableau[i][num_cols - 1] > EPSILON)
                throw runtime_error("Problema este infezabila.");
        }

        solution.assign(num_orig_vars, 0.0);
        for (int i = 0; i < num_rows - 1; ++i) {
            if (basic_variables[i] < num_orig_vars) solution[basic_variables[i]] = tableau[i][num_cols - 1];
        }
        optimal_value = is_max ? tableau[num_rows - 1][num_cols - 1] : -tableau[num_rows - 1][num_cols - 1];
    }
};

void validateSolution(const vector<vector<double>>& A, const vector<double>& b, const vector<int>& signs, const vector<double>& c, const vector<double>& x) {
    cout << "\n=== VALIDAREA REZULTATELOR ===\n";
    double z_calc = 0;
    for (int j = 0; j < (int)c.size(); ++j) z_calc += c[j] * x[j];
    cout << "1. Functia obiectiv calculata: " << fixed << setprecision(2) << z_calc << endl;

    bool all_ok = true;
    for (int i = 0; i < (int)A.size(); ++i) {
        double lhs = 0;
        for (int j = 0; j < (int)A[i].size(); ++j) lhs += A[i][j] * x[j];
        bool ok = (signs[i] == -1) ? (lhs <= b[i] + EPSILON) : (signs[i] == 1) ? (lhs >= b[i] - EPSILON) : (abs(lhs - b[i]) < EPSILON);
        if (!ok) all_ok = false;
    }
    cout << "2. Status restrictii: " << (all_ok ? "TOATE OK" : "EROARE") << endl;
}

int main() {
    int n, m, opt;
    cout << "=== CONFIGURARE SIMPLEX ===\n";
    cout << "Numar variabile: "; if(!(cin >> n)) return 0;
    cout << "Numar restrictii: "; cin >> m;

    vector<vector<double>> A(m, vector<double>(n));
    vector<double> b(m), c(n);
    vector<int> signs(m);

    cout << "Coeficienti functie obiectiv:\n";
    for (int j = 0; j < n; ++j) { cout << "c" << j + 1 << ": "; cin >> c[j]; }

    cout << "Restrictii (Coeficienti, Semn [-1 pt <=, 0 pt =, 1 pt >=], b):\n";
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) cin >> A[i][j];
        cin >> signs[i] >> b[i];
    }

    cout << "Tip (1-MAX, 2-MIN): "; cin >> opt;

    try {
        UniversalSimplex solver(A, b, c, signs, (opt == 1));
        vector<double> solution;
        double optimal_value;
        solver.solve(solution, optimal_value);

        cout << "\n=== REZULTAT FINAL ===\n";
        cout << "Valoare optima Z = " << optimal_value << endl;
        for (int i = 0; i < n; ++i) cout << "x" << i + 1 << " = " << solution[i] << endl;
        validateSolution(A, b, signs, c, solution);
    } catch (const exception& e) {
        cerr << "\nEroare: " << e.what() << endl;
    }

    return 0;
}
