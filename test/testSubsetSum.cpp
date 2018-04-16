#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <mgr++/MatIO.h>

typedef Eigen::VectorXd Set;
typedef Eigen::MatrixXd SetList;
typedef Eigen::VectorXi Indexes;

using namespace std;


Set GenRandomSet(int n, double lower, double upper) {
    Set s = Set::Random(n);
    s += Set::Ones(n);
    s /= 2.0;
    s *= abs(upper - lower);
    s += lower * Set::Ones(n);
    return s;
}

SetList GenRandomSetList(int n, int m, double lower, double upper) {
    SetList S(n,m);
    for (int i = 0; i < m; i++) {
        S.col(i) = GenRandomSet(n, lower, upper);
    }
    return S;
}

Indexes GenRandomSolution(SetList& S) {
    int rows = S.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Random(rows);
    Indexes sol = (rows * r).cast<int>().cwiseAbs(); // Sample from [0, rows)
    return sol;
}

double CalcSum(SetList& S, Indexes& sel) {
    double sum = 0;
    for (int i = 0; i < sel.size(); i++) {
        sum += S(sel[i], i);
    }
    return sum;
}

Indexes DoSubsetSum(SetList& S, double sum, double c) {
    vector<double> B(1,0); // Bins
    vector<double> Q;      // Middle sums
    int n = S.rows();
    int m = S.cols();
    for (int mi = 0; mi < m; mi++) {
        // Calc middle sums
        Q.clear();
        Q.reserve(B.size() * n);
        for (int i = 0; i < B.size(); i++)
            for (int j = 0; j < n; j++)
                Q.push_back(B[i] + S(j, mi));
        sort(Q.begin(), Q.end());
        // Discretize into bins
        B.clear();
        B.push_back(Q[0]);
        double y = B.back();
        for (int i = 1; i < Q.size(); i++) {
            if (y + c <= Q[i] && Q[i] <= sum) {
                B.push_back(Q[i]);
                y = B.back();
            }
        }
    }
    for (int i = 0; i < B.size(); i++) {
        if (abs(B[i] - sum) <= c)
            cout << "valid" << endl;
    }

    return Indexes(S.cols());
}

int main(int argc, char* argv[]) {
    SetList S = GenRandomSetList(10, 10, 0, 1);
    // cout << S << endl;
    Indexes sol = GenRandomSolution(S);
    double sum_sol = CalcSum(S, sol);
    Indexes sol_alg = DoSubsetSum(S, sum_sol, 0.01);
}