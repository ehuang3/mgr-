#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <Eigen/Dense>
#include <mgr++/MatIO.h>
#include <complex>


using namespace std;

class Node;
typedef std::shared_ptr<Node> NodePtr;

struct DEdge {
    NodePtr src;
    NodePtr dest;
};

class Node {
public:
    std::vector<DEdge> in;
    std::vector<DEdge> out;
};

void createDEdge(NodePtr src, NodePtr dest) {
    DEdge d;
    d.src = src;
    d.dest = dest;
    src->out.push_back(d);
    dest->in.push_back(d);
}

class Complex : public Node {
public:
    std::complex<double> value;
    double abs() { return std::abs(value); }
    double arg() { return std::arg(value); }
    double real() { return value.real(); }
    double imag() { return value.imag(); }
};

typedef std::shared_ptr<Complex> ComplexPtr;

class Set;
typedef std::shared_ptr<Set> SetPtr;

class Elem : public Complex {
public:
    SetPtr set;
    double x;
    double h;
    int index;
    bool ignore;
};

typedef std::shared_ptr<Elem> ElemPtr;

class Set {
public:
    std::vector<ElemPtr> elems;
    void ignore(int i, bool flag) {
        elems[i]->ignore = flag;
    }
};

class SetList {
public:
    std::vector<SetPtr> sets;
};

typedef std::shared_ptr<SetList> SetListPtr;

class Sum : public Complex {
};

typedef std::shared_ptr<Sum> SumPtr;
typedef std::vector<SumPtr> SumList;
typedef std::vector<std::vector<SumPtr> > SumsList;

class Solution {
public:
    SetListPtr sets;
    SumsList sigma;
    ComplexPtr root;
};

typedef std::shared_ptr<Solution> SolutionPtr;

SolutionPtr DoSubsetSum(SetListPtr sets, ComplexPtr target, double c, bool useAbs=true) {
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