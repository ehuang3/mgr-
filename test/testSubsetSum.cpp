#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <Eigen/Dense>
#include <complex>
#include "gtest/gtest.h"


using namespace std;

class Node;
typedef shared_ptr<Node> NodePtr;

class Complex;
typedef shared_ptr<Complex> ComplexPtr;

struct DiEdge;
typedef shared_ptr<DiEdge> DiEdgePtr;

class Node {
public:
    vector<DiEdgePtr> in;
    vector<DiEdgePtr> out;
};

struct Graph {
    vector<ComplexPtr> V;
    vector<DiEdgePtr> E;
    ComplexPtr root;
};

class Complex : public Node {
public:
    complex<double> value;
    Complex() : value(0,0) {}
    Complex(double r, double i) : value(r, i) {}
    Complex(const complex<double>& val) : value(val) {}
    double abs() { return std::abs(value); }
    double arg() { return std::arg(value); }
    double real() { return value.real(); }
    void real(double r) { value.real(r); }
    double imag() { return value.imag(); }
    void imag(double i) { value.imag(i); }
    friend ostream& operator<<(ostream& stream, const Complex& c) {
        stream << c.value;
        return stream;
    }
};

struct DiEdge : enable_shared_from_this<DiEdge> {
    ComplexPtr src;
    ComplexPtr dest;
    ComplexPtr incr;

    static DiEdgePtr Connect(ComplexPtr src, ComplexPtr dest, ComplexPtr incr) {
        DiEdgePtr e = make_shared<DiEdge>(src, dest, incr);
        src->out.push_back(e);
        dest->in.push_back(e);
        return e;
    }

    DiEdge(ComplexPtr src_, ComplexPtr dest_, ComplexPtr incr_)
    : src(src_), dest(dest_), incr(incr_) {
    }
};

class PointSet;
typedef shared_ptr<PointSet> PointSetPtr;

class Point : public Complex {
public:
    PointSetPtr set;
    Eigen::Vector2d pt;
    int index;
    bool ignore;
    Point(PointSetPtr set_, int i, double x, double y)
    : set(set_), index(i), pt(x,y), ignore(false) {
        // Do nothing.
    }
    void recompute(const Eigen::Vector2d& c, const Eigen::Vector2d& s, double f) {
        double dist = (pt-s).norm() + (c-pt).norm();
        complex<double> phase(0, 2*M_PI/f*dist);
        this->value = exp(phase);
    }
};

typedef shared_ptr<Point> PointPtr;

class PointSet {
public:
    vector<PointPtr> pts;
    void add(PointPtr pt) {
        pts.push_back(pt);
    }
    void ignore(int i, bool flag) {
        pts[i]->ignore = flag;
    }
    void recompute(const Eigen::Vector2d& c, const Eigen::Vector2d& s, double f) {
        for (int i = 0; i < pts.size(); i++) {
            pts[i]->recompute(c,s,f);
            assert(pts[i]->abs() > 0);
        }
    }
};

class HeightField {
public:
    vector<PointSetPtr> sets;

    HeightField(const Eigen::Vector2d& c, const Eigen::Vector2d& s,
              const Eigen::VectorXd& x, const Eigen::VectorXd& y, double f) {
        int nx = x.size();
        int ny = y.size();
        sets.resize(nx);
        for (int i = 0; i < nx; i++) {
            sets[i] = make_shared<PointSet>();
            for (int j = 0; j < ny; j++) {
                sets[i]->add(make_shared<Point>(sets[i], j, x[i], y[j]));
            }
        }
        recompute(c, s, f);
    }
    void recompute(const Eigen::Vector2d& c, const Eigen::Vector2d& s, double f) {
        for (int i = 0; i < sets.size(); i++) {
            sets[i]->recompute(c, s, f);
        }
    }
};

typedef shared_ptr<HeightField> HeightFieldPtr;

Graph DoSubsetSum(HeightFieldPtr S, ComplexPtr target, double eps, bool useAbs=true) {
    // Calc real/imag min and max sums.
    vector<PointSetPtr> sets = S->sets;
    double inf = numeric_limits<double>::infinity();
    ComplexPtr cmin = make_shared<Complex>();
    ComplexPtr cmax = make_shared<Complex>();
    for (int i = 0; i < sets.size(); i++) {
        double rmin, imin = inf;
        double rmax, imax = -inf;
        for (int j = 0; j < sets[i]->pts.size(); j++) {
            if (sets[i]->pts[j]->real() < rmin)
                rmin = sets[i]->pts[j]->real();
            if (rmax < sets[i]->pts[j]->real())
                rmax = sets[i]->pts[j]->real();
            if (sets[i]->pts[j]->imag() < imin)
                imin = sets[i]->pts[j]->imag();
            if (imax < sets[i]->pts[j]->imag())
                imax = sets[i]->pts[j]->imag();
        }
        assert(rmin < rmax);
        assert(imin < imax);
        cmin->real(cmin->real() + rmin); // FIXME computes total sum min, not intermediate sum min
        cmin->imag(cmin->imag() + imin);
        cmax->real(cmax->real() + rmax);
        cmax->imag(cmax->imag() + imax);
    }
    assert(cmin->real() < cmax->real());
    assert(cmin->imag() < cmax->imag());
    cout << "cmin: " << *cmin << endl;
    cout << "cmax: " << *cmax << endl;

    Graph G;
    vector<vector<ComplexPtr> > grid;
    vector<ComplexPtr> sums;
    vector<ComplexPtr> sums_prev;
    // Initialize with zero sum.
    ComplexPtr zero = make_shared<Complex>();
    G.V.push_back(zero);
    sums_prev.push_back(zero);
    for (int si = 0; si < sets.size(); si++) {
        // Calc intermediate sums.
        sums.clear();
        int nr = round((cmax->real()-cmin->real())/eps)+1;
        int nc = round((cmax->imag()-cmin->imag())/eps)+1;
        grid.clear();
        grid.resize(nr, vector<ComplexPtr>(nc));
        std::vector<PointPtr>& pts = sets[si]->pts;
        for (int ei = 0; ei < pts.size(); ei++) {
            for (int pi = 0; pi < sums_prev.size(); pi++) {
                complex<double> sum = pts[ei]->value + sums_prev[pi]->value;
                // Calc index into grid.
                int ri = round((sum.real()-cmin->real())/eps);
                int ci = round((sum.imag()-cmin->imag())/eps);
                assert(0 <= ri && ri < nr);
                cout << "ci " << ci << endl;
                cout << "nc " << nc << endl;
                cout << "cmax " << cmax->imag() << endl;
                cout << "sum " << sum.imag() << endl;

                assert(0 <= ci && ci < nc);

                if (grid[ri][ci].get() == nullptr) {
                    grid[ri][ci] = make_shared<Complex>(cmin->real()+ri*eps,cmin->imag()+ci*eps);
                    sums.push_back(grid[ri][ci]);
                }
                // Connect.
                G.E.push_back(DiEdge::Connect(sums_prev[pi], grid[ri][ci], pts[ei]));
            }
        }
        // Add sums to graph.
        G.V.insert(G.V.end(), sums.begin(), sums.end());
        // Replace previous sums.
        sums_prev = sums;
    }

    // Check solution.
    for (int i = 0; i < sums.size(); i++) {
        if (abs(sums[i]->abs()-target->abs()) < eps/2*sets.size()) {
            cout << "Success" << endl;
            cout << *sums[i] << endl;
            cout << *target << endl;
            cout << eps/2*sets.size() << endl;
        }
    }

    return G;
}

Eigen::MatrixXd initCamera(int nc, double w, double h) {
    Eigen::MatrixXd C(2,nc);
    C.row(0).setLinSpaced(-w/2, w/2);
    C.row(1) = h * Eigen::VectorXd::Ones(nc);
    return C;
}

Eigen::Vector2d initLaser() {
    Eigen::Vector2d s;
    s << 500.0, 500.0;
    return s;
}

Eigen::MatrixXd randomHeightMap(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    Eigen::MatrixXd h(2, x.size());
    h.row(0) = x;
    for (int i = 0; i < x.size(); i++) {
        Eigen::MatrixXd r = y.size() * Eigen::MatrixXd::Random(1,1);
        int ri = r.cwiseAbs().cast<int>()(0,0);
        h(1,i) = y[ri];
    }
    return h;
}

std::vector<ComplexPtr> calcE(const Eigen::MatrixXd& C, const Eigen::Vector2d& s, 
                              const Eigen::MatrixXd& h, double f) {
    int nc = C.cols();
    int nh = h.cols();
    std::vector<ComplexPtr> E(nc);
    for (int i = 0; i < nc; i++) {
        E[i] = make_shared<Complex>();
        for (int j = 0; j < nh; j++) {
            double dist = (h.block<2,1>(0,j)-s).norm() + 
                (C.block<2,1>(0,i) - h.block<2,1>(0,j)).norm();
            complex<double> phase(0, 2*M_PI/f*dist);
            E[i]->value += exp(phase);
        }
    }
    return E;
}

TEST(CalcE, 2pts) {
    // Init camera pixels and laser.
    Eigen::MatrixXd C = initCamera(200, 400, 1000);
    Eigen::Vector2d s = initLaser();
    // Init height field.
    int nx = 5;
    int ny = 2;
    double f = 1;
    Eigen::VectorXd x(nx);
    x.setLinSpaced(-1,1);
    Eigen::VectorXd y(ny);
    y.setLinSpaced(0,f/2);
    HeightFieldPtr H = make_shared<HeightField>(C.block<2,1>(0,0), s, x, y, f);
    // Init solution.
    Eigen::MatrixXd h = randomHeightMap(x, y);
    // Calc E.
    std::vector<ComplexPtr> E = calcE(C, s, h, f);
    // Run solver.
    double eps = 0.01;
    DoSubsetSum(H, E[0], eps);
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}