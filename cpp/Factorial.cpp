/*
Author : Vegard Gjeldvik Jervell
Contains : See Factorial.h for documentation and comments.
*/

#include "Factorial.h"

// Returns the exponential of an integer base as a Product.
Product ipow(int base, int expo){
    if (expo == 0){
        return 1; // Following pow(0, 0) = 1 convention
    }
    if (expo < 0){
        double val = pow(base, expo);
        return Product{val};
    }
    if (base == -1){
        if (expo % 2 == 0){
            return Product{1};
        }
        return Product{-1};
    }
    Product r{base};
    for (int i = 1; i < expo; i++){
        r *= base;
    }
    return r;
}

Fac::Fac(int v) : val{v} {}

long long Fac::eval(){    
    if (val == 0 || val == 1){
        return 1;
    }
    else {
        long long r = 1;
        for (int x = 1; x <= val; x++){
            r *= x;
        }
        return r;
    }
}

double Fac::eval_d(){
    if (val == 0 || val == 1){
        return 1;
    }
    double r = 1;
    for (int x = 1; x <= val; x++){
        r *= x;
    }
    return r;
}

Product::Product(){
    isize = 0;
    d = 0.0;
}

Product::Product(int const& i) : isize{1}, d{1.}{
    ilist[0] = i;
}

Product::Product(double const& d) : isize{0}, d{d}{}

Product::Product(Fac const& f) : isize{f.val - 1}, d{1.}{
    if (f.val == 0 || f.val == 1){
        isize = 0;
    }
    else{
        for (int i = 2; i <= f.val; i++){
            ilist[i-2] = i; 
        }
    }
}

double Product::eval(){
    double r = 1.0;
    for (int i = 0; i < isize; i++){
        r *= ilist[i];
    }
    return r * d;
}

Product Product::operator*=(const Product& rhs){
    for (int i = 0; i < rhs.isize; i++){
        ilist[isize + i] = rhs.ilist[i];
    }
    isize += rhs.isize;
    d *= rhs.d;
    return *this;
}

Frac::Frac(Product num, Product den) : numerator{num}, denominator{den} {}
Frac::Frac(Product num) : numerator{num}, denominator{1} {}
Frac::Frac(int i) : numerator{i}, denominator{1} {}
Frac::Frac(double d) : numerator{d}, denominator{1} {}
Frac::Frac(Fac f) : numerator{f}, denominator{1} {}

// Uses a slightly funky interation method to improve performance:
// The outer loop iterates over the integers in the numerator, and the inner loop iterates over the denominator
// Because most of the numbers in both integer arrays come from factorials, it is very common that if one number
// cancels, the next will also cancel. Therefore, if a "hit" is detected (eg. a factor in the numerator and
// denominator are equal), the numerator index is incremented, and the inner loop continues. Thus, a numerator and
// denominator of [1, 2, 3, 4, 5], [1, 2, 3, 4, 5] will only require 5 equality checks, rather than 25, because once
// the first factors are detected to cancel, the inner loop continues incrementing the outer loop while cancelling
// the remaining factors. Similarly, the numerator/denominator of [1, 2, 3, 5, 5, 5], [1, 2, 3, 5, 5, 5] will require 9
// equality checks, rather than the naive 36.
double Frac::eval(){
    for (int ni = 0; ni < numerator.isize; ni++){
        for (int di = 0; di < denominator.isize; di++){
            if (numerator.ilist[ni] == denominator.ilist[di]){
                numerator.ilist[ni] = 1;
                denominator.ilist[di] = 1;
                if (++ni == numerator.isize) break;
            }
        }
    }
    return numerator.eval() / denominator.eval();
}


Frac& Frac::operator*=(const Frac& rhs){
    numerator *= rhs.numerator;
    denominator *= rhs.denominator;
    return *this;
}

Frac& Frac::operator/=(const Frac& rhs){
    numerator *= rhs.denominator;
    denominator *= rhs.numerator;
    return *this;
}

Frac operator*(const Frac& lhs, const Frac& rhs){
    Frac r{lhs};
    return r *= rhs;
}

Frac operator/(const Frac& lhs, const Frac& rhs){
    Frac r{lhs};
    return r /= rhs;
}

double operator+(Frac& lhs, Frac& rhs){
    return lhs.eval() + rhs.eval();
}

double operator+=(double& lhs, const Frac& rhs){
    Frac r{rhs};
    lhs += r.eval();
    return lhs;
}

int factorial_tests(){
    int ia{1}, ib{2}, ic{3}, id{4};
    double da{1.5}, db{2.25}, dc{3.5}, dd{4.1};
    double tmp;

    Fac f1{ia}, f2{ic};
    Product p1{ib};
    Product p2{dc};

    Frac pt = f1 * f2;

    if (fabs(pt.eval() - f1.eval() * f2.eval()) > FLTEPS ){
        return 1;
    }
    pt = f1 * ib;
    if (fabs(pt.eval() - f1.eval() * ib) > FLTEPS){
        return 2;
    }
    pt = ic * f2;
    if (fabs(pt.eval() - f2.eval() * ic) > FLTEPS){
        return 3;
    }
    pt = p1 * p2;
    if (fabs(pt.eval() - p1.eval() * p2.eval()) > FLTEPS ){
        return 4;
    }
    tmp = pt.eval();
    pt *= p1;
    if (fabs(pt.eval() - p1.eval() * tmp) > FLTEPS){
        return 5;
    }
    pt = p1 * ib;
    if (fabs(pt.eval() - p1.eval() * ib) > FLTEPS){
        return 6;
    }
    pt = ic * p2;
    if (fabs(pt.eval() - p2.eval() * ic) > FLTEPS){
        return 7;
    }
    tmp = pt.eval();
    pt *= id;
    if (fabs(pt.eval() - tmp * id) > FLTEPS){
        return 8;
    }
    pt = p1 * da;
    if (fabs(pt.eval() - p1.eval() * da) > FLTEPS){
        return 9;
    }
    pt = da * p1;
    if (fabs(pt.eval() - p1.eval() * da) > FLTEPS){
        return 10;
    }
    tmp = pt.eval();
    pt *= dd;
    if (fabs(pt.eval() - tmp * dd) > FLTEPS){
        return 11;
    }

    pt = ipow(2,3);
    if (fabs(pt.eval() - pow(2,3)) > FLTEPS){
        return 12;
    }

    pt = ipow(-1, 16);
    if (fabs(pt.eval() - 1) > FLTEPS){
        return 13;
    }
    pt = ipow(-1, 17);
    if (fabs(pt.eval() + 1) > FLTEPS){
        return 14;
    }
    pt = ipow(8, 0);
    if (fabs(pt.eval() - 1) > FLTEPS){
        return 15;
    }

    Frac frac1{p1, p2};
    if (fabs(frac1.eval() - p1.eval() / p2.eval()) > FLTEPS){
        return 16;
    }
    Fac fac1{499}, fac2{500}, fac3{28}, fac4{30};
    Frac prod1 = fac2 * fac4;
    Frac prod2 = fac1 * fac3;
    Frac frac2 = prod1 / prod2;
    double t_val = frac2.eval() - 500 * 30 * 29;
    if (fabs(t_val) > FLTEPS){
        std::printf("Frac : %E\n", fabs(t_val));
        return 17;
    }
    Frac frac3 = db / p1;
    if (fabs(frac3.eval() - db / p1.eval()) > FLTEPS){
        return 18;
    }

    Frac frac4 = p2 / dd;
    if (fabs(frac4.eval() - p2.eval() / dd) > FLTEPS){
        return 19;
    }

    prod1 = ipow(10, 100);
    prod2 = ipow(10, 98);
    frac1 = prod1 / prod2;
    if (fabs(frac1.eval() - 100) > FLTEPS){
        return 20;
    }

    frac1 = p1 / ic;
    if (fabs(frac1.eval() - p1.eval() / ic) > FLTEPS){
        return 21;
    }
    frac1 = id / p1;
    if (fabs(frac1.eval() - id/p1.eval()) > FLTEPS){
        return 22;
    }
    prod1 = Fac{3} * 3.5;
    prod2 = Fac{3} * 5;
    frac1 = prod2 / prod1;
    frac2 = frac1 * 3.5;
    if (fabs(frac2.eval() - 5) > FLTEPS){
        return 23;
    }
    frac2 = frac1 * 7;
    if (fabs(frac2.eval() - 10) > FLTEPS){
        return 24;
    }

    prod1 = ipow(2,3);
    prod2 = Fac(4);
    frac1 = prod2 / prod1;
    tmp = 4.0;
    tmp += frac1;
    if (fabs(tmp - 7.0) > FLTEPS){
        return 25;
    }
    Frac frac5{5};
    Frac frac6{3};
    frac1 = frac6 / frac5;
    if (fabs(frac1.eval() - 3.0 / 5.0) > FLTEPS){
        return 26;
    }
    return 0;
}

double partialfactorial(int start, int stop){
    // Evaluate (stop! / start!)
    double fac = 1.;
    start = (start == 0) ? 1 : start;
    for (int i = start + 1; i <= stop; i++){
        fac *= i;
    }
    return fac;
}

double binom(int n, int k) {
    if (k > n) return 0;
    if ( (k == 0) || (k == n) ) return 1;
    if ( (k == 1) || (k == n - 1) ) return n;
    return (Fac(n) / (Fac(k) * Fac(n - k))).eval();             
}

std::vector<std::vector<int>> get_partitions(int N, int m){
    // Find all unique integer partitions of the number N, with largest value smaller than or equal to m
    // See: Memo on derivatives
    if (N == 0) return std::vector<std::vector<int>>(1, std::vector<int>());
    if (N == 1) return std::vector<std::vector<int>>(1, std::vector<int>(1, 1));
    if (m < 0) m = N;

    std::vector<std::vector<int>> partitions;
    int min_k = (N - m < 0) ? 0 : N - m;
    for (int k = min_k; k < N; k++){
        int n = N - k;
        std::vector<std::vector<int>> sub_partitions = get_partitions(k, n);
        for (std::vector<int>& part : sub_partitions){
            part.insert(part.begin(), n);
            partitions.push_back(part);
        }
    }
    return partitions;
}

long long partition_multiplicity(const std::vector<int>& partition){
    // The "Multiplicity" of a partition is the number of ways a set of N elements can be subdivided into subsets of k_1, k_2, ... elements
    // For N = k_1 + k_2 + ...
    // See: Memo on derivatives
    int n = 0;
    for (const int p : partition){
        n += p;
    }

    Fac num{n};
    Product denom{1};
    for (int p : partition){
        denom *= Fac(p);
    }

    std::vector<int> counted;
    for (int p : partition){
        if (std::count(counted.begin(), counted.end(), p) > 0) continue;
        denom *= Fac(std::count(partition.begin(), partition.end(), p));
        counted.push_back(p);
    }

    return static_cast<long long>((num / denom).eval() + 0.5); // + 0.5 To ensure correct flooring in case of funky floating point error.
}


Polynomial::Polynomial(int k_min, int k_max, std::vector<double> coeff, int k_step) 
    : k_min{k_min}, k_max{k_max}, k_step{k_step},
    _is_linear{((k_min == 0) || (k_min == 1)) && (k_max == 1)}, 
    _is_constant{(k_min == 0) && (k_max == 0)}, 
    coeff(coeff)
    {
        if ((k_max - k_min) % k_step != 0) {
            throw std::out_of_range("Polynomial: k_min - k_max is not a multiple of k_step : (" + std::to_string(k_min) + ", " + std::to_string(k_max) + ", " + std::to_string(k_step) + ")");
        }
        if (((k_max - k_min) / k_step) + 1 != coeff.size()) {
            throw std::out_of_range("Polynomial: Wrong number of coefficients (Expected " + std::to_string(((k_max - k_min) / k_step) + 1) + ", got " + std::to_string(coeff.size())
                                    + ")\n(min, max, step) = (" + std::to_string(k_min) + ", " + std::to_string(k_max) + ", " + std::to_string(k_step) + ")");
        }
    }

double Polynomial::operator()(double x) const {
    if (_is_constant) return coeff[0];
    double p = 0;
    size_t C_idx = 0;
    for (int k = k_min; k <= k_max; k += k_step){
        p += coeff[C_idx++] * pow(x, k);
    }
    return p;
}

double Polynomial::derivative(double x, int n) const {
    if (n == 0) return this->operator()(x);
    if (_is_constant) return 0;
    if (_is_linear && (n > 1)) return 0;

    double p = 0;
    int prefactor = ((n % 2 == 0) ? 1 : -1);
    int max_k_negative = (k_max < 0) ? k_max : -1;
    size_t C_idx = 0;
    int k = k_min;
    for (; k <= max_k_negative; k += k_step){
        p += prefactor * partialfactorial(- k - 1, n - k - 1) * coeff[C_idx++] * pow(x, k - n);
    }
    for (; k <= k_max; k += k_step){
        if (k < n){C_idx++; continue;}
        p += partialfactorial(k - n, k) * coeff[C_idx++] * pow(x, k - n);
     }
    return p;
}

std::ostream& operator<<(std::ostream& strm, const Polynomial& p){
    size_t C_idx = 0;
    for (int k = p.k_min; k < p.k_max; k += p.k_step){
        strm << "(" << p.coeff[C_idx++] << ", " << k << ") + ";
    }
    strm << "(" << p.coeff[C_idx++] << ", " << p.k_max << ")";
    return strm;
}

PolyExp::PolyExp(Polynomial pref, Polynomial expo)
    : pref{pref}, expo{expo}
    {}

PolyExp::PolyExp(Polynomial pref, double expo)
    : pref{pref}, expo{Polynomial::constant(expo)}
    {}

PolyExp::PolyExp(double pref, Polynomial expo)
    : pref{Polynomial::constant(pref)}, expo{expo}
    {}

double PolyExp::operator()(double x) const {
    return pref(x) * exp(expo(x));
}

double PolyExp::derivative(double x, int n) const {
    if (expo._is_constant){
        if (expo.coeff[0] == 0) return pref.derivative(x, n);
        return pref.derivative(x, n) * exp(expo.coeff[0]);
    }
    if (n == 0) {
        double p = pref(x);
        double E = exp(expo(x));
        return p * E;
    }

    vector1d df(n + 1, 0);
    vector1d dg(n + 1, 0);
    for (int k = 0; k <= n; k++){
        df[k] = pref.derivative(x, k);
        dg[k] = expo.derivative(x, k);
    }
    double val = 0;
    for (int k = 0; k <= n; k++){
        double Gk = get_Gk(x, k, dg);
        val += binom(n, k) * Gk * df[n - k];
    }
    return val * exp(dg[0]);
}

double PolyExp::get_Gk(double x, int k, const vector1d& dg) const {
    if (k == 0) return 1;
    if (expo._is_linear) return pow(expo.derivative(x, 1), k);

    // If expo only contains positive powers, we do not need to count partitions containing numbers larger than the largest exponent
    int max_dg_order = ((expo.k_min >= 0) && (expo.k_max < k)) ? expo.k_max : k;

    std::vector<std::vector<int>> partitions = get_partitions(k, max_dg_order);
    double G = 0;
    for (const std::vector<int>& partition : partitions){
        double tmp = 1;
        for (int l : partition){
            tmp *= dg[l];
        }
        if (tmp != 0){
            long long Z = partition_multiplicity(partition);
            G += Z * tmp;
        }
    }
    return G;
}

std::ostream& operator<<(std::ostream& strm, const PolyExp& p){
    strm << "[ " << p.pref << " ] * exp[ " << p.expo << " ]";
    return strm;
}