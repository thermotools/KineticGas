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
#pragma endregion

#pragma region // Fac definition
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
#pragma endregion

#pragma region // Product definition
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
#pragma endregion

#pragma region // Frac definition
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

#pragma endregion

#pragma region // Operators
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
#pragma endregion

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