/*
Author: Vegard Gjeldvik Jervell
Contains: Classes to facilitate analytic simplification of fractions containing factorials to prevent overflow.
            The Fac, Product and Frac classes represent analytic factorials, products and fractions. They overload
            the relevant arithmetic operators (* and /), and are simplified analytically before numerical evaluation when
            evaluated with <expression.eval()>

          Frac can be initialized from both Fac and Product, as a fraction with denominator = 1. Thus, all arithmetic
          is only defined for the Frac class, and when eg. computing Fac * Product, the result is implicitly converted
          to a fraction.

          The function 'ipow(int, int) is also defined. This function is an exponential that returns a Product (not
          evaluated to an int), to ensure that Fractions containing powers can be simplified before evaluation.

Working principle:
        A factorial is represented as a list of integers, up until the point of evaluation. Similarly, a Product is
        represented as a list of integers, and a single double, up until evaluation. A fraction is represented as
        a numerator and a denominator that are both Products.

        Upon evaluating a fraction (using Frac::eval()) common integers in the numerator and denominator are cancelled
        (in practice: set to 1), before both the numerator and denominator are evaluated (using Product::eval()) and the
        results are divided by each other, returning a double.

        This makes it possible to eg. evaluate Fac(100) / Fac(99) without any fear of overflow, and with minimal
        performance issues.

Note: The max number of integers in a Product or a factorial (Fac) is 1000. Exceeding this number leads to undefined
        behaviour (likely a segfault).
      The 'Frac' class cannot contain a sum of products in the numerator or the denominator. Both the numerator and the
      denominator must be pure products.
*/
#pragma once
#include "KineticGas.h"
#include "global_params.h"
#include "math.h"
#include "cmath"
#include <array>
#include <cstdio>

class Fac{
    public:
    int val;

    Fac(int v);
    long long eval();
    double eval_d();
};

class Product{
    public:
    int isize;
    std::array<int, 1000> ilist{};
    double d;

    Product();
    Product(int const& i);
    Product(double const& d);
    Product(Fac const& f);

    double eval();

    Product operator*=(const Product& rhs);
};

class Frac{
    public:
    Product numerator;
    Product denominator;

    Frac(){};
    Frac(int i);
    Frac(double d);
    Frac(Fac f);
    Frac(Product num, Product den);
    Frac(Product num);

    double eval();

    Frac& operator*=(const Frac& rhs);
    Frac& operator/=(const Frac& rhs);
};

Frac operator*(const Frac& lhs, const Frac& rhs);
Frac operator/(const Frac& lhs, const Frac& rhs);
double operator+(Frac& lhs, Frac& rhs);
double operator+=(double& lhs, const Frac& rhs);

Product ipow(int base, int expo);
int factorial_tests();

double partialfactorial(int start, int stop);

std::vector<std::vector<int>> get_partitions(int N, int m);
size_t partition_multiplicity(const std::vector<int>& partition);

class Term{
public:
    virtual double operator()(double x){
        return derivative(x, 0);
    }

    virtual double derivative(double x, size_t n) = 0;
}

class Polynomial : public Term {
public:
    Polynomial(int k_min, int k_max, std::vector<double> coeff, int k_step=1) 
        : k_min{k_min}, k_max{k_max}, 
        _is_linear{(k_min == 0) && (k_max == 1)}, 
        _is_constant{(k_min == 0) && (k_max == 0)}, 
        coeff(coeff)
    {}

    double operator()(double x) override;
    double derivative(double x, size_t n) override;
    inline bool is_linear(){return _is_linear;}
    inline bool is_constant(){return _is_constant;}

private:
    int k_min, k_max;
    bool _is_linear, _is_constant;
    std::vector<double> coeff;
};

class PolyExp{
public:
    PolyExp(Polynomial pref, Polynomial expo);

    double operator()(double x) override;
    double derivative(double x, size_t n) override;
    double get_Gn(double x, size_t n);

private:
    Polynomial pref;
    Polynomial expo;

    std::vector<double> dg;
};