/*
Author: Vegard Gjeldvik Jervell
Contains: 
        - Classes to facilitate analytic simplification of fractions containing factorials to prevent overflow.
        - Classes to compute analytical derivatives of functions of the form t(x) = f(x) exp[g(x)]

    -------------------------------------------    Factorials    ------------------------------------------------------
        The Fac, Product and Frac classes represent analytic factorials, products and fractions. They overload
        the relevant arithmetic operators (* and /), and are simplified analytically before numerical evaluation when
        evaluated with <expression.eval()>

        Frac can be initialized from both Fac and Product, as a fraction with denominator = 1. Thus, all arithmetic
        is only defined for the Frac class, and when eg. computing Fac * Product, the result is implicitly converted
        to a fraction.

        The function 'ipow(int, int) is also defined. This function is an exponential that returns a Product (not
        evaluated to an int), to ensure that Fractions containing powers can be simplified before evaluation.

        Working principle (Factorials):
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

    -------------------------------------------    Derivatives    ------------------------------------------------------
        The Polynomial and PolyExp classes represent polynomials with arbitrary (both positive and negative) integer
        exponents and "t-functions" respectively, where "t-functions" refer to functions of the form t(x) = f(x) exp[g(x)],
        with f(x) and g(x) polynomials.

        Both the Polynomial and PolyExp class inherit from the "Term" abstract class, such that one can create 
        pointers to "Term" in order to hold combinations of Polynomials and t-functions. However, note that the 
        PolyExp class is implemented with specific handling of the case where g(x) = 0, such that there should be
        no performance difference between using the objects

        Polynomial f(<params>); // Just a polynomial
        PolyExp t(f, Polynomial::zero()); // Exactly the same as f (because g(x) = 0).

        The functions get_partition(int, int) and get_partition_multiplicity(const std::vector<int>&) are defined,
        the first finds all integer partitions of a given number N, with largest element smaller than or equal to m.
        The second computes the "multiplicity" of a given partition (see derivative memo for details).
*/
#pragma once
#include "KineticGas.h"
#include "global_params.h"
#include "math.h"
#include "cmath"
#include <array>
#include <cstdio>
#include <iostream>

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

inline double partialfactorial(int start, int stop){
    // Evaluate (stop! / start!)
    double fac = 1.;
    for (; stop > start; stop--){
        fac *= stop;
    }
    return fac;
}

double binom(int n, int k);

std::vector<std::vector<int>> get_partitions(int N, int m);
std::vector<std::vector<std::vector<int>>> build_partitions(int N, int maxval=-1);
double partition_multiplicity(const std::vector<int>& partition);

class Term{
public:
    virtual double operator()(double x) const {
        return derivative(x, 0);
    }

    virtual double derivative(double x, int n) const = 0;
};

class PolyExp;
class Polynomial : public Term {
public:
    friend class PolyExp;

    /*
        Represents a polynomial 
        f(x) = \sum_{k = k_{min}}^{k_{max}} C_k x^k,
        where k_min and k_max can be any integer.
    */
    Polynomial(int k_min, int k_max, std::vector<double> coeff, int k_step=1);
    Polynomial(double val);
    static Polynomial constant(double v) {return Polynomial(0, 0, {v});}
    static Polynomial zero() {return Polynomial::constant(0);}
    static Polynomial unity() {return Polynomial::constant(1);} 
    static Polynomial linear(double a, double b) {return Polynomial(0, 1, {a, b});}
    static Polynomial linear(double a) {return Polynomial(1, 1, {a});}

    double operator()(double x) const override;
    double derivative(double x, int n) const override;
    inline bool is_linear() const {return _is_linear;}
    inline bool is_constant() const {return _is_constant;}

    friend std::ostream& operator<<(std::ostream& strm, const Polynomial& p);
// private:
    const int k_min, k_max, k_step;
    const bool _is_linear, _is_constant;
    const std::vector<double> coeff;
};

class PolyExp : public Term{
public:
    /*
        Represents the function t(x) = f(x) exp[g(x)], where
        f is "pref" and g is "expo".
    */
    PolyExp(Polynomial pref, Polynomial expo);
    PolyExp(Polynomial pref, double expo);
    PolyExp(double pref, Polynomial expo);

    double operator()(double x) const override;
    double derivative(double x, int n) const override;
    double get_Gk(double x, int k, const vector1d& dg, const std::vector<std::vector<int>>& partitions) const;
    double get_Gk(double x, int k, const vector1d& dg) const; // Less efficient version, used if partitions are not known beforehand.

    friend std::ostream& operator<<(std::ostream& strm, const PolyExp& p);
// private:
    const Polynomial pref;
    const Polynomial expo;
};