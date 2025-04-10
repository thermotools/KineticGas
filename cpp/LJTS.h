/*Author: Johannes Salomonsen Løken
Puropse: Implementation of Lennard-Jones truncated and shifted, with r_cut = 2.5*sigma
only for dilute gas transport properties.
*/

#pragma once
#include "Spherical.h"
#include "LJSpline.h"

class LJTS : public Spherical {
    public:
    double rc; 
    
    LJTS(std::vector<double> mole_weights, std::vector<std::vector<double>> sigmaij, std::vector<std::vector<double>> eps, bool is_idealgas, bool is_singlecomp)
    : Spherical(mole_weights, sigmaij, eps, true, true), rc{2.5*sigmaij[0][0]} 
    {
        // lag eos//LJs_bh bh_eos{"Default",1.0};
        // GenericEoS ljs_eos{ThermoWrapper(std::move(bh_eos))};
        //this -> set_eos(std::move(ljs_eos));
        // eos = std::make_unique<GenericEoS>(Ideal());
        if ((sigmaij.size() > 2) | (mole_weights.size() > 2) | (eps.size() > 2)) 
        {
            throw std::invalid_argument("The Lennard-Jones/spline is not implemented for multicomponent systems (yet)!");
        }
    }   
    dual2 LJ(int i, int j, dual2 r) {
        return 4 * eps[0][0] * (pow(sigma[0][0] / r, 12) - pow(sigma[0][0] / r, 6));
    }

    dual2 potential(int i, int j, dual2 r) override {
        if (r <= rc) {
            return LJ(i,j,r) - LJ(i,j, static_cast<dual2>(rc));
        }
        else {
            return static_cast<dual2>(0);
        }
    }

    double potential(int i, int j, double r) override{
        return potential(i, j, static_cast<dual2>(r)).val.val;
    }

    double potential_derivative_r(int i, int j, double r) override {
        if (r <= rc) {
            return -24*eps[0][0]/r*(2*pow(sigma[0][0] / r, 12) - pow(sigma[0][0] / r, 6));
        }
        else {
            return 0;
        }
    }
    double potential_dblderivative_rr(int i, int j, double r) override {
        if (r <= rc) {
            return 24*eps[0][0]/pow(r,2)*(26*pow(sigma[0][0] / r, 12) - 7*pow(sigma[0][0] / r, 6));
        }
        else {
            return 0;
        }
    }
    inline std::vector<std::vector<double>> model_rdf(double rho, double T, const std::vector<double>& x) override {
        return {{1.,1.},{1.,1.}};
    }
};