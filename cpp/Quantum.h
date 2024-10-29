#include "KineticGas.h"

class Quantum : public KineticGas {
public:
    Quantum(std::string comps) : KineticGas(comps, false) {
        half_spin = std::vector<size_t>(Ncomps, 0);
        for (size_t i = 0; i < Ncomps; i++){
            half_spin[i] = static_cast<size_t>(static_cast<double>(compdata[i]["spin"]) * 2 + 0.5);
        }
    }

    double potential(int i, int j, double r);
    double potential_r(int i, int j, double r);
    double potential_rr(int i, int j, double r);

    vector2d wave_function(int i, int j, int l, double E, double r_end, double dr);
    double phase_shift(int i, int j, int l, double E);
    double cross_section_A(int n, int l, size_t k);
    double cross_section_kernel(int i, int j, double n, double l, double E);
    double cross_section(int i, int j, int n, double E);
    double omega(int i, int j, int n, int s, double T);

    std::vector<size_t> half_spin; // Spin of each particle multiplied by two 
protected:
    vector2d model_rdf(double rho, double T, const vector1d& mole_fracs) override {throw std::runtime_error("Method model_rdf not implemented for Quantum!");}
    vector2d model_mtl(double rho, double T, const vector1d& x) override {throw std::runtime_error("Method model_mtl not implemented for Quantum!");}
    vector2d model_etl(double rho, double T, const vector1d& x) override {throw std::runtime_error("Method model_etl not implemented for Quantum!");}

private:
    

};