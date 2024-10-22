/*
Various utility structures and enums

References:
    (I) The Kinetic Gas theory of Mie fluids, V. G. Jervell, Norwegian university of science and technology (2022)
    (II) Revised Enskog theory for Mie fluids: Prediction of diffusion coefficients, thermal diffusion coefficients,
    viscosities and thermal conductivities, V. G. Jervell and Ø. Wilhelmsen, J. Chem. Phys. (2023)
*/

/*
   To avoid unneccesary evaluations of the collision integrals, this struct is used to represent a point in 
   The five-dimensional (i, j, l, r, T)-space where the collision integral has been evaluated.

   NB: Resolution along the T-axis is 0.1 K, as (by experience) the collision integrals are a weak enough
   function of T to justify using T-values rounded to the nearest .1 K to improve speed, with little cost to precision.
*/
struct StatePoint{
    int T_dK;
    double rho;
    StatePoint(double T) : T_dK{static_cast<int>((T * 100.) + 0.5)} {}
    StatePoint(double T, double rho) : T_dK{static_cast<int>((T * 100.) + 0.5)}, rho{rho}{}

    bool operator<(const StatePoint& other) const {
        if (T_dK < other.T_dK) return true;
        else if (T_dK == other.T_dK){
            if (rho < other.rho) return true;
        }
        return false;
    }
};

struct OmegaPoint{
    int i, j, l, r, T_dK;
    double rho;
    OmegaPoint(int i, int j, int l, int r, double T, double rho) : i{i}, j{j}, l{l}, r{r}, rho{rho} {
         T_dK = (int) ((T * 10.0) + 0.5); // Temperature in dK (10^-1 K)
    };

    OmegaPoint(int i, int j, int l, int r, double T) : OmegaPoint(i, j, l, r, T, 0){}

    bool operator<(const OmegaPoint& other) const {
        if (i < other.i) return true;
        else if (i == other.i){
            if (j < other.j) return true;
            else if (j == other.j){
                if (l < other.l) return true;
                else if (l == other.l){
                    if (r < other.r) return true;
                    else if (r == other.r){
                        if (T_dK < other.T_dK) return true;
                        else if (T_dK == other.T_dK){
                            if (rho < other.rho) return true;
                        }
                    }
                }
            }
        }
        return false;
    }
};

enum FrameOfReference{
    CoM,
    CoN,
    CoV,
    solvent,
    zarate,
    zarate_x,
    zarate_w
};

enum TransferLengthModel{
    DEFAULT = -1, // Default model
    collision_diameter, // Model presented in Refs. (I) and (II) (see top of file)
    EWCA, // Unpublished model for momentum- and energy transfer lengths (MTL and ETL)
    correlation, // Unpublished correlation for Argon transfer lengths
    INVALID // Used to terminate loops over all model identifiers
};