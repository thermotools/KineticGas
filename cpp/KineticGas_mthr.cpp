/*
Author: Vegard Gjeldvik Jervell
Contains: Methods to precompute collision integrals on multiple threads.
            The idea is that once the Enskog approximation order is known, the collision integrals can be computed
            on multiple threads and stored in the omega_map, before computations continue on a single thread.
Usage: Set the constexpr int Ncores to the desired number of threads to split the computations among. Note that
        the maximum number of computations on a thread is the limiting factor for speed, so if 24 computations are
        required, using 10 cores will not be more advantageous than 8 cores, while 12 will give an increase in speed.
*/

#include "KineticGas.h"

constexpr size_t Ncores{8}; // Number of cores to use for multithreading collision integrals (save one for the transfer lengths)

// Function to distribute the nvals values of vec as evenly as possible among Ncores different vectors
inline std::vector<std::vector<int>> slice_by_ncores(const std::vector<int>& vec, int n_vals){
    std::vector<std::vector<int>> slices(Ncores);
    int vals_per_thread = n_vals / Ncores; // Minimum number of values computed on a thread
    int val_idx{0};
    for (int core_idx = 0; core_idx < Ncores; core_idx++){
        for (;val_idx < vals_per_thread * (core_idx + 1); val_idx++){
            slices[core_idx].push_back(vec[val_idx]);
        }
    }
    int core_idx{0};
    for (; val_idx < n_vals; val_idx++){
        slices[core_idx].push_back(vec[val_idx]);
        core_idx++;
    }
    return slices;
}

// Computes the omega integral for each (i, j, l, r) in the vectors
// If i = [1, 2], j = [3, 4], l = [5, 6], r = [7, 8], the method will compute the omega integrals
// (1, 3, 5, 7) and (2, 4, 6, 8)
// This is the method that is activated on different threads simultaneously in order to parrellelise the evaluation
// of omega integrals. Ensuring that no integrals are computed twice must be done beforehand by not passing the same
// combination of (i, j, l, r) to different threads. If the same combination is passed twice, this may lead
// to undefined behaviour if two threads simultaneously try to write to omega_map.
void KineticGas::precompute_omega(const std::vector<int>& i_vec, const std::vector<int>& j_vec,
                        const std::vector<int>& l_vec, const std::vector<int>& r_vec, double T)
     {
        for (int idx = 0; idx < i_vec.size(); idx++){
            omega(i_vec[idx], j_vec[idx], l_vec[idx], r_vec[idx], T); // Computed values are stored in omega_map
        }
     }

/*
 The following methods identify which omega integrals are required for the computation of a given property at
 a given Enskog approximation order, then distribute those properties among different threads that run
 the 'precompute_omega' method. The distribution is done using the 'slice_by_ncores' function.

 First, generate vectors of all the (i, j, l, r) values for which the omega integrals must be computed,
 then, distribute the values among Ncores slices, where the index of the values still matches after slicing, such
 that i = [1, 2, 3, 4], j = [5, 6, 7, 8] can be sliced to i_slices = [[1, 2], [3, 4]], j_slices = [[5, 6], [7, 8]]
 Finally, pass the slices to individual threads that compute omega_integrals for their slice.
 */
void KineticGas::precompute_conductivity(int N, double T, double rho, bool precompute_etl){
    int true_Ncomps = (is_singlecomp) ? 1 : Ncomps;
    std::vector<int> i_vec, j_vec, l_vec, r_vec;
    int n_vals{0};
    for (int i = 0; i < true_Ncomps; i++){
        for (int j = i; j < true_Ncomps; j++){
            for (int l = 1; l <= N; l++){
                for (int r = l; r <= 2 * (N - 1) + 2 - l; r++){
                    i_vec.push_back(i);
                    j_vec.push_back(j);
                    l_vec.push_back(l);
                    r_vec.push_back(r);
                    n_vals++;
                }
            }
        }
    }
    std::vector<std::vector<int>> i_slices{slice_by_ncores(i_vec, n_vals)},
                                  j_slices{slice_by_ncores(j_vec, n_vals)},
                                  l_slices{slice_by_ncores(l_vec, n_vals)},
                                  r_slices{slice_by_ncores(r_vec, n_vals)};

    std::vector<std::thread> threads;
    for (int core_idx = 0; core_idx < Ncores; core_idx++){
        if (i_slices[core_idx].size() == 0) continue;
        threads.push_back(std::thread(&KineticGas::precompute_omega, this, std::ref(i_slices[core_idx]), std::ref(j_slices[core_idx]),
                           std::ref(l_slices[core_idx]), std::ref(r_slices[core_idx]), std::ref(T)));
    }
    if (precompute_etl && !is_idealgas) get_etl(rho, T, vector1d(Ncomps, 0));
    for (int core_idx = 0; core_idx < Ncores; core_idx++){
        if (i_slices[core_idx].size() == 0) continue;
        threads[core_idx].join();
    }
}
void KineticGas::precompute_diffusion(int N, double T, double rho){
    precompute_conductivity(N, T, rho, false); // Don't need ETL for diffusion
}
void KineticGas::precompute_th_diffusion(int N, double T, double rho){
    precompute_conductivity(N, T, rho, true);
}
void KineticGas::precompute_viscosity(int N, double T, double rho){
    int true_Ncomps = (is_singlecomp) ? 1 : Ncomps;
    std::vector<int> i_vec, j_vec, l_vec, r_vec;
    int n_vals{0};
    for (int i = 0; i < true_Ncomps; i++){
        for (int j = i; j < true_Ncomps; j++){
            for (int l = 1; l <= N + 1; l++){
                if ((l % 2 == 1) && (is_singlecomp)) continue;
                for (int r = l; r <= 2 * (N - 1) + 4 - l; r++){
                    i_vec.push_back(i);
                    j_vec.push_back(j);
                    l_vec.push_back(l);
                    r_vec.push_back(r);
                    n_vals++;
                }
            }
        }
    }
    std::vector<std::vector<int>> i_slices{slice_by_ncores(i_vec, n_vals)},
                                  j_slices{slice_by_ncores(j_vec, n_vals)},
                                  l_slices{slice_by_ncores(l_vec, n_vals)},
                                  r_slices{slice_by_ncores(r_vec, n_vals)};

    std::vector<std::thread> threads;
    for (int core_idx = 0; core_idx < Ncores; core_idx++){
        if (i_slices[core_idx].size() == 0) continue;
        threads.push_back(std::thread(&KineticGas::precompute_omega, this, std::ref(i_slices[core_idx]), std::ref(j_slices[core_idx]),
                           std::ref(l_slices[core_idx]), std::ref(r_slices[core_idx]), std::ref(T)));
    }
    if (!is_idealgas) get_mtl(rho, T, vector1d(Ncomps, 0));
    for (int core_idx = 0; core_idx < Ncores; core_idx++){
        if (i_slices[core_idx].size() == 0) continue;
        threads[core_idx].join();
    }
}

void KineticGas::precompute_ideal_diffusion(vector1d T){
    double rho = AVOGADRO;
    int N = 1;
    std::vector<std::thread> threads;
    size_t ti = 0;
    size_t tj = 0;
    do {
        for (; ti < std::max(Ncores, T.size()); ti++){
            threads.push_back(std::thread(&KineticGas::precompute_diffusion, this, std::ref(N), std::ref(T[ti]), std::ref(rho)));
        }
        for (; tj < std::max(Ncores, T.size()); tj++){
            threads[tj].join();
        }
    } while (ti < T.size());

}