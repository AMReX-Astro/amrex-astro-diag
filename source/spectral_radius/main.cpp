//
// loop over zones and compute the spectral radius of the reaction
// network Jacobian and then output the thermodynamic state of the
// zone with the largest spectral radius.
//
#include <iostream>
// #include <stringstream>
#include <regex>
#include <string>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <amrex_astro_util.H>

#include <extern_parameters.H>

#include <network.H>
#include <eos.H>
#include <burn_type.H>
#include <rkc_type.H>
#include <integrator_data.H>
#include <circle_theorem.H>

using namespace amrex;


int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv);

    // initialize the runtime parameters

    init_extern_parameters();

    // initialize C++ Microphysics

    eos_init(diag_rp::small_temp, diag_rp::small_dens);
    network_init();

    // timer for profiling

    BL_PROFILE_VAR("main()", pmain);

    // read the plotfile metadata

    PlotFileData pf(diag_rp::plotfile);

    int fine_level = pf.finestLevel();
    const int dim = pf.spaceDim();

    // find variable indices -- we want density, temperature, and species.
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    // the plotfile can store either (rho X) or just X alone.  Here we'll assume
    // that we have just X alone

    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = get_dens_index(var_names_pf);
    int temp_comp = get_temp_index(var_names_pf);
    int spec_comp = get_spec_index(var_names_pf);

    amrex::Real sprad_max{-1.0};
    burn_t burn_state_max;

    for (int ilev = 0; ilev <= fine_level; ++ilev) {

        IntVect ratio{pf.refRatio(ilev)};
        for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
            ratio[idim] = 1;
        }

        const MultiFab& lev_data_mf = pf.get(ilev);

        for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.ok()) {
                const auto& fab = lev_data_mf.array(mfi);
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                for (int k = lo.z; k <= hi.z; ++k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {

                            burn_t burn_state;

                            burn_state.rho = fab(i,j,k,dens_comp);
                            burn_state.T = fab(i,j,k,temp_comp);
                            for (int n = 0; n < NumSpec; ++n) {
                                burn_state.xn[n] = fab(i,j,k,spec_comp+n);
                            }

                            // we also need internal energy
                            eos(eos_input_rt, burn_state);
                            burn_state.e_scale = burn_state.e;

                            // we need to load up the integrator type
                            constexpr int int_neqs = integrator_neqs<burn_t>();
                            rkc_t<int_neqs> rstate;
                            burn_to_integrator(burn_state, rstate);
                            rstate.t = 0.0;

                            amrex::Real sprad{};
                            circle_theorem_sprad(rstate.t, burn_state, rstate, sprad);

                            if (std::abs(sprad) > sprad_max) {
                                sprad_max = std::abs(sprad);
                                burn_state_max = burn_state;
                            }
                        }
                    }
                }

            } // bx.ok()

        } // MFIter

    } // level loop

    std::cout << "maximum spectral radius = " << sprad_max << std::endl;
    std::cout << "T = " << burn_state_max.T << std::endl;
    std::cout << "rho = " << burn_state_max.rho << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "X(" << short_spec_names_cxx[n] << ") = " << burn_state_max.xn[n] << std::endl;
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}

