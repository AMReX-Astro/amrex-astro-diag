#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>

#include <extern_parameters.H>

#include <network.H>
#include <eos.H>

#include <amrex_astro_util.H>

using namespace amrex;

void main_main()
{

    std::string pltfile(diag_rp::plotfile);

    if (pltfile.empty()) {
        std::cout << "no plotfile specified" << std::endl;
        std::cout << "use: diag.plotfile=plt00000 (for example)" << std::endl;
        amrex::Error("no plotfile");
    }

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    std::string outfile = "conduction." +
        std::filesystem::path(pltfile).filename().string();


    PlotFileData pf(pltfile);

    const int ndims = pf.spaceDim();
    AMREX_ALWAYS_ASSERT(ndims <= AMREX_SPACEDIM);

    const int nlevs = pf.finestLevel() + 1;

    Vector<std::string> varnames;
    varnames = pf.varNames();

    // find variable indices -- we want density, temperature, and species.
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    // the plotfile can store either (rho X) or just X alone.  Here we'll assume
    // that we have just X alone

    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = get_dens_index(var_names_pf);
    int energy_comp = get_energy_index(var_names_pf);
    int pres_comp = get_pres_index(var_names_pf);
    int ux_comp = get_ux_index(var_names_pf);
    int uy_comp = get_uy_index(var_names_pf);
    int diff_comp = get_diff_index(var_names_pf);


    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;
    gvarnames.push_back("enthalpy");
    gvarnames.push_back("enthalpy_flux");
    gvarnames.push_back("conduction_flux");

    // interpret the boundary conditions

    BCRec bcr_default;
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    IntVect ng(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim < ndims) {
            bcr_default.setLo(idim, BCType::hoextrapcc);
            bcr_default.setHi(idim, BCType::hoextrapcc);
        } else {
            bcr_default.setLo(idim, BCType::int_dir);
            bcr_default.setHi(idim, BCType::int_dir);
            is_periodic[idim] = 1;
            ng[idim] = 0;
        }
    }

    // get center if spherical

    Array<Real, AMREX_SPACEDIM> center;
    auto const probLo = pf.probLo();
    auto const probHi = pf.probHi();

    // we need both rho, T, P, u_x and u_y constructed with ghost cells

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), static_cast<int>(gvarnames.size()), 0);

        Vector<BCRec> bcr{bcr_default};
        auto is_per = is_periodic;

        Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                         pf.coordSys(), is_per);
        geom.push_back(vargeom);

        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
            (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

        // fill the density and temperature mfs with ghost cells
        // we also need all of the species


        MultiFab    dens_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab  energy_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab    pres_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab      ux_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab      uy_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab    diff_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);

        if (ilev == 0) {

            // density
            {
                MultiFab smf = pf.get(ilev, var_names_pf[dens_comp]);
                FillPatchSingleLevel(dens_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // energy
            {
                MultiFab smf = pf.get(ilev, var_names_pf[energy_comp]);
                FillPatchSingleLevel(energy_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // pressure
            {
                MultiFab smf = pf.get(ilev, var_names_pf[pres_comp]);
                FillPatchSingleLevel(pres_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // u_x
            {
                MultiFab smf = pf.get(ilev, var_names_pf[ux_comp]);
                FillPatchSingleLevel(ux_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // u_y
            {
                MultiFab smf = pf.get(ilev, var_names_pf[uy_comp]);
                FillPatchSingleLevel(uy_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // diff term
            {
                MultiFab smf = pf.get(ilev, var_names_pf[diff_comp]);
                FillPatchSingleLevel(diff_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                        0, 0, 1, vargeom, physbcf, 0);
            }

        } else {
            auto* mapper = (Interpolater*)(&cell_cons_interp);

            IntVect ratio(pf.refRatio(ilev-1));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }

            Geometry cgeom(pf.probDomain(ilev-1), RealBox(pf.probLo(),pf.probHi()),
                           pf.coordSys(), is_per);
            PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cphysbcf
                (cgeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

            // density
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[dens_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[dens_comp]);
                FillPatchTwoLevels(dens_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                    {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                    cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // energy
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[energy_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[energy_comp]);
                FillPatchTwoLevels(energy_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // pressure
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[pres_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[pres_comp]);
                FillPatchTwoLevels(pres_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // ux
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[ux_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[ux_comp]);
                FillPatchTwoLevels(ux_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                    {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                    cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // uy
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[uy_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[uy_comp]);
                FillPatchTwoLevels(uy_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                    {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                    cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // diff term
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[diff_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[diff_comp]);
                FillPatchTwoLevels(diff_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                    {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                    cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

        }

        auto const& dx = pf.cellSize(ilev);

        const MultiFab& lev_data_mf = pf.get(ilev);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(dens_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();

            // output storage
            auto const& out = gmf[ilev].array(mfi);

            // temperature and pressure with ghost cells
            auto const& rho  = dens_mf.const_array(mfi);
            auto const& rhoE    = energy_mf.const_array(mfi);
            auto const& p    = pres_mf.const_array(mfi);
            auto const& ux   = ux_mf.const_array(mfi);
            auto const& uy   = uy_mf.const_array(mfi);
            auto const& diff = diff_mf.const_array(mfi);

            // all of the data without ghost cells, needed for eos calculations
            // const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                if (ndims == 2) {
                    out(i,j,k,0) = rhoE(i,j,k)/rho(i,j,k) + p(i,j,k) / rho(i,j,k);

                    Real h_flux_hi = (rho(i+1,j,k) * ux(i+1,j,k) * out(i+1,j,k,0) - rho(i,j,k) * ux(i,j,k) * out(i,j,k,0))/dx[0] \
                                      + (rho(i,j+1,k) * uy(i,j+1,k) * out(i,j+1,k,0) - rho(i,j,k) * uy(i,j,k) * out(i,j,k,0))/dx[1];

                    Real h_flux_lo = (rho(i,j,k) * ux(i,j,k) * out(i,j,k,0) - rho(i-1,j,k) * ux(i-1,j,k) * out(i-1,j,k,0))/dx[0] \
                                      + (rho(i,j,k) * uy(i,j,k) * out(i,j,k,0) - rho(i,j-1,k) * uy(i,j-1,k) * out(i,j-1,k,0))/dx[1];

                    out(i,j,k,1) = 0.5_rt * (h_flux_hi + h_flux_lo);

                    out(i,j,k,2) = diff(i,j,k);
                } else {
                    amrex::Error("This routine is only used in 2D");
                }
            });
        }
    }

    Vector<int> level_steps;
    Vector<IntVect> ref_ratio;
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        level_steps.push_back(pf.levelStep(ilev));
        if (ilev < pf.finestLevel()) {
            ref_ratio.push_back(IntVect(pf.refRatio(ilev)));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ref_ratio[ilev][idim] = 1;
            }
        }
    }

    WriteMultiLevelPlotfile(outfile, nlevs, GetVecOfConstPtrs(gmf), gvarnames,
                            geom, pf.time(), level_steps, ref_ratio);
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv);

    // initialize the runtime parameters

    init_extern_parameters();

    // initialize C++ Microphysics

    eos_init(diag_rp::small_temp, diag_rp::small_dens);
    network_init();

    main_main();
    amrex::Finalize();
}
