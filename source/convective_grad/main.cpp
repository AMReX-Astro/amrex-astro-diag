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

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile(diag_rp::plotfile);

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    std::string outfile = "convgrad." +
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

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "density");
    if (idx == var_names_pf.cend()) {
        // MAESTROeX uses "rho"
        idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "rho");
        if (idx == var_names_pf.cend()) {
            amrex::Error("Error: could not find the density component");
        }
    }
    int dens_comp = std::distance(var_names_pf.cbegin(), idx);

    idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "Temp");
    if (idx == var_names_pf.cend()) {
        // for MAESTROeX, we'll default to "tfromp"
        idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "tfromp");
        if (idx == var_names_pf.cend()) {
            amrex::Error("Error: could not find the temperature component");
        }
    }
    int temp_comp = std::distance(var_names_pf.cbegin(), idx);

    idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "pressure");
    if (idx == var_names_pf.cend()) {
        // for MAESTROeX, we'll default to "tfromp"
        idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "p0");
        if (idx == var_names_pf.cend()) {
            amrex::Error("Error: could not find the pressure component");
        }
    }
    int pres_comp = std::distance(var_names_pf.cbegin(), idx);

    std::string first_spec_name = "X(" + short_spec_names_cxx[0] + ")";
    idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), first_spec_name);
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find the first species");
    }
    int spec_comp = std::distance(var_names_pf.cbegin(), idx);

    // safety check -- make sure the species in the plotfile are identical to
    // those defined in the network we built this tool with.

    for (int n = 0; n < NumSpec; ++n) {
        std::string current_spec_name = "X(" + short_spec_names_cxx[n] + ")";
        if (current_spec_name != var_names_pf[spec_comp+n]) {
            std::cout << current_spec_name << std::endl;
            std::cout << var_names_pf[spec_comp+n] << std::endl;
            amrex::Error("Error: species don't match");
        }
    }

    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;
    gvarnames.push_back("del");
    gvarnames.push_back("del_ad");

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

    // we need both T and P constructed with ghost cells

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), gvarnames.size(), 0);

        Vector<BCRec> bcr{bcr_default};
        auto is_per = is_periodic;

        Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                         pf.coordSys(), is_per);
        geom.push_back(vargeom);

        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
            (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

        // fill the pressure and temperature mfs with ghost cells

        MultiFab temp_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab pres_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);

        if (ilev == 0) {

            // temperature
            {
                MultiFab smf = pf.get(ilev, var_names_pf[temp_comp]);
                FillPatchSingleLevel(temp_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // pressure
            {
                MultiFab smf = pf.get(ilev, var_names_pf[pres_comp]);
                FillPatchSingleLevel(pres_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
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

            // temperature
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[temp_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[temp_comp]);
                FillPatchTwoLevels(temp_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // pressure
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[temp_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[temp_comp]);
                FillPatchTwoLevels(pres_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

        }

        auto const& dx = pf.cellSize(ilev);

        const MultiFab& lev_data_mf = pf.get(ilev);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(temp_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();

            // output storage
            auto const& ga = gmf[ilev].array(mfi);

            // temperature and pressure with ghost cells
            auto const& T = temp_mf.const_array(mfi);
            auto const& P = pres_mf.const_array(mfi);

            // all of the data without ghost cells
            const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                if (ndims == 1) {
                    // x is the vertical
                    ga(i,j,k,0) = (T(i+1,j,k) - T(i-1,j,k)) /
                                  (P(i+1,j,k) - P(i-1,j,k)) *
                                  (P(i,j,k) / T(i,j,k));

                } else if (ndims == 2) {
                    // y is the vertical
                    ga(i,j,k,0) = (T(i,j+1,k) - T(i,j-1,k)) /
                                  (P(i,j+1,k) - P(i,j-1,k)) *
                                  (P(i,j,k) / T(i,j,k));

                } else {
                    // z is the vertical
                    ga(i,j,k,0) = (T(i,j,k+1) - T(i,j,k-1)) /
                                  (P(i,j,k+1) - P(i,j,k-1)) *
                                  (P(i,j,k) / T(i,j,k));
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