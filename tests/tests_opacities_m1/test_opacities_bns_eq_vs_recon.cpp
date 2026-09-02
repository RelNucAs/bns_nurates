//=================================================
// bns-nurates neutrino opacities code
// ================================================
//! \file  test_opacities_bns_eq_vs_recon.cpp
//  \brief Reproduce the Chiesa+25 (arXiv:2412.04570) Fig. 8 comparison along a
//         1-D BNS profile: for each radial zone, compute the grey (energy-
//         integrated) neutrino opacities with (a) the neutrino distribution
//         RECONSTRUCTED from the M1 moments (n, J, chi) and (b) the EQUILIBRIUM
//         (Fermi-Dirac) distribution, and print them side by side.
//
//  Output per zone, per species [nue, anue, nux, anux]:
//     eta_0     : energy-integrated NUMBER emissivity      [cm^-3 s^-1]
//     eta       : energy-integrated ENERGY emissivity      [MeV cm^-3 s^-1]
//     kappa_0_a : energy-averaged STIMULATED opacity (number) [cm^-1]
//     kappa_a   : energy-averaged STIMULATED opacity (energy) [cm^-1]
//  ...once for reconstructed ("_re") and once for equilibrium ("_eq").
//
//  NOTE on heavy-lepton lumping: bns_nurates treats id_nux/id_anux as a SINGLE
//  heavy-lepton species. The production runs multiply the nux/anux emissivities
//  by 2 (mu + tau) in the AthenaK wrapper; that factor is NOT applied here so
//  the eq-vs-recon comparison is at the raw library level (the x2 cancels in the
//  ratio anyway). Absorption is per-neutrino and never gets the x2.
//
//  Usage:  ./test_opacities_bns_eq_vs_recon [profile.txt]
//  Default profile: inputs/BNS/1d/DD2_M12980-12980_M1_LR_t_50ms_x_axis.txt

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../tests.hpp"   // pulls in Kokkos, bns_nurates.hpp, distribution.hpp,
                          // integration.hpp, m1_opacities.hpp and the Dev*/Host
                          // typedefs used below.

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
static constexpr int   MAX_ZONES = 4096; // generous upper bound for the profile
static constexpr int   N_QUAD    = 10;   // Gauss-Legendre points (2*N_QUAD<=BS_N_MAX)

// Unit conversions for OUTPUT only (library works in nm; convert to cgs):
//   number emissivity  nm^-3 s^-1     -> cm^-3 s^-1      : x 1e21
//   energy emissivity  MeV nm^-3 s^-1 -> MeV cm^-3 s^-1  : x 1e21
//   opacity (kappa)    nm^-1          -> cm^-1           : x 1e7
static constexpr double C_ETA = 1e21;
static constexpr double C_KAP = 1e7;

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
        const char* profile =
            (argc > 1) ? argv[1]
                       : "inputs/BNS/1d/DD2_M12980-12980_M1_LR_t_50ms_x_axis.txt";
        // Optional 2nd arg: "neps" -> inelastic scattering only; "iso" ->
        // iso-energetic scattering only. Default: all reactions.
        const bool neps_only = (argc > 2 && std::strcmp(argv[2], "neps") == 0);
        const bool iso_only  = (argc > 2 && std::strcmp(argv[2], "iso") == 0);

        // -------------------------------------------------------------------
        // Read the profile on the host (skip lines beginning with '#').
        // Columns (1-indexed) used: 2:r 3:rho 4:T 5:Ye 6:Yn 7:Yp
        //   10:mu_e 11:mu_n 12:mu_p 17:Ynue 18:Yanue 19:Ynux
        //   20:Enue 21:Eanue 22:Enux 23:chinue 24:chianue 25:chinux
        // -------------------------------------------------------------------
        FILE* fp = fopen(profile, "r");
        if (!fp)
        {
            fprintf(stderr, "ERROR: cannot open profile '%s'\n", profile);
            Kokkos::finalize();
            return 1;
        }

        Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_r("r", MAX_ZONES),
            h_rho("rho", MAX_ZONES), h_T("T", MAX_ZONES), h_Ye("Ye", MAX_ZONES),
            h_Yn("Yn", MAX_ZONES), h_Yp("Yp", MAX_ZONES),
            h_mue("mue", MAX_ZONES), h_mun("mun", MAX_ZONES),
            h_mup("mup", MAX_ZONES);
        // per-species (nue, anue, nux) number fraction, energy/baryon, chi
        Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_Ynu("Ynu", MAX_ZONES, 3),
            h_Enu("Enu", MAX_ZONES, 3), h_chi("chi", MAX_ZONES, 3);

        char line[2048];
        int nz = 0;
        while (fgets(line, sizeof(line), fp))
        {
            // skip comments / blank lines
            const char* p = line;
            while (*p == ' ' || *p == '\t') ++p;
            if (*p == '#' || *p == '\n' || *p == '\0') continue;

            double c[25];
            int got = std::sscanf(
                line,
                "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8],
                &c[9], &c[10], &c[11], &c[12], &c[13], &c[14], &c[15], &c[16],
                &c[17], &c[18], &c[19], &c[20], &c[21], &c[22], &c[23], &c[24]);
            if (got < 25) continue;
            if (nz >= MAX_ZONES) { fprintf(stderr, "too many zones\n"); break; }

            h_r(nz)   = c[1];   // radius [km]
            h_rho(nz) = c[2];   // matter density [g cm^-3]
            h_T(nz)   = c[3];   // temperature [MeV]
            h_Ye(nz)  = c[4];
            h_Yn(nz)  = c[5];
            h_Yp(nz)  = c[6];
            h_mue(nz) = c[9];   // mu_e [MeV] (rest mass included)
            h_mun(nz) = c[10];  // mu_n [MeV]
            h_mup(nz) = c[11];  // mu_p [MeV]
            h_Ynu(nz, 0) = c[16]; h_Ynu(nz, 1) = c[17]; h_Ynu(nz, 2) = c[18];
            h_Enu(nz, 0) = c[19]; h_Enu(nz, 1) = c[20]; h_Enu(nz, 2) = c[21];
            h_chi(nz, 0) = c[22]; h_chi(nz, 1) = c[23]; h_chi(nz, 2) = c[24];
            ++nz;
        }
        fclose(fp);
        printf("# read %d zones from %s\n", nz, profile);
        if (nz == 0) { Kokkos::finalize(); return 1; }

        // -------------------------------------------------------------------
        // Quadrature (Gauss-Legendre, matches the other opacity tests).
        // -------------------------------------------------------------------
        MyQuadrature my_quad = {.type   = kGauleg,
                                .alpha  = -42.,
                                .dim    = 1,
                                .nx     = N_QUAD,
                                .ny     = 1,
                                .nz     = 1,
                                .x1     = 0.,
                                .x2     = 1.,
                                .y1     = -42.,
                                .y2     = -42.,
                                .z1     = -42.,
                                .z2     = -42.,
                                .points = {0},
                                .w      = {0}};
        if (my_quad.nx * 2 > BS_N_MAX)
        {
            fprintf(stderr, "2*nx exceeds BS_N_MAX\n");
            Kokkos::finalize();
            return 1;
        }
        GaussLegendre(&my_quad);

        // -------------------------------------------------------------------
        // Copy inputs + quadrature to device.
        // -------------------------------------------------------------------
        auto d_r   = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_r);
        auto d_rho = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_rho);
        auto d_T   = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_T);
        auto d_Ye  = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_Ye);
        auto d_Yn  = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_Yn);
        auto d_Yp  = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_Yp);
        auto d_mue = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_mue);
        auto d_mun = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_mun);
        auto d_mup = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_mup);
        auto d_Ynu = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_Ynu);
        auto d_Enu = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_Enu);
        auto d_chi = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_chi);

        Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_nq("nq", 1);
        Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_w("w", my_quad.nx),
            h_pt("pt", my_quad.nx);
        h_nq(0) = my_quad.nx;
        for (int i = 0; i < my_quad.nx; ++i) { h_w(i) = my_quad.w[i]; h_pt(i) = my_quad.points[i]; }
        auto d_nq = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_nq);
        auto d_w  = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_w);
        auto d_pt = Kokkos::create_mirror_view_and_copy(DevMemSpace(), h_pt);

        // Output views: [zone][ {eta0,eta,kap0,kap} x {re,eq} ][species]
        Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace>
            d_eta0_re("eta0_re", nz, 4), d_eta_re("eta_re", nz, 4),
            d_kap0_re("kap0_re", nz, 4), d_kap_re("kap_re", nz, 4),
            d_eta0_eq("eta0_eq", nz, 4), d_eta_eq("eta_eq", nz, 4),
            d_kap0_eq("kap0_eq", nz, 4), d_kap_eq("kap_eq", nz, 4),
            d_kaps_re("kaps_re", nz, 4), d_kaps_eq("kaps_eq", nz, 4),
            d_eta0_eqc("eta0_eqc", nz, 4), d_eta_eqc("eta_eqc", nz, 4);

        // Reaction channels. Default all ON; with "neps" arg, only inelastic
        // scattering (NEPS) is active -- to isolate the NEPS contribution.
        const int fl_absem = (neps_only || iso_only) ? 0 : 1,
                  fl_pair  = (neps_only || iso_only) ? 0 : 1,
                  fl_brem  = (neps_only || iso_only) ? 0 : 1,
                  fl_inel  = iso_only ? 0 : 1,   // NEPS: on for all & neps
                  fl_iso   = neps_only ? 0 : 1;  // iso:  on for all & iso
        // Microphysics corrections: OFF for a clean baseline (they are identical
        // in the eq and recon runs, so they cancel in the comparison). To match
        // the production runs, set these to 1 AND wire eos_pars.dU / dm_eff from
        // the profile's single-particle potentials (cols 13,14) and effective
        // masses (cols 15,16) -- mind the sign convention (dU = U_n - U_p).
        const bool op_dU = false, op_dmeff = false, op_WMab = false,
                   op_WMsc = false, op_decay = false, op_NNmed = false,
                   op_noblock = false;

        // -------------------------------------------------------------------
        // Per-zone computation.
        // -------------------------------------------------------------------
        Kokkos::parallel_for(
            "bns_eq_vs_recon", Kokkos::RangePolicy<>(DevExeSpace(), 0, nz),
            KOKKOS_LAMBDA(const int& i) {
                MyQuadrature quad;
                quad.nx = d_nq(0);
                for (int k = 0; k < quad.nx; ++k) { quad.w[k] = d_w(k); quad.points[k] = d_pt(k); }

                GreyOpacityParams gp;
                gp.opacity_flags = {.use_abs_em          = fl_absem,
                                    .use_pair            = fl_pair,
                                    .use_brem            = fl_brem,
                                    .use_inelastic_scatt = fl_inel,
                                    .use_iso             = fl_iso};
                gp.opacity_pars = {.use_dU              = op_dU,
                                   .use_dm_eff          = op_dmeff,
                                   .use_WM_ab           = op_WMab,
                                   .use_WM_sc           = op_WMsc,
                                   .use_decay           = op_decay,
                                   .brem_implementation = BREM_HR98,
                                   .use_NN_medium_corr  = op_NNmed,
                                   .neglect_blocking    = op_noblock};

                const BS_REAL nb = d_rho(i) / kBS_Mu * 1e-21; // [nm^-3]
                gp.eos_pars.mu_e = d_mue(i);
                gp.eos_pars.mu_p = d_mup(i);
                gp.eos_pars.mu_n = d_mun(i);
                gp.eos_pars.temp = d_T(i);
                gp.eos_pars.yp   = d_Yp(i);
                gp.eos_pars.yn   = d_Yn(i);
                gp.eos_pars.ye   = d_Ye(i);
                gp.eos_pars.nb   = nb;
                gp.eos_pars.dU   = 0.; // unused while op_dU=0 (see note above)

                // ---- (a) RECONSTRUCTED distribution from M1 moments ----
                // per-species number/energy density [nm^-3],[MeV nm^-3]:
                //   n = Y_nu * nb ,  J = E_perbaryon * nb
                // profile carries a single heavy species -> anux := nux
                gp.m1_pars.n[id_nue]   = d_Ynu(i, 0) * nb;
                gp.m1_pars.n[id_anue]  = d_Ynu(i, 1) * nb;
                gp.m1_pars.n[id_nux]   = d_Ynu(i, 2) * nb;
                gp.m1_pars.n[id_anux]  = d_Ynu(i, 2) * nb;
                gp.m1_pars.J[id_nue]   = d_Enu(i, 0) * nb;
                gp.m1_pars.J[id_anue]  = d_Enu(i, 1) * nb;
                gp.m1_pars.J[id_nux]   = d_Enu(i, 2) * nb;
                gp.m1_pars.J[id_anux]  = d_Enu(i, 2) * nb;
                gp.m1_pars.chi[id_nue]  = d_chi(i, 0);
                gp.m1_pars.chi[id_anue] = d_chi(i, 1);
                gp.m1_pars.chi[id_nux]  = d_chi(i, 2);
                gp.m1_pars.chi[id_anux] = d_chi(i, 2);

                // capture the ACTUAL (field) densities before the eq overwrite
                BS_REAL nfield0[4], nfield1[4];
                for (int s = 0; s < 4; ++s)
                { nfield0[s] = gp.m1_pars.n[s]; nfield1[s] = gp.m1_pars.J[s]; }

                gp.distr_pars = CalculateDistrParamsFromM1(&gp.m1_pars, &gp.eos_pars);
                M1Opacities re = ComputeM1Opacities(&quad, &quad, &gp);

                // ---- (b) EQUILIBRIUM distribution ----
                gp.distr_pars = NuEquilibriumParams(&gp.eos_pars);
                ComputeM1DensitiesEq(&gp.eos_pars, &gp.distr_pars, &gp.m1_pars);
                BS_REAL neq0[4], neq1[4];
                for (int s = 0; s < 4; ++s)
                { neq0[s] = gp.m1_pars.n[s]; neq1[s] = gp.m1_pars.J[s]; }
                for (int s = 0; s < total_num_species; ++s) gp.m1_pars.chi[s] = 1. / 3.;
                M1Opacities eq = ComputeM1Opacities(&quad, &quad, &gp);

                for (int s = 0; s < 4; ++s)
                {
                    d_eta0_re(i, s) = re.eta_0[s];   d_eta_re(i, s) = re.eta[s];
                    d_kap0_re(i, s) = re.kappa_0_a[s]; d_kap_re(i, s) = re.kappa_a[s];
                    d_eta0_eq(i, s) = eq.eta_0[s];   d_eta_eq(i, s) = eq.eta[s];
                    d_kap0_eq(i, s) = eq.kappa_0_a[s]; d_kap_eq(i, s) = eq.kappa_a[s];
                    d_kaps_re(i, s) = re.kappa_s[s]; d_kaps_eq(i, s) = eq.kappa_s[s];
                    // proposed fix: scale eq emissivity by occupation ratio
                    // (capped at 1). In NEPS-only mode the whole emissivity IS
                    // the non-thermal part, so this is the exact fix prototype.
                    BS_REAL f0 = (neq0[s] > 0.0)
                                     ? Kokkos::fmin(nfield0[s] / neq0[s], 1.0) : 1.0;
                    BS_REAL f1 = (neq1[s] > 0.0)
                                     ? Kokkos::fmin(nfield1[s] / neq1[s], 1.0) : 1.0;
                    d_eta0_eqc(i, s) = eq.eta_0[s] * f0;
                    d_eta_eqc(i, s)  = eq.eta[s]   * f1;
                }
            });
        Kokkos::fence();

        // -------------------------------------------------------------------
        // Copy back and print.
        // -------------------------------------------------------------------
        auto eta0_re = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_eta0_re);
        auto eta_re  = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_eta_re);
        auto kap0_re = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_kap0_re);
        auto kap_re  = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_kap_re);
        auto eta0_eq = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_eta0_eq);
        auto eta_eq  = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_eta_eq);
        auto kap0_eq = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_kap0_eq);
        auto kap_eq  = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_kap_eq);
        auto kaps_re = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_kaps_re);
        auto kaps_eq = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_kaps_eq);
        auto eta0_eqc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_eta0_eqc);
        auto eta_eqc  = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_eta_eqc);

        // Header. Species order: nue, anue, nux, anux.
        printf("# mode=%s\n",
               neps_only ? "NEPS-only" : (iso_only ? "iso-only" : "all-reactions"));
        printf("# r[km] rho[g/cc] T[MeV]");
        const char* tag[12] = {"eta0_re", "eta_re", "kap0_re", "kap_re",
                               "eta0_eq", "eta_eq", "kap0_eq", "kap_eq",
                               "kaps_re", "kaps_eq", "eta0_eqc", "eta_eqc"};
        const char* sp[4] = {"nue", "anue", "nux", "anux"};
        for (int t = 0; t < 12; ++t)
            for (int s = 0; s < 4; ++s) printf(" %s_%s", tag[t], sp[s]);
        printf("\n");
        printf("# eta0,eta scaled x%.0e (cm^-3 s^-1); kappa scaled x%.0e (cm^-1)\n",
               C_ETA, C_KAP);

        for (int i = 0; i < nz; ++i)
        {
            printf("%.8e %.8e %.8e", h_r(i), h_rho(i), h_T(i));
            for (int s = 0; s < 4; ++s) printf(" %.8e", eta0_re(i, s) * C_ETA);
            for (int s = 0; s < 4; ++s) printf(" %.8e", eta_re(i, s) * C_ETA);
            for (int s = 0; s < 4; ++s) printf(" %.8e", kap0_re(i, s) * C_KAP);
            for (int s = 0; s < 4; ++s) printf(" %.8e", kap_re(i, s) * C_KAP);
            for (int s = 0; s < 4; ++s) printf(" %.8e", eta0_eq(i, s) * C_ETA);
            for (int s = 0; s < 4; ++s) printf(" %.8e", eta_eq(i, s) * C_ETA);
            for (int s = 0; s < 4; ++s) printf(" %.8e", kap0_eq(i, s) * C_KAP);
            for (int s = 0; s < 4; ++s) printf(" %.8e", kap_eq(i, s) * C_KAP);
            for (int s = 0; s < 4; ++s) printf(" %.8e", kaps_re(i, s) * C_KAP);
            for (int s = 0; s < 4; ++s) printf(" %.8e", kaps_eq(i, s) * C_KAP);
            for (int s = 0; s < 4; ++s) printf(" %.8e", eta0_eqc(i, s) * C_ETA);
            for (int s = 0; s < 4; ++s) printf(" %.8e", eta_eqc(i, s) * C_ETA);
            printf("\n");
        }
    }
    Kokkos::finalize();
    return 0;
}
