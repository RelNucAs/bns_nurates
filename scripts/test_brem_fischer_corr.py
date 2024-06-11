import numpy as np
import bns_nurates as bns
import matplotlib.pyplot as plt

c   = 2.997924562E10  # speed of light [cm s-1]
MeV = 1.602176634E-6  # MeV to erg conversion
mb = 922.293227699518 # baryonic mass to convert from trho to nb (taken from DD2 stellar collapse table)

##############################################################
# Reproduce behaviour of NN bremsstrahlung kernel correction #
# due to medium effects as in Fig. 2 in Fischer2016          #
##############################################################

# 1D profile taken from BNS simulation + DD2 stellar collapse table
# DD2_M12980-12980_M1_LR (rl = 3, x axis, t = 23.75 ms post-merger)
rho_array = [692091315217157.5, 669166828429746.1,
                   580438042821992.8, 437877988168366.0,
                   291698155663313.2, 198158314550266.25,
                   129660029205976.39, 78418507867115.28,
                   41111324737689.45, 20104589964605.434,
                   10112305551424.934, 5120554937516.626,
                   3931905610863.119, 3101221362480.4575,
                   2841883715228.641, 1904652507888.0417,
                   1447581967724.244, 1097844665855.3694,
                   809848197032.8032, 678400239111.4371,
                   614996689981.9164, 563975144576.5432,
                   419488592639.89264, 238174687523.39374,
                   216916023317.7148, 203505898623.57486,
                   190632803840.05273, 176332848687.99338,
                   165532917912.70142, 155375850190.61267,
                   144118742708.48807, 133957242598.19185,
                   125037750325.66563, 98439655745.14645,
                   68791687501.1082, 62684501661.18724,
                   59893786664.206474, 57361958246.56711,
                   54823911319.927765, 51824720479.58496,
                   48077899849.72525, 44332786347.6009, 40396031300.29652,
                   33649486626.958153, 24228844845.778717,
                   18954553023.575504, 17208094425.638344,
                   16109014650.672195, 14903239392.246632,
                   13698855289.327766, 12561761202.727837,
                   11484927847.485151, 10563390454.121267,
                   9853829832.527843, 8825305120.243359,
                   6772165303.441048, 4895666579.401289,
                   4152282118.143988, 3923217612.1331115,
                   3772196875.1917953, 3618882322.4968257,
                   3422427412.1773977, 3226191652.7353063,
                   3082488470.558036, 2895081008.098228,
                   2670198745.6138015, 2423878919.397487,
                   2167962808.1431403, 1951153197.2504299,
                   1810100945.2727358, 1734820836.078078,
                   1690880262.327674]
temp_array = [12.18228530883789, 14.105046272277832,
                   23.41455841064453, 34.46710968017578,
                   33.40366744995117, 23.287614822387695,
                   16.319942474365234, 11.190557479858398,
                   7.722237586975098, 5.983924388885498,
                   5.908384799957275, 7.3729248046875, 7.308124542236328,
                   7.221435070037842, 6.408595085144043, 6.67584753036499,
                   6.7033796310424805, 7.227153778076172,
                   7.3906331062316895, 7.476637363433838,
                   7.607393741607666, 7.712704658508301,
                   7.450197219848633, 5.528027057647705,
                   5.2435102462768555, 5.084035873413086,
                   4.917929172515869, 4.7169976234436035,
                   4.555782318115234, 4.422448635101318,
                   4.293033123016357, 4.183816432952881,
                   4.167706489562988, 3.8147246837615967,
                   3.108804225921631, 2.9898409843444824,
                   2.9495182037353516, 2.909569025039673,
                   2.868633508682251, 2.830944776535034,
                   2.7984256744384766, 2.781472682952881,
                   2.776703357696533, 2.670297145843506,
                   2.3788483142852783, 2.2135019302368164,
                   2.1888999938964844, 2.198323965072632,
                   2.202788829803467, 2.2058680057525635,
                   2.2142844200134277, 2.2269446849823,
                   2.2468981742858887, 2.2837469577789307,
                   2.301438808441162, 2.202881097793579,
                   2.0407941341400146, 1.979393482208252,
                   1.9789927005767822, 1.9739326238632202,
                   1.9582414627075195, 1.932193636894226,
                   1.9075164794921875, 1.903809905052185,
                   1.886286973953247, 1.8593289852142334,
                   1.821901559829712, 1.7723792791366577,
                   1.7273764610290527, 1.697688341140747,
                   1.6849199533462524, 1.6812838315963745]
ye_array = [0.07206907868385315, 0.07042089849710464,
                   0.06655445694923401, 0.06538038700819016,
                   0.06490999460220337, 0.0604977086186409,
                   0.05783260986208916, 0.056365665048360825,
                   0.05538692697882652, 0.05387461557984352,
                   0.0505199208855629, 0.06126539781689644,
                   0.06218162178993225, 0.06856312602758408,
                   0.056119561195373535, 0.0719020888209343,
                   0.08237294107675552, 0.10490885376930237,
                   0.12769991159439087, 0.14383022487163544,
                   0.15403234958648682, 0.16205169260501862,
                   0.1636190563440323, 0.16106148064136505,
                   0.16213330626487732, 0.16218270361423492,
                   0.16071555018424988, 0.15745918452739716,
                   0.15315869450569153, 0.1485365927219391,
                   0.14341485500335693, 0.13799211382865906,
                   0.1360314041376114, 0.13574004173278809,
                   0.1332605928182602, 0.13094979524612427,
                   0.12898363173007965, 0.127104252576828,
                   0.12559671700000763, 0.12517079710960388,
                   0.1264037787914276, 0.12937739491462708,
                   0.13386492431163788, 0.13922850787639618,
                   0.1444927453994751, 0.14905206859111786,
                   0.15370570123195648, 0.15895096957683563,
                   0.16489309072494507, 0.1718137562274933,
                   0.17991381883621216, 0.18859565258026123,
                   0.1969536542892456, 0.2061663120985031,
                   0.2201358526945114, 0.23926915228366852,
                   0.25589221715927124, 0.2647881507873535,
                   0.2680870592594147, 0.26852211356163025,
                   0.2685197591781616, 0.2682563066482544,
                   0.2677820026874542, 0.2674899399280548,
                   0.2675788700580597, 0.2679418623447418,
                   0.26817047595977783, 0.2685425281524658,
                   0.2686082422733307, 0.26839321851730347,
                   0.26807400584220886, 0.26779672503471375]
xn_array = [0.92793092, 0.9295791 , 0.93344554, 0.93461961, 0.93507999,
       0.93329893, 0.9258373 , 0.91614851, 0.8799617 , 0.85189414,
       0.87530481, 0.88251547, 0.88693936, 0.88155007, 0.89744375,
       0.88493792, 0.87773189, 0.86017774, 0.8416432 , 0.82784598,
       0.81946109, 0.81310667, 0.81607829, 0.81823622, 0.81678719,
       0.81655958, 0.81802005, 0.82143359, 0.82581046, 0.83083975,
       0.83668495, 0.84301591, 0.84634596, 0.84713078, 0.84541495,
       0.84692124, 0.84918174, 0.85139899, 0.85309503, 0.8540109 ,
       0.85371738, 0.85193011, 0.84892642, 0.84325418, 0.83299337,
       0.82344878, 0.81935266, 0.81664634, 0.81300552, 0.8082741 ,
       0.80251167, 0.7961029 , 0.78988785, 0.7827797 , 0.77047785,
       0.75188742, 0.73398594, 0.72499355, 0.72256564, 0.72249919,
       0.72248079, 0.72262627, 0.72321423, 0.72415774, 0.72441023,
       0.72352645, 0.72287652, 0.72132927, 0.71943496, 0.71926247,
       0.71934847, 0.71987908]
xp_array = [0.07206908, 0.0704209 , 0.06655446, 0.06538039, 0.06490358,
       0.05699355, 0.0492716 , 0.04261913, 0.02509262, 0.00921773,
       0.00858385, 0.02186617, 0.0251778 , 0.03099271, 0.0223421 ,
       0.0379025 , 0.04960642, 0.0745126 , 0.09997873, 0.11773718,
       0.12934814, 0.13872521, 0.14435458, 0.14156573, 0.1423075 ,
       0.14221945, 0.14077358, 0.13770176, 0.13352772, 0.12932437,
       0.12490025, 0.12033988, 0.11960423, 0.11978175, 0.11358023,
       0.11050939, 0.10880744, 0.1072342 , 0.10589265, 0.10590035,
       0.10796259, 0.11198439, 0.11778704, 0.12275585, 0.12307465,
       0.12260782, 0.12773444, 0.13541417, 0.14355571, 0.15256793,
       0.1629082 , 0.17377063, 0.18418843, 0.19542956, 0.21099654,
       0.23062263, 0.24594529, 0.25472624, 0.25888026, 0.25967638,
       0.25964992, 0.25926584, 0.25890082, 0.25925109, 0.25967453,
       0.25951326, 0.25931785, 0.25851295, 0.25674873, 0.25614342,
       0.25558894, 0.2555624]
mu_e_array = [187.72394055, 183.25068812, 164.36552894, 133.54417358,
       111.44468074, 103.34950608,  92.53320281,  79.91153108,
        65.06574737,  50.82075988,  38.5432493 ,  30.06248582,
        26.93478067,  25.37266734,  23.23941738,  21.225626  ,
        19.74669929,  18.52056055,  17.14849435,  16.42822125,
        15.89047188,  15.35206668,  13.15715074,  11.9389234 ,
        11.82740271,  11.66446343,  11.45750342,  11.15798219,
        10.86304984,  10.5127558 ,  10.08389052,   9.6296404 ,
         9.19359941,   8.54292902,   7.97723758,   7.70473889,
         7.51469687,   7.34608546,   7.18296008,   7.00604524,
         6.81120159,   6.61636295,   6.39989644,   6.06239646,
         5.58501839,   5.20334945,   5.04631487,   4.9396538 ,
         4.80628361,   4.67864782,   4.54797656,   4.40328638,
         4.25820998,   4.12203718,   3.96892945,   3.6719754 ,
         3.33798146,   3.1521631 ,   3.05478831,   2.98101633,
         2.92133495,   2.84833156,   2.76574012,   2.68009031,
         2.5850812 ,   2.48222588,   2.36691484,   2.25441317,
         2.15159874,   2.07553392,   2.02689523,   1.98964012]
mu_n_array = [1229.693481  , 1210.60682051, 1139.94826083, 1045.74465221,
        984.35256375,  966.23683227,  957.40774254,  952.47153259,
        950.00857696,  948.27909614,  945.18346311,  937.73846784,
        936.150217  ,  934.75124676,  936.91081226,  933.38430976,
        931.43720512,  927.24913257,  924.1646306 ,  922.34115962,
        920.87159003,  919.6223369 ,  918.80919981,  925.64874707,
        926.70473739,  927.25119893,  927.85813259,  928.61510958,
        929.24344049,  929.7250842 ,  930.15532339,  930.48702285,
        930.30984072,  931.38637525,  934.22933539,  934.64971154,
        934.75715604,  934.87808394,  934.99356425,  935.05607284,
        935.04251065,  934.9101291 ,  934.67095899,  934.82125737,
        935.82973706,  936.32734794,  936.26545467,  936.05611935,
        935.83935138,  935.62127397,  935.35480424,  935.04544036,
        934.69958986,  934.24622114,  933.82514109,  933.92572677,
        934.45755781,  934.5974816 ,  934.48108759,  934.44270169,
        934.49083764,  934.59654212,  934.68756101,  934.63738359,
        934.66220704,  934.73981418,  934.87447451,  935.09393873,
        935.29423056,  935.42458746,  935.46463256,  935.45238075]
mu_p_array = [1019.47213129, 1002.44484433,  938.18639242,  851.90147807,
        813.14719223,  828.35915045,  850.34385934,  875.4058322 ,
        896.69708561,  907.09403174,  909.82574673,  905.85191462,
        906.25567923,  907.28178984,  910.0142341 ,  909.81919836,
        909.96484521,  907.627476  ,  906.68273896,  906.11245916,
        905.22370151,  904.41441438,  904.40932244,  914.54154354,
        916.14054606,  916.97303526,  917.81845114,  918.80812956,
        919.56548344,  920.1271072 ,  920.62268645,  920.97896876,
        920.79474578,  922.58067984,  926.6638408 ,  927.24566849,
        927.3846841 ,  927.53949122,  927.69410737,  927.82937688,
        927.9376956 ,  927.94999281,  927.87164175,  928.37640793,
        929.99345306,  930.84008753,  930.92244586,  930.83390108,
        930.74584155,  930.66895481,  930.54881042,  930.37852552,
        930.14777392,  929.79341388,  929.55893283,  930.03631448,
        930.94240124,  931.24392835,  931.16584006,  931.13831707,
        931.20163154,  931.32977529,  931.4409528 ,  931.39430493,
        931.43915085,  931.54788834,  931.72330051,  931.99360126,
        932.23425715,  932.3920608 ,  932.44104406,  932.43073554]

# Neutrino energy
e_nu = 10. # MeV

# Populate quad structure
quad = bns.MyQuadrature()
quad.nx = 50
quad.dim = 1
quad.x1 = 0
quad.x2 = 1
bns.GaussLegendreMultiD(quad)

# bns_nurates structures (to be properly populated at each grid point)
eos_pars = bns.MyEOSParams()
m1_pars = bns.M1Quantities()
grey_pars = bns.GreyOpacityParams()

# Reactions (turn on only NN bremmstrahlung reaction)
opacity_flags = {'use_abs_em': 0, 'use_pair': 0, 'use_brem': 1, 'use_inelastic_scatt': 0, 'use_iso': 0}
grey_pars.opacity_flags = opacity_flags

# Corrections (switch on 'use_NN_medium_corr' when evaluating corrected rates)
opacity_pars = {'use_dU': False, 'use_dm_eff': False, 'use_WM_ab': False, 'use_WM_sc': False, 'use_decay': False, 'neglect_blocking': False, 'use_BRT_brem': False, 'use_NN_medium_corr': False}
grey_pars.opacity_pars = opacity_pars

# Output (spectral emissivity as a function of grid point and nu flavour)
j_brem_std  = np.zeros((len(rho_array),bns.total_num_species))
j_brem_corr = np.zeros((len(rho_array),bns.total_num_species))

idx = 0
for rho, temp, ye, xn, xp, mu_e, mu_n, mu_p in zip(rho_array, temp_array, ye_array, xn_array, xp_array, mu_e_array, mu_n_array, mu_p_array):
    # popualate EOS structure
    nb = rho * c * c / (mb * MeV)
    eos_pars.nb = nb
    eos_pars.temp = temp
    eos_pars.ye = ye
    eos_pars.yn = xn
    eos_pars.yp = xp 
    eos_pars.mu_n = mu_n
    eos_pars.mu_p = mu_p
    eos_pars.mu_e = mu_n
    grey_pars.eos_pars = eos_pars

    # reconstruct equilibrium distribution function for neutrinos
    distr_pars = bns.NuEquilibriumParams(eos_pars)
    grey_pars.distr_pars = distr_pars

    nu_dens = bns.ComputeM1DensitiesEq(eos_pars, distr_pars)
    nvals = [nu_dens.integrand[0], nu_dens.integrand[2], nu_dens.integrand[4]]
    Jvals = [nu_dens.integrand[1], nu_dens.integrand[3], nu_dens.integrand[5]]
    chivals = [1./3.] * bns.total_num_species
    if (bns.total_num_species == 4):
        nvals = nvals + [nvals[-1]]
        Jvals = Jvals + [Jvals[-1]]
    m1_pars.n = nvals
    m1_pars.J = Jvals
    m1_pars.chi = chivals
    grey_pars.m1_pars = m1_pars

    # evaluate non corrected rates
    opacity_pars['use_NN_medium_corr'] = False
    grey_pars.opacity_pars = opacity_pars
    tmp = bns.ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, quad, grey_pars)
    j_brem_std[idx] = tmp["j"]

    # evaluate corrected rates
    opacity_pars["use_NN_medium_corr"] = True
    grey_pars.opacity_pars = opacity_pars
    tmp = bns.ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, quad, grey_pars)
    j_brem_corr[idx] = tmp["j"]

    idx = idx + 1

# compute NN*/NN ratio (same for all flavours)
ratio = (j_brem_corr / j_brem_std)[:,0]

print("Uncorrected spectral NN bremsstrahlung emissivity for each neutrino flavour:")
print(j_brem_std)
print("")

print("Corrected spectral NN bremsstrahlung emissivity for each neutrino flavour:")
print(j_brem_corr)
print("")

print("Ratio:")
print(ratio)
print("")

# Reproduce Fig.2 in Fischer2016
fig = plt.figure()
plt.title("Fig. 2 in Fischer (2016)")
plt.plot(rho_array, ratio, color="tab:orange", ls="-", label="NN*/NN")
plt.xlabel(r"$\rho~[{\rm g}~{\rm cm}^{-3}]$")
plt.xscale('log')
plt.xlim((1.0E+10,1.0E+15))
plt.ylim((0.0,1.4))
plt.axhline(y = 1.     , ls="--", color="tab:grey", alpha=0.5)
plt.axvline(x = 2.5E+14, ls="--", color="tab:grey", alpha=0.5)
plt.legend()
plt.show()