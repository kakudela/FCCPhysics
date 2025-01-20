
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 1
processList = {
    'kkmcee_ee_mumu_ecm240_ISRenhanced':     {'fraction':fraction},
    'kkmcee_ee_mumu_ecm240_ISRsuppressed':     {'fraction':fraction},
    'kkmcee_ee_mumu_ecm240':       {'fraction':fraction},
}





inputDir = "test/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/ewk/radiative_return/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 64

# scale the histograms with the cross-section and integrated luminosity
doScale = False
intLumi = 1

# define histograms
bins_p_mu = (4000, 0, 400) # 100 MeV bins
bins_m_ll = (4000, 0, 400) # 100 MeV bins
bins_p_ll = (200, 0, 200) # 1 GeV bins
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (20000, 120, 140) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, 0, 5)
bins_phi = (500, -5, 5)
bins_aco = (400, -4, 4)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 5)
bins_dR = (1000, 0, 10)

bins_cat = (10, 0, 10)
bins_resolution = (10000, 0.95, 1.05)

bins_cos = (100, -1, 1)
bins_aco = (1000,-360,360)

def build_graph(df, dataset):

    results = []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    

    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Lepton0", "Muon#0.index")

    results.append(df.Histo1D(("evts_initial", "", *bins_count), "weight"))


    # reco muons
    df = df.Define("leps_all", "FCCAnalyses::ReconstructedParticle::get(Lepton0, ReconstructedParticles)")
    df = df.Define("leps_all_p", "FCCAnalyses::ReconstructedParticle::get_p(leps_all)")
    df = df.Define("leps_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leps_all)")
    df = df.Define("leps_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leps_all)")
    df = df.Define("leps_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(leps_all)")
    df = df.Define("leps_all_no", "FCCAnalyses::ReconstructedParticle::get_n(leps_all)")


    # Cut 0: number of muons
    results.append(df.Histo1D(("evts_initial", "", *bins_count), "weight"))
    results.append(df.Histo1D(("leps_all_no_none", "", *bins_count), "leps_all_no"))
    df = df.Filter("leps_all_no == 2")
    results.append(df.Histo1D(("evts_cut0", "", *bins_count), "weight"))

    # construct Lorentz vectors of the leptons
    df = df.Define("leps_tlv", "FCCAnalyses::makeLorentzVectors(leps_all)")
    df = df.Define("m_inv", "(leps_tlv[0] + leps_tlv[1]).M()")
    df = df.Define("missingEnergy_vec", "FCCAnalyses::missingEnergy(91., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_vec[0].energy")
    df = df.Define("visibleEnergy", "FCCAnalyses::visibleEnergy(ReconstructedParticles)") 

    # Cut 1: cosTheta
    df = df.Define("CosTheta", "FCCAnalyses::Vec_f ret; for(auto & theta: leps_all_theta) ret.push_back(cos(theta)); return ret;")
    df = df.Define("negative_muon", "(leps_all_q[0] < leps_all_q[1]) ? CosTheta[0] : CosTheta[1]")

    #df = df.Filter("abs(CosTheta[0]) <= .9 && abs(CosTheta[1]) <= .9") #LEP
    results.append(df.Histo1D(("CosTheta_none", "", *bins_cos), "negative_muon"))
    df = df.Filter("abs(CosTheta[0]) <= .98 && abs(CosTheta[1]) <= .98") #FCC
    results.append(df.Histo1D(("evts_cut1", "", *bins_count), "weight"))


    # Cut 2: Maximum momentum
    df = df.Define("maximum_momentum","auto max_iter = std::max_element(leps_all_p.begin(), leps_all_p.end()); return *max_iter;")
    results.append(df.Histo1D(("maximum_momentum_none", "", *bins_p_mu), "maximum_momentum"))

    #df =df.Filter("maximum_momentum > 27.36") #LEP
    df = df.Filter("maximum_momentum > 36") #FCC
    df = df.Define("maxp_overE","maximum_momentum / 45.6") # 45.6 beam energy
    results.append(df.Histo1D(("evts_cut2", "", *bins_count), "weight"))


    # Cut 3: Acolinearity
    df = df.Define("acolinearity_rad", "FCCAnalyses::acolinearity(leps_all)")
    df = df.Define("acolinearity_d", "acolinearity_rad * 57.2958")
    results.append(df.Histo1D(("acolinearity_d_none", "", *bins_aco), "acolinearity_d"))

    df = df.Filter("acolinearity_d < 90")
    results.append(df.Histo1D(("evts_cut3", "", *bins_count), "weight"))

    results.append(df.Histo1D(("leps_all_p", "", *bins_p_mu), "leps_all_p"))
    results.append(df.Histo1D(("leps_all_theta", "", *bins_theta), "leps_all_theta"))
    results.append(df.Histo1D(("leps_all_phi", "", *bins_phi), "leps_all_phi"))
    results.append(df.Histo1D(("leps_all_q", "", *bins_charge), "leps_all_q"))
    results.append(df.Histo1D(("leps_all_no", "", *bins_count), "leps_all_no"))

    results.append(df.Histo1D(("m_inv", "", *bins_m_ll), "m_inv"))
    results.append(df.Histo1D(("missingEnergy", "", *bins_m_ll), "missingEnergy"))
    results.append(df.Histo1D(("visibleEnergy", "", *bins_m_ll), "visibleEnergy"))

    df = df.Define("theta_plus", "(leps_all_q[0] > 0) ? leps_all_theta[0] : leps_all_theta[1]")
    df = df.Define("theta_minus", "(leps_all_q[0] < 0) ? leps_all_theta[0] : leps_all_theta[1]")
    df = df.Define("cos_theta_plus", "cos(theta_plus)")
    df = df.Define("cos_theta_minus", "cos(theta_minus)")
    df = df.Define("cosThetac", "(sin(theta_plus-theta_minus))/(sin(theta_plus)+sin(theta_minus))")

    results.append(df.Histo1D(("theta_plus", "", *bins_theta), "theta_plus"))
    results.append(df.Histo1D(("theta_minus", "", *bins_theta), "theta_minus"))
    results.append(df.Histo1D(("cos_theta_plus", "", *bins_cos), "cos_theta_plus"))
    results.append(df.Histo1D(("cos_theta_minus", "", *bins_cos), "cos_theta_minus"))

    results.append(df.Histo1D(("cosThetac", "", *bins_cos), "cosThetac"))
    results.append(df.Histo1D(("evts_final", "", *bins_count), "weight"))
    results.append(df.Histo1D(("CosTheta", "", *bins_cos), "negative_muon"))
    results.append(df.Histo1D(("acolinearity_d", "", *bins_aco), "acolinearity_d"))
    results.append(df.Histo1D(("maximum_momentum", "", *bins_p_mu), "maximum_momentum"))
    results.append(df.Histo1D(("maxp_overE", "", *bins_resolution), "maxp_overE"))


    return results, weightsum