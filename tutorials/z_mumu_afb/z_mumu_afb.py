
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all guns
processList = {
    'kkmcee_ee_mumu_ecm91p2': {'fraction':1},
    'wzp6_ee_mumu_ecm91p2': {'fraction':1},
    'p8_ee_Zmumu_ecm91': {'fraction':1},
}



inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/tutorials/z_mumu_afb/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 128

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 1 # normalize to 1 pb-1

# define histograms
bins_p_mu = (20000, 0, 200) # 10 MeV bins
bins_m_ll = (10000, 0, 100) # 10 MeV bins
bins_p_ll = (20000, 0, 200) # 10 MeV bins

bins_theta = (500, -5, 5)
bins_phi = (500, -5, 5)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_cos = (100, -1, 1)


def build_graph(df, dataset):

    hists = []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muons", "Muon#0.index")


    # reco muons
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muons, ReconstructedParticles)")

    # select muons with at least 10 GeV momentum
    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(10, 1000)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    # construct Lorentz vectors of the leptons
    df = df.Define("leps_tlv", "FCCAnalyses::makeLorentzVectors(muons)")

    # apply some basic filters (number of muons, opposite sign, invariant mass)
    df = df.Filter("muons_no == 2")
    df = df.Filter("(muons_q[0] + muons_q[1]) == 0")
    df = df.Define("m_inv", "(leps_tlv[0]+leps_tlv[1]).M()")
    df = df.Filter("m_inv >= 50")

    df = df.Define("visibleEnergy", "FCCAnalyses::visibleEnergy(ReconstructedParticles)")
    hists.append(df.Histo1D(("visibleEnergy", "", *bins_m_ll), "visibleEnergy"))
    hists.append(df.Histo1D(("m_inv", "", *bins_m_ll), "m_inv"))

    hists.append(df.Histo1D(("muons_p", "", *bins_p_mu), "muons_p"))
    hists.append(df.Histo1D(("muons_theta", "", *bins_theta), "muons_theta"))
    hists.append(df.Histo1D(("muons_phi", "", *bins_phi), "muons_phi"))
    hists.append(df.Histo1D(("muons_q", "", *bins_charge), "muons_q"))


    # select the corresponding gen-level muons
    df = df.Define("muons_gen", "FCCAnalyses::getRP2MC(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle)")
    df = df.Define("muons_gen_tlv", "FCCAnalyses::makeLorentzVectors(muons_gen)")

    df = df.Define("gen_muons_p", "FCCAnalyses::MCParticle::get_p(muons_gen)")
    df = df.Define("gen_muons_theta", "FCCAnalyses::MCParticle::get_theta(muons_gen)")
    df = df.Define("gen_muons_phi", "FCCAnalyses::MCParticle::get_phi(muons_gen)")
    df = df.Define("gen_muons_no", "FCCAnalyses::MCParticle::get_n(muons_gen)")
    df = df.Define("gen_muons_q", "FCCAnalyses::MCParticle::get_charge(muons_gen)")

    hists.append(df.Histo1D(("gen_muons_p", "", *bins_p_mu), "gen_muons_p"))
    hists.append(df.Histo1D(("gen_muons_theta", "", *bins_theta), "gen_muons_theta"))
    hists.append(df.Histo1D(("gen_muons_phi", "", *bins_phi), "gen_muons_phi"))
    hists.append(df.Histo1D(("gen_muons_no", "", *bins_count), "gen_muons_no"))


    # calculate the angles (reco)
    df = df.Define("theta_plus", "(muons_q[0] > 0) ? muons_theta[0] : muons_theta[1]")
    df = df.Define("theta_minus", "(muons_q[0] < 0) ? muons_theta[0] : muons_theta[1]")
    df = df.Define("cos_theta_plus", "cos(theta_plus)")
    df = df.Define("cos_theta_minus", "cos(theta_minus)")
    df = df.Define("cosThetac", "(sin(theta_plus-theta_minus))/(sin(theta_plus)+sin(theta_minus))")

    hists.append(df.Histo1D(("theta_plus", "", *bins_theta), "theta_plus"))
    hists.append(df.Histo1D(("theta_minus", "", *bins_theta), "theta_minus"))
    hists.append(df.Histo1D(("cos_theta_plus", "", *bins_cos), "cos_theta_plus"))
    hists.append(df.Histo1D(("cos_theta_minus", "", *bins_cos), "cos_theta_minus"))
    hists.append(df.Histo1D(("cosThetac", "", *bins_cos), "cosThetac"))


    # calculate the angles (gen)
    df = df.Define("gen_theta_plus", "(gen_muons_q[0] > 0) ? gen_muons_theta[0] : gen_muons_theta[1]")
    df = df.Define("gen_theta_minus", "(gen_muons_q[0] < 0) ? gen_muons_theta[0] : gen_muons_theta[1]")
    df = df.Define("gen_cos_theta_plus", "cos(gen_theta_plus)")
    df = df.Define("gen_cos_theta_minus", "cos(gen_theta_minus)")
    df = df.Define("gen_cosThetac", "(sin(gen_theta_plus-gen_theta_minus))/(sin(gen_theta_plus)+sin(gen_theta_minus))")

    hists.append(df.Histo1D(("gen_theta_plus", "", *bins_theta), "gen_theta_plus"))
    hists.append(df.Histo1D(("gen_theta_minus", "", *bins_theta), "gen_theta_minus"))
    hists.append(df.Histo1D(("gen_cos_theta_plus", "", *bins_cos), "gen_cos_theta_plus"))
    hists.append(df.Histo1D(("gen_cos_theta_minus", "", *bins_cos), "gen_cos_theta_minus"))
    hists.append(df.Histo1D(("gen_cosThetac", "", *bins_cos), "gen_cosThetac"))

    return hists, weightsum
