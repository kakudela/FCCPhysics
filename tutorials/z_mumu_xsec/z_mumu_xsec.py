
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all guns
processList = {
    'wzp6_ee_mumu_ecm91p2': {'fraction':1},
    'wzp6_ee_tautau_ecm91p2': {'fraction':1},
    'p8_ee_gaga_mumu_ecm91p2': {'fraction':1},
}



inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/tutorials/z_mumu_xsec/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 128

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 44.84 # 44.84 pb-1 = LEP, 100e6=100 ab-1 = FCCee

# define histograms
bins_p_mu = (20000, 0, 200) # 10 MeV bins
bins_m_ll = (20000, 0, 200) # 10 MeV bins
bins_p_ll = (20000, 0, 200) # 10 MeV bins

bins_theta = (500, -5, 5)
bins_phi = (500, -5, 5)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_cos = (100, -1, 1)
bins_norm = (200, 0, 2)
bins_nparticles = (200, 0, 200)
bins_aco = (800,-4,4)

def build_graph(df, dataset):

    hists = []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muons", "Muon#0.index")


    hists.append(df.Histo1D(("evts_initial", "", *bins_count), "weight"))

    # reco muons
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muons, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_costheta", "FCCAnalyses::get_costheta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    # define cos(theta) of the muons
    hists.append(df.Histo1D(("muons_all_costheta", "", *bins_cos), "muons_all_costheta"))
    df = df.Define("muons", "FCCAnalyses::sel_range(0, 0.97, true)(muons_all, muons_all_costheta)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")

    # select 2 opposite charge muons
    df = df.Filter("muons_no == 2 && (muons_q[0] + muons_q[1]) == 0")

    # max normalized muon momentum > 0.6
    df = df.Define("muon_max_p", "(muons_p[0] > muons_p[1]) ? muons_p[0] : muons_p[1]")
    df = df.Define("muon_max_p_norm", "muon_max_p/45.6")
    hists.append(df.Histo1D(("muon_max_p_norm", "", *bins_norm), "muon_max_p_norm"))
    df = df.Filter("muon_max_p_norm > 0.6")

    # acolinearity between 2 muons
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(muons)")
    hists.append(df.Histo1D(("acolinearity", "", *bins_aco), "acolinearity"))

    # plot invariant mass of both muons
    df = df.Define("leps_tlv", "FCCAnalyses::makeLorentzVectors(muons)")
    df = df.Define("invariant_mass", "(leps_tlv[0]+leps_tlv[1]).M()")
    hists.append(df.Histo1D(("invariant_mass", "", *bins_m_ll), "invariant_mass"))


    hists.append(df.Histo1D(("evts_final", "", *bins_count), "weight"))

    return hists, weightsum
