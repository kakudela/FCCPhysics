
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 0.05
processList = {
    'wzp6_ee_nunuH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_nunuH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_nunuH_Hss_ecm240':          {'fraction':fraction},
}





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/clustering_tagging/histmaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 24

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000

# define histograms
bins_m = (250, 0, 250) # 100 MeV bins
bins_score = (100, 0, 1)


################################################################################
## load modules and files for jet clustering and flavor tagging
################################################################################
## latest particle transformer model, trainied on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc_v1"

# model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

# model files locally stored on /eos
model_dir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

# get local file, else download from url
def get_file_path(url, filename):
    import os
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        import urllib.request
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)

weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import ExclusiveJetClusteringHelper
from examples.FCCee.weaver.config import collections, njets

################################################################################


def build_graph(df, dataset):

    hists, cols = [], []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")


    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")


    # run jet clustering
    njets = 2
    jetClusteringHelper = ExclusiveJetClusteringHelper("ReconstructedParticles", njets)
    df = jetClusteringHelper.define(df)


    # run flavor tagging
    jetFlavourHelper = JetFlavourHelper(collections, jetClusteringHelper.jets, jetClusteringHelper.constituents, "")
    df = jetFlavourHelper.define(df) # define variables
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df) # run inference

    # jetFlavourHelper adds new columns for each jet flavor (recojet_isB/C/S/...)
    # each column is a vector of probabilities per jet to be that specific flavor
    hists.append(df.Histo1D(("recojet_isB", "", *bins_score), "recojet_isB"))
    hists.append(df.Histo1D(("recojet_isC", "", *bins_score), "recojet_isC"))
    hists.append(df.Histo1D(("recojet_isS", "", *bins_score), "recojet_isS"))

    # compute the invariant mass of the jets
    df = df.Define("jet_p4", f"JetConstituentsUtils::compute_tlv_jets({jetClusteringHelper.jets})")
    df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
    hists.append(df.Histo1D(("event_invariant_mass", "", *bins_m), "event_invariant_mass"))
    return hists, weightsum


