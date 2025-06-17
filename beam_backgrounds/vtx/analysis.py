
import sys, os, glob, math
import ROOT
import logging
import time

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

ROOT.EnableImplicitMT() # use all cores
#ROOT.DisableImplicitMT() # single core

# load libraries
ROOT.gInterpreter.Declare('#include "functions.h"')


layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

bins_theta = (36, 0, 180)
bins_phi = (72, -180, 180)
bins_z = (int(max_z/2), -max_z, max_z)
bins_layer = (10, 0, 10)

def analysis(input_files, output_file):

    hists = []
    df = ROOT.RDataFrame("events", input_files)


    df = df.Define("hits_x", "getSimHitPosition_x(VertexBarrelCollection)")
    df = df.Define("hits_y", "getSimHitPosition_y(VertexBarrelCollection)")
    df = df.Define("hits_z", "getSimHitPosition_z(VertexBarrelCollection)")
    df = df.Define("hits_r", "getSimHitPosition_r(VertexBarrelCollection)")
    df = df.Define("hits_theta", "getSimHitPosition_theta(VertexBarrelCollection, true)")
    df = df.Define("hits_phi", "getSimHitPosition_phi(VertexBarrelCollection, true)")
    df = df.Define("hits_layer", "getSimHitLayer(hits_r, {14, 36, 58})")

    df = df.Define("hits_layer0", "hits_layer == 0")
    df = df.Define("hits_z_layer0", "hits_z[hits_layer0]")
    df = df.Define("hits_r_layer0", "hits_r[hits_layer0]")
    df = df.Define("hits_theta_layer0", "hits_theta[hits_layer0]")
    df = df.Define("hits_phi_layer0", "hits_phi[hits_layer0]")

    hists.append(df.Histo1D(("z", "", *bins_z), "hits_z"))
    hists.append(df.Histo1D(("theta", "", *bins_theta), "hits_theta"))
    hists.append(df.Histo1D(("phi", "", *bins_phi), "hits_phi"))

    hists.append(df.Histo1D(("z_layer0", "", *bins_z), "hits_z_layer0"))
    hists.append(df.Histo1D(("theta_layer0", "", *bins_theta), "hits_theta_layer0"))
    hists.append(df.Histo1D(("phi_layer0", "", *bins_phi), "hits_phi_layer0"))

    hists.append(df.Histo2D(("z_phi", "", *(bins_z + bins_phi)), "hits_z", "hits_phi"))
    hists.append(df.Histo2D(("theta_phi", "", *(bins_theta + bins_phi)), "hits_theta", "hits_phi"))

    hists.append(df.Histo2D(("z_phi_layer0", "", *(bins_z + bins_phi)), "hits_z_layer0", "hits_phi_layer0"))
    hists.append(df.Histo2D(("theta_phi_layer0", "", *(bins_theta + bins_phi)), "hits_theta_layer0", "hits_phi_layer0"))


    hists.append(df.Histo1D(("layer", "", *bins_layer), "hits_layer"))

    fout = ROOT.TFile(output_file, "RECREATE")
    for h in hists:
        h.Write()
    fout.Close()



if __name__ == "__main__":

    input_dir, output_file = "/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z256", "CLD_o2_v05_FCCee_Z_4IP_04may23_FCCee_Z256.root"

    logger.info(f"Start analysis")
    input_files = glob.glob(f"{input_dir}/*.root")
    logger.info(f"Found {len(input_files)} input files")

    logger.info(f"Start analysis")
    time0 = time.time()
    analysis(input_files, output_file)
    time1 = time.time()

    logger.info(f"Analysis done in {int(time1-time0)} seconds")
    logger.info(f"Output saved to {output_file}")

