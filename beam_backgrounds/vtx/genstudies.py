from podio import root_io
import glob
import pickle
import argparse
import functions
import math
import ROOT
import os
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument('--plots', help="Plot the energy deposits", action='store_true')
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()


##########################################################################################
# this file is for plotting the number of hits in a 2D map of phi and z and purely as a
# function of phi and theta
##########################################################################################

# folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_o2_v05/guineaPig_andrea_June2024_v23/"
# folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23/"

# dataset_tag, folder = ("256_mapped", "/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05_2T/FCCee_Z_4IP_04may23_FCCee_Z256_edge/")
dataset_tag, folder = ("256_vtx000", "/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05_2T/FCCee_Z_4IP_04may23_FCCee_Z256_vtx000/")

#folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23_vtx000/"
files = glob.glob(f"{folder}/*.root")

outdir_base = "genstudies"
os.makedirs(outdir_base, exist_ok=True)

# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

#me = 0.000511  # Electron mass in GeV



if args.calculate:

    hist_theta_pt_mc = ROOT.TH2D("theta_pt_mc", "", 500, -4, 1, 400, -4, 0)
    hist_theta_pt_hits = ROOT.TH2D("theta_pt_hits", "", 500, -4, 1, 400, -4, 0)
    hist_z_r_mc = ROOT.TH2D("z_r_mc", "", 220, -110, 110, 400, 0, 50)
    hist_z_r_hits = ROOT.TH2D("z_r_hits", "", 220, -110, 110, 400, 0, 50)

    nEvents = 0
    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            nEvents += 1

            for genp in event.get("MCParticles"):
                pdg = genp.getPDG()
                status = genp.getGeneratorStatus()
                if pdg in [11,-11] and status == 1:
                    px = genp.getMomentum().x
                    py = genp.getMomentum().y
                    pz = genp.getMomentum().z
                    x = genp.getVertex().x
                    y = genp.getVertex().y
                    z = genp.getVertex().z
                    r = (x*x + y*y)**0.5
                    energy = genp.getEnergy()
                    lv = ROOT.TLorentzVector(px, py, pz, energy)
                    pt = lv.Pt()
                    theta = lv.Theta()
                    if theta > math.pi/2:
                        theta = math.pi - theta

                    log_theta = math.log10(theta)
                    log_pt = math.log10(pt)

                    hist_theta_pt_mc.Fill(log_theta, log_pt)
                    hist_z_r_mc.Fill(z, r)
                    # print(x,y,z,r)
            # loop over first hits in the first layer
            for hit in event.get("VertexBarrelCollection"):
                hit_r = hit.rho()
                if hit_r > 13.0000 and hit_r < 14.2850: # layer 1
                    isCreatedBySecondary = hit.isProducedBySecondary()
                    if isCreatedBySecondary == False:
                        genp = hit.getMCParticle()
                        pdg = genp.getPDG()
                        status = genp.getGeneratorStatus()
                        if pdg in [11,-11] and status == 1:

                            px = genp.getMomentum().x
                            py = genp.getMomentum().y
                            pz = genp.getMomentum().z
                            x = genp.getVertex().x
                            y = genp.getVertex().y
                            z = genp.getVertex().z
                            r = (x*x + y*y)**0.5
                            energy = genp.getEnergy()
                            lv = ROOT.TLorentzVector(px, py, pz, energy)
                            pt = lv.Pt()
                            theta = lv.Theta()
                            if theta > math.pi/2:
                                theta = math.pi - theta

                            log_theta = math.log10(theta)
                            log_pt = math.log10(pt)
                            hist_theta_pt_hits.Fill(log_theta, log_pt)
                            hist_z_r_hits.Fill(z, r)

        if i > args.maxFiles:
            break

    # normalize the histograms over the number of events, to get the average number of hits per event
    # hist_theta_pt_mc.Scale(1./nEvents)
    # hist_theta_pt_hits.Scale(1./nEvents)

    fOut_name = os.path.join(outdir_base, f"{dataset_tag}.root")
    fOut = ROOT.TFile(fOut_name, "RECREATE")
    hist_theta_pt_mc.Write()
    hist_theta_pt_hits.Write()
    hist_z_r_mc.Write()
    hist_z_r_hits.Write()
    fOut.Close()


if args.plots:

    fIn_name = os.path.join(outdir_base, f"{dataset_tag}.root")
    fIn = ROOT.TFile(fIn_name)

    outdir = outdir_base

    hist_theta_pt_mc = fIn.Get('theta_pt_mc')
    hist_theta_pt_hits = fIn.Get('theta_pt_hits')

    scale1 = 1e4/hist_theta_pt_mc.Integral()
    scale2 = 1e4/hist_theta_pt_hits.Integral()

    functions.plot_2dhist(hist_theta_pt_mc, f"{outdir}/theta_pt_mc_{dataset_tag}.png", f"{dataset_tag}: theta vs pt (MC)", xMin=-4, xMax=0.2, yMin=-4, yMax=0, xLabel=r'$\log_{10}(\theta)\ \mathrm{(rad)}$', yLabel=r'$\log_{10}(p_{T})\ \mathrm{(GeV)}$', scale=scale1)
    functions.plot_2dhist(hist_theta_pt_hits, f"{outdir}/theta_pt_hits_{dataset_tag}.png", f"{dataset_tag}: theta vs pt (hits)", xMin=-4, xMax=0.2, yMin=-4, yMax=0, xLabel=r'$\log_{10}(\theta)\ \mathrm{(rad)}$', yLabel=r'$\log_{10}(p_{T})\ \mathrm{(GeV)}$', scale=scale2)

    hist_z_r_mc = fIn.Get('z_r_mc')
    hist_z_r_hits = fIn.Get('z_r_hits')

    functions.plot_2dhist(hist_z_r_mc, f"{outdir}/z_r_mc_{dataset_tag}.png", f"{dataset_tag}: z vs r (MC)", xMin=-110, xMax=110, yMin=0, yMax=50, xLabel=r'z (mm)', yLabel=r'r (mm)', scale=1)
    functions.plot_2dhist(hist_z_r_hits, f"{outdir}/z_r_hits_{dataset_tag}.png", f"{dataset_tag}: z vs r (hits)", xMin=-110, xMax=110, yMin=0, yMax=50, xLabel=r'z (mm)', yLabel=r'r (mm)', scale=1)



