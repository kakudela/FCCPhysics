import numbers
import os
import glob
import math
import argparse
import logging
import shutil
import time
import ROOT

from functions import (ensure_index_php, ensure_public_dir, save_hist1d, save_hist2d, save_hist2d_with_box, hist_mean_and_p95,
                       sensor_area_mm2, pixel_area_mm2, write_text, filter_good_files, weed_files_on_open,)

logging.basicConfig(format="%(levelname)s: %(message)s")
log = logging.getLogger("fcclogger")
log.setLevel(logging.INFO)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(1110)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH2.SetDefaultSumw2(True)

ROOT.gInterpreter.Declare('#include "functions.h"')
ROOT.gInterpreter.Declare('#include "digitizer_simHit.h"')

index_php_source = "/home/submit/kudela/public_html/fccee/beam_background/index.php"
public_ddsim_base = "/home/submit/kudela/public_html/fccee/beam_background/ddsim"

# (input_dir, bfield_T)
# note that the detector B-field is only used for calculating the VTX acceptance box on the theta-pT plot in this script
default_input_dirs = [
    ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z512_2T_grids8", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z256_2T_grids8", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z320_2T_grids8", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z64_2T_grids8", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z32_2T_grids8", 2.0),

    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids2", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids3", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids4", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids5", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids6", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids7", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids2_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids3_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids4_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids5_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids6_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids7_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_nbf", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_noBoundary", 2.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_noBoundary", 2.0),

    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids2", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids3", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids4", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids5", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids6", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids7", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids2_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids3_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids4_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids5_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids6_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids7_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_nbf", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_theta0", 0.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0p1T_grids8", 0.1),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0p2T_grids8", 0.2),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0p3T_grids8", 0.3),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0p4T_grids8", 0.4),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0p5T_grids8", 0.5),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0p75T_grids8", 0.75),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_1T_grids8", 1.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_3T_grids8", 3.0),
    # ("/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_4T_grids8", 4.0),

    # ("/ceph/submit/data/group/fcc/kudela/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_vtx000", 2.0),
    # ("/ceph/submit/data/group/fcc/kudela/ddsim/guineapig_orig/CLD_o2_v07_2T/FCCee_Z_4IP_04may23_FCCee_Z128", 2.0),
]

# geometry selections (mm unless stated)
barrel_l1_inner_r_min_mm = 13.0
barrel_l1_inner_r_max_mm = 14.0
barrel_l1_z_extent_mm = 109.5
barrel_l1_r_mm = 14.0
box_r_mm = 13.0

endcap_l1_z_min_mm = 158.0
endcap_l1_z_inner_hi_mm = 160.0

dedup_radius_um = 25.0
dedup_radius_mm = dedup_radius_um * 1e-3

# occupancy scaling
cluster_size = 5.0
safety_factor = 3.0
pix_pitch_x_um = 25.0
pix_pitch_y_um = 25.0

l1_r_min_mm, l1_r_max_mm = 13.0, 15.0
l2_r0_mm, l3_r0_mm = 36.0, 58.0
l23_tol_mm = 4.0

bins_phi_360 = (200, 0.0, 200.0)
bins_theta_deg = (200, 0.0, 200.0)
bins_costheta = (50, -1.0, 1.0)

bins_z_barrel = (110, -120.0, 120.0)
bins_r_zoom = (200, 13.0, 15.0)

bins_edep_mev = (250, 0.0, 0.5)

bins_hits_per_evt = (500, 0.0, 500.0)
bins_mult = (22, 0.0, 22.0)

bins_abs_dz_mm = (250, 0.0, 10.0)
bins_abs_dphi_deg = (250, 0.0, 20.0)

bins_incidence_deg = (180, 0.0, 180.0)

endcap_xy_extent_cm = 12.0
bins_ec_xy_cm = (240, -endcap_xy_extent_cm, endcap_xy_extent_cm)

bins_hit_pcomp_mev_abs = (300, 0.0, 50.0)
bins_hit_pt_mev_abs = (300, 0.0, 50.0)
bins_hit_p_mev_abs = (300, 0.0, 50.0)

bins_mc_x_mm = (200, -12.0, 12.0)
bins_mc_y_mm = (200, -12.0, 12.0)
bins_mc_z_mm = (220, -35.0, 35.0)
bins_mc_r_mm = (200, 0.0, 15.0)

bins_mc_p_mev_abs = (100, 0.0, 10.0)
bins_mc_pcomp_mev_abs = (100, 0.0, 10.0)
bins_mc_e_mev = (200, 0.0, 30.0)

# mc energy in gev for log-y plot
bins_mc_e_gev = (600, 0.0, 50.0)

bins_mc_beta = (200, 0.0, 1.1)
bins_mc_beta_comp = (240, 0.0, 1.2)

log_theta_min, log_theta_max = -4.5, 0.2
log_pt_min, log_pt_max = -4.0, 0.0
bins_log_theta = (235, log_theta_min, log_theta_max)
bins_log_pt = (200, log_pt_min, log_pt_max)

bins_global_digis = (500, 0.0, 500.0)
bins_global_occ_u6 = (500, 0.0, 500.0)

occ_configs = [
    ("pix_25x25", 25, 25, True, (400, 0.0, 100.0), (100, 0.0, 10.0)),
    ("pix_50x50", 50, 50, True, (400, 0.0, 100.0), (100, 0.0, 10.0)),
    ("pix_100x100", 100, 100, True, (400, 0.0, 100.0), (100, 0.0, 10.0)),
]


def collect_files(input_dir, maxfiles):
    roots = sorted(glob.glob(os.path.join(input_dir, "*.root")))
    pairs = []
    pairs0 = []
    for f in roots:
        b = os.path.basename(f)
        if b.startswith("ddsim0_") and b.endswith(".root"):
            pairs0.append(f)
        elif b.startswith("output_") and b.endswith("sim.root"):
            pairs.append(f)

    if maxfiles > 0:
        pairs = pairs[:maxfiles]
        pairs0 = pairs0[:maxfiles]
    return pairs, pairs0


def cyclotron_pt_threshold_gev(b_t, r_mm):
    r_m = r_mm * 1e-3
    return 0.3 * b_t * r_m


def theta_cut_from_z_extent_rad(r_mm, z_half_mm):
    return math.atan2(r_mm, z_half_mm)


def tparam_i(name, value):
    return ROOT.TParameter("int")(name, int(value))


def tparam_d(name, value):
    return ROOT.TParameter("double")(name, float(value))


def write_params(tf, params):
    tf.cd()
    for p in params:
        p.Write()

def read_tparam(tf, name, default=None):
    obj = tf.Get(name)
    if not obj:
        return default
    try:
        return obj.GetVal()
    except Exception:
        try:
            return obj.GetValue()
        except Exception:
            return default

def make_plot_labels(occ_configs, pix_pitch_x_um, pix_pitch_y_um):
    def occ_dr_label(cfg_name, n_r, n_c):
        win_x_um = float(n_r) * float(pix_pitch_x_um)
        win_y_um = float(n_c) * float(pix_pitch_y_um)
        return f"{cfg_name} dR={int(round(win_x_um))}x{int(round(win_y_um))}um"

    occ_dr = {}
    for cfg_name, n_r, n_c, _, _, _ in occ_configs:
        occ_dr[cfg_name] = occ_dr_label(cfg_name, n_r, n_c)

    xlabel_1d = {
        "phi_barrel_l1": "hit phi (degrees)",
        "theta_barrel_l1": "hit theta (degrees)",
        "costheta_barrel_l1": "cos(theta)",
        "z_barrel_l1": "hit z (mm)",
        "edep_barrel_l1_mev": "edep (mev)",

        "phi_endcap_l1": "hit phi (degrees)",
        "theta_endcap_l1": "hit theta (degrees)",
        "costheta_endcap_l1": "cos(theta)",
        "edep_endcap_l1_mev": "edep (mev)",

        "simhits_barrel_l1_per_event": "simhits per event",
        "dedup25um_simhits_barrel_l1_per_event": "unique simhits per event (25 um merge)",
        "unique_mc_barrel_l1_per_event": "unique MC with >=1 hit per event",

        "simhits_endcap_l1_per_event": "simhits per event",
        "dedup25um_simhits_endcap_l1_per_event": "unique simhits per event (25 um merge)",
        "unique_mc_endcap_l1_per_event": "unique MC with >=1 hit per event",

        "r_zoom_barrel_13to15": "hit radius r (mm)",

        "mult_per_mc_barrel_l1": "simhits per mc (L1 inner barrel)",
        "mult_per_mc_barrel_l1_dedup25um": "simhits per mc (L1 inner barrel, 25 um merge)",

        "abs_dz_seq_mm": "abs(dz) consecutive hits (mm)",
        "abs_dphi_seq_deg": "abs(dphi) consecutive hits (degrees)",

        "hit_px_mev_abs_barrel_l1": "abs(px) (mev)",
        "hit_py_mev_abs_barrel_l1": "abs(py) (mev)",
        "hit_pz_mev_abs_barrel_l1": "abs(pz) (mev)",
        "hit_pt_mev_abs_barrel_l1": "abs(pt) (mev)",
        "hit_p_mev_abs_barrel_l1": "abs(p) (mev)",

        "incidence_angle_deg_barrel_l1": "incidence angle wrt radial normal (degrees)",
        "hits_per_layer_barrel_3layers": "barrel layer",

        "global_hits_total": "digis per event (L1 inner barrel)",
        "global_occ_u6": "global occupancy (x1e-6)",

        "mc_x_mm": "mc vertex x (mm)",
        "mc_y_mm": "mc vertex y (mm)",
        "mc_z_mm": "mc vertex z (mm)",
        "mc_r_mm": "mc vertex r (mm)",

        "mc_p_mev_abs": "abs(p) (mev)",
        "mc_px_mev_abs": "abs(px) (mev)",
        "mc_py_mev_abs": "abs(py) (mev)",
        "mc_pz_mev_abs": "abs(pz) (mev)",
        "mc_pt_mev_abs": "abs(pt) (mev)",
        "mc_e_mev": "mc energy (mev)",
        "mc_e_gev": "mc energy (gev)",

        "mc_beta": "beta",
        "mc_bx": "abs(beta_x)",
        "mc_by": "abs(beta_y)",
        "mc_bz": "abs(beta_z)",

        "mc_logtheta": "log10(theta_rad)",
        "mc_logpt": "log10(pt/gev)",
        "mc_logtheta_l1hit": "log10(theta_rad) (mc with L1 hit)",
        "mc_logpt_l1hit": "log10(pt/gev) (mc with L1 hit)",
    }

    for cfg_name, _, _, _, _, _ in occ_configs:
        xlabel_1d[f"local_max_hits_{cfg_name}"] = f"local max hits ({occ_dr[cfg_name]})"
        xlabel_1d[f"local_max_occ_{cfg_name}_u6"] = f"local max occupancy (x1e-6) ({occ_dr[cfg_name]})"

    xylabel_2d = {
        "z_phi_barrel_l1": ("hit z (mm)", "hit phi (degrees)"),
        "theta_phi_barrel_l1": ("hit theta (degrees)", "hit phi (degrees)"),
        "costheta_phi_barrel_l1": ("cos(theta)", "hit phi (degrees)"),

        "endcap_l1_xy_cm": ("hit x (cm)", "hit y (cm)"),

        "mc_xy_mm": ("mc vertex x (mm)", "mc vertex y (mm)"),
        "mc_xz_mm": ("mc vertex x (mm)", "mc vertex z (mm)"),
        "mc_yz_mm": ("mc vertex y (mm)", "mc vertex z (mm)"),
        "mc_zr_mm": ("mc vertex z (mm)", "mc vertex r (mm)"),

        "mc_logtheta_logpt": ("log10(theta_rad)", "log10(pt/gev)"),
        "mc_logtheta_logpt_l1hit": ("log10(theta_rad)", "log10(pt/gev)"),
    }

    def plot_title(nm):
        if "barrel_l1" in nm or nm.endswith("_barrel_l1") or "barrel_l1_" in nm:
            region = "L1 inner barrel"
        elif "endcap_l1" in nm or nm.endswith("_endcap_l1") or "endcap_l1_" in nm:
            region = "L1 endcap"
        else:
            region = ""

        if nm == "r_zoom_barrel_13to15":
            return "L1 inner barrel r zoom"
        if nm.startswith("mc_"):
            return nm

        for cfg_name, _, _, _, _, _ in occ_configs:
            if nm == f"local_max_hits_{cfg_name}":
                return f"local max hits ({occ_dr[cfg_name]})"
            if nm == f"local_max_occ_{cfg_name}_u6":
                return f"local max occupancy (x1e-6) ({occ_dr[cfg_name]})"

        if nm == "global_hits_total":
            return "global hits total (L1 inner barrel)"
        if nm == "global_occ_u6":
            return "global occupancy (x1e-6) (L1 inner barrel)"

        core = nm.replace("_", " ")
        if region:
            return f"{core} ({region})"
        return core

    return xlabel_1d, xylabel_2d, plot_title, occ_dr


def plot_from_root(out_root_vtx, vtx_dir, gen_dir, gp_tag, stream, input_dir, bfield_t):
    tf = ROOT.TFile.Open(out_root_vtx, "READ")
    if not tf or tf.IsZombie():
        log.error("cannot open root file for plotting: %s", out_root_vtx)
        return

    n_events = read_tparam(tf, "n_events", default=0)
    try:
        n_events = int(n_events)
    except Exception:
        n_events = 0
    if n_events <= 0:
        log.warning("n_events missing/invalid in %s; defaulting to 1", out_root_vtx)
        n_events = 1
    scale_per_event = 1.0 / float(n_events)

    xlabel_1d, xylabel_2d, plot_title, occ_dr = make_plot_labels(occ_configs, pix_pitch_x_um, pix_pitch_y_um)

    # --- 1D plots ---
    for nm, xlab in xlabel_1d.items():
        h = tf.Get(nm)
        if not h:
            continue

        out_dir = gen_dir if nm.startswith("mc_") else vtx_dir
        ttl = plot_title(nm)

        if nm == "hits_per_layer_barrel_3layers":
            htmp = h.Clone()
            htmp.SetDirectory(0)
            htmp.GetXaxis().SetBinLabel(1, "1")
            htmp.GetXaxis().SetBinLabel(2, "2")
            htmp.GetXaxis().SetBinLabel(3, "3")
            save_hist1d(
                htmp,
                os.path.join(out_dir, nm + ".png"),
                ttl,
                xlab,
                "entries per event",
                scale=scale_per_event,
                xbuffer=1.0,
                logy=False,
            )
            continue

        if nm == "mc_e_gev":
            # handled separately as logy
            continue

        save_hist1d(
            h,
            os.path.join(out_dir, nm + ".png"),
            ttl,
            xlab,
            "entries per event",
            scale=scale_per_event,
            xbuffer=1.0,
            logy=False,
        )

    # mc_e_gev logy
    h_egev = tf.Get("mc_e_gev")
    if h_egev:
        save_hist1d(
            h_egev,
            os.path.join(gen_dir, "mc_e_gev_logy.png"),
            "mc e gev logy",
            "mc energy (gev)",
            "entries per event",
            scale=scale_per_event,
            xbuffer=1.0,
            logy=True,
        )

    # --- 2D plots ---
    for nm, (xlab, ylab) in xylabel_2d.items():
        h = tf.Get(nm)
        if not h:
            continue
        out_dir = gen_dir if nm.startswith("mc_") else vtx_dir
        ttl = plot_title(nm)
        save_hist2d(
            h,
            os.path.join(out_dir, nm + ".png"),
            ttl,
            xlab,
            ylab,
            scale=scale_per_event,
            xbuffer=1.0,
            ybuffer=1.0,
        )

    pt_thr = cyclotron_pt_threshold_gev(bfield_t, box_r_mm)
    log_pt_thr = math.log10(max(1e-12, pt_thr))

    theta_min_rad = theta_cut_from_z_extent_rad(box_r_mm, barrel_l1_z_extent_mm)
    log_theta_cut = math.log10(max(1e-12, theta_min_rad))

    h_box = tf.Get("mc_logtheta_logpt")
    if h_box:
        save_hist2d_with_box(
            h_box,
            os.path.join(gen_dir, "mc_logtheta_logpt_with_box.png"),
            "mc logtheta logpt with box",
            "log10(theta_rad)",
            "log10(pt/gev)",
            log_theta_cut, log_pt_thr, log_theta_max, log_pt_max,
            scale=scale_per_event,
        )

    h_box_hit = tf.Get("mc_logtheta_logpt_l1hit")
    if h_box_hit:
        save_hist2d_with_box(
            h_box_hit,
            os.path.join(gen_dir, "mc_logtheta_logpt_l1hit_with_box.png"),
            "mc logtheta logpt l1hit with box",
            "log10(theta_rad)",
            "log10(pt/gev)",
            log_theta_cut, log_pt_thr, log_theta_max, log_pt_max,
            scale=scale_per_event,
        )

    ensure_index_php(vtx_dir, index_php_source)
    ensure_index_php(gen_dir, index_php_source)

    tf.Close()
    log.info("replotted from root: %s (%s %s)", out_root_vtx, gp_tag, stream)


def plot_one_stream(input_dir, stream, bfield_t):
    gp_tag = os.path.basename(os.path.normpath(input_dir))
    vtx_dir = ensure_public_dir(public_ddsim_base, "vtx", gp_tag, stream, index_php_source)
    gen_dir = ensure_public_dir(public_ddsim_base, "gen", gp_tag, stream, index_php_source)

    out_root_vtx = os.path.join(vtx_dir, f"{gp_tag}_ddsim_analysis.root")
    if not os.path.exists(out_root_vtx):
        log.warning("plot-only: missing %s -> skip", out_root_vtx)
        return

    plot_from_root(out_root_vtx, vtx_dir, gen_dir, gp_tag, stream, input_dir, bfield_t)


def run_one_stream(input_dir, files, stream, threads, bfield_t):
    if not files:
        return

    t0 = time.perf_counter()
    gp_tag = os.path.basename(os.path.normpath(input_dir))

    vtx_dir = ensure_public_dir(public_ddsim_base, "vtx", gp_tag, stream, index_php_source)
    gen_dir = ensure_public_dir(public_ddsim_base, "gen", gp_tag, stream, index_php_source)

    # prune bad files
    # note that this is only implemented because for some reason ddsim outputs are sometimes corrupted (only Katie's runs :( need to fix submit.py)
    # ideally we would not need this step at all, but better safe than sorry (comment these two lines out if Jan's runs)
    files = filter_good_files(files, treename="events")
    files = weed_files_on_open(files, treename="events")

    if not files:
        log.warning("no good root files after pruning for %s %s", gp_tag, stream)
        return

    log.info("start: %s %s files=%d", gp_tag, stream, len(files))

    tprobe0 = time.perf_counter()
    for f in files[:10]:
        t0 = time.perf_counter()
        tf = ROOT.TFile.Open(f, "READ")
        ok = (tf and not tf.IsZombie())
        if tf:
            tf.Close()
        log.info("TIMING: TFile.Open %s: %.3f s (ok=%s)",
                os.path.basename(f), time.perf_counter() - t0, ok)
    log.info("TIMING: open-probe total: %.3f s", time.perf_counter() - tprobe0)

    # now do your normal RDataFrame timing
    t_open0 = time.perf_counter()
    df = ROOT.RDataFrame("events", files)
    count_events = df.Count()
    t_open1 = time.perf_counter()
    log.info("TIMING: RDataFrame+Count booking (open/metadata): %.3f s", t_open1 - t_open0)

    # mc selection: primary e-/e+ only
    df = df.Define("mc_pdg_abs", "getMCPDGID(MCParticles, true)")
    df = df.Define("mc_gen", "getMCGeneratorStatus(MCParticles)")
    df = df.Define("mc_mask", "mc_pdg_abs==11 && mc_gen==1")
    df = df.Define("mc_sel", "MCParticles[mc_mask]")
    df = df.Define("mc_sel_idx", "ROOT::VecOps::Nonzero(mc_mask)")

    df = df.Define("mc_x", "getMCVertex_x(mc_sel)")
    df = df.Define("mc_y", "getMCVertex_y(mc_sel)")
    df = df.Define("mc_z", "getMCVertex_z(mc_sel)")
    df = df.Define("mc_r", "getMCVertex_r(mc_sel)")

    df = df.Define("mc_outside_bp_mask", "mc_r > 10.0f")
    df = df.Define("n_mc_sel_outside_bp", "static_cast<int>(ROOT::VecOps::Sum(mc_outside_bp_mask))")

    df = df.Define("mc_px_gev", "getMCMomentum_px(mc_sel)")
    df = df.Define("mc_py_gev", "getMCMomentum_py(mc_sel)")
    df = df.Define("mc_pz_gev", "getMCMomentum_pz(mc_sel)")
    df = df.Define("mc_pt_gev", "getMCMomentum_pt(mc_sel)")
    df = df.Define("mc_p_gev", "getMCMomentum_p(mc_sel)")
    df = df.Define("mc_e_gev", "getMCEnergy(mc_sel)")

    df = df.Define("mc_px_mev_abs", "vecAbs(mc_px_gev*1000.0f)")
    df = df.Define("mc_py_mev_abs", "vecAbs(mc_py_gev*1000.0f)")
    df = df.Define("mc_pz_mev_abs", "vecAbs(mc_pz_gev*1000.0f)")
    df = df.Define("mc_pt_mev_abs", "vecAbs(mc_pt_gev*1000.0f)")
    df = df.Define("mc_p_mev_abs", "vecAbs(mc_p_gev*1000.0f)")
    df = df.Define("mc_e_mev", "mc_e_gev*1000.0f")

    df = df.Define("mc_beta", "getMCBeta(mc_sel)")
    df = df.Define("mc_bx", "vecAbs(vecDiv(mc_px_gev, mc_e_gev))")
    df = df.Define("mc_by", "vecAbs(vecDiv(mc_py_gev, mc_e_gev))")
    df = df.Define("mc_bz", "vecAbs(vecDiv(mc_pz_gev, mc_e_gev))")

    df = df.Define("mc_theta_rad", "getMCMomentum_theta(mc_sel, false)")
    df = df.Define("mc_theta_log10", "log10Vec(mc_theta_rad)")
    df = df.Define("mc_pt_log10", "log10Vec(mc_pt_gev)")

    # vtx barrel hits
    df = df.Define("b_x", "getSimHitPosition_x(VertexBarrelCollection)")
    df = df.Define("b_y", "getSimHitPosition_y(VertexBarrelCollection)")
    df = df.Define("b_z", "getSimHitPosition_z(VertexBarrelCollection)")
    df = df.Define("b_r", "getSimHitPosition_r(VertexBarrelCollection)")
    df = df.Define("b_absz", "vecAbs(b_z)")

    df = df.Define("b_phi_deg", "wrapVecDeg360(getSimHitPosition_phi(VertexBarrelCollection, true))")
    df = df.Define("b_theta_deg", "getSimHitPosition_theta(VertexBarrelCollection, true)")
    df = df.Define("b_theta_rad", "getSimHitPosition_theta(VertexBarrelCollection, false)")
    df = df.Define("b_costh", "cosFromRadians(b_theta_rad)")

    df = df.Define("b_edep_gev", "getEnergyDeposition(VertexBarrelCollection)")
    df = df.Define("b_edep_mev", "b_edep_gev * 1000.0f")

    df = df.Define("b_is_secondary", "isProducedBySecondary(VertexBarrelCollection)")
    df = df.Define("b_is_primary_hit", "b_is_secondary == false")

    df = df.Define("b_mc", "getMCParticle(VertexBarrelCollection, MCParticles, _VertexBarrelCollection_particle.index)")
    df = df.Define("b_mc_gen", "getMCGeneratorStatus(b_mc)")
    df = df.Define("b_mc_pdg_abs", "getMCPDGID(b_mc, true)")
    df = df.Define("b_is_elec", "b_mc_pdg_abs==11")
    df = df.Define("b_is_mc_primary", "b_mc_gen==1")

    df = df.Define("b_l1inner_geom", f"(b_r>={barrel_l1_inner_r_min_mm}f) & (b_r<{barrel_l1_inner_r_max_mm}f) & (b_absz<={barrel_l1_z_extent_mm}f)",)
    df = df.Define("b_l1inner_mask", "b_l1inner_geom & b_is_primary_hit & b_is_elec & b_is_mc_primary")

    df = df.Define("b_l1_phi", "b_phi_deg[b_l1inner_mask]")
    df = df.Define("b_l1_theta", "b_theta_deg[b_l1inner_mask]")
    df = df.Define("b_l1_costh", "b_costh[b_l1inner_mask]")
    df = df.Define("b_l1_z", "b_z[b_l1inner_mask]")
    df = df.Define("b_l1_edep", "b_edep_mev[b_l1inner_mask]")

    df = df.Define("b_px_mev_abs_all", "vecAbs(getSimHitMomentum_x(VertexBarrelCollection)*1000.0f)")
    df = df.Define("b_py_mev_abs_all", "vecAbs(getSimHitMomentum_y(VertexBarrelCollection)*1000.0f)")
    df = df.Define("b_pz_mev_abs_all", "vecAbs(getSimHitMomentum_z(VertexBarrelCollection)*1000.0f)")
    df = df.Define("b_pt_mev_abs_all", "vecAbs(getSimHitMomentum_pt(VertexBarrelCollection)*1000.0f)")
    df = df.Define("b_p_mev_abs_all", "vecAbs(getSimHitMomentum_p(VertexBarrelCollection)*1000.0f)")

    df = df.Define("b_l1_px_mev_abs", "b_px_mev_abs_all[b_l1inner_mask]")
    df = df.Define("b_l1_py_mev_abs", "b_py_mev_abs_all[b_l1inner_mask]")
    df = df.Define("b_l1_pz_mev_abs", "b_pz_mev_abs_all[b_l1inner_mask]")
    df = df.Define("b_l1_pt_mev_abs", "b_pt_mev_abs_all[b_l1inner_mask]")
    df = df.Define("b_l1_p_mev_abs", "b_p_mev_abs_all[b_l1inner_mask]")

    df = df.Define("b_l1_px_gev_dir", "getSimHitMomentum_x(VertexBarrelCollection)[b_l1inner_mask]")
    df = df.Define("b_l1_py_gev_dir", "getSimHitMomentum_y(VertexBarrelCollection)[b_l1inner_mask]")
    df = df.Define("b_l1_pz_gev_dir", "getSimHitMomentum_z(VertexBarrelCollection)[b_l1inner_mask]")
    df = df.Define("b_l1_inc_deg", "getBarrelIncidenceAngleDeg(b_x[b_l1inner_mask], b_y[b_l1inner_mask], b_l1_px_gev_dir, b_l1_py_gev_dir, b_l1_pz_gev_dir)")

    df = df.Define("n_hits_b_l1", "static_cast<int>(b_l1_z.size())")

    # r zoom
    df = df.Define("b_rzoom_mask", f"(b_r>={l1_r_min_mm}f) & (b_r<={l1_r_max_mm}f) & (b_absz<={barrel_l1_z_extent_mm}f) & b_is_primary_hit & b_is_elec & b_is_mc_primary")
    df = df.Define("b_r_zoom", "b_r[b_rzoom_mask]")

    # 3-layer mapping for hits-per-layer plot
    df = df.Define(
        "b_layer3",
        f"ROOT::VecOps::Where((b_r>={l1_r_min_mm}f) & (b_r<{l1_r_max_mm}f), 0, "
        f"ROOT::VecOps::Where((abs(b_r-{l2_r0_mm}f)<{l23_tol_mm}f), 1, "
        f"ROOT::VecOps::Where((abs(b_r-{l3_r0_mm}f)<{l23_tol_mm}f), 2, -1)))")

    df = df.Define("b_layer3_sel", "b_layer3[b_is_primary_hit & b_is_elec & b_is_mc_primary]")
    df = df.Define("b_layer3_sel_valid", "b_layer3_sel[b_layer3_sel>=0]")

    # multiplicity on barrel l1 inner
    df = df.Define("b_mcidx_all", "_VertexBarrelCollection_particle.index")
    df = df.Define("b_l1_mcidx", "b_mcidx_all[b_l1inner_mask]")
    df = df.Define("b_l1_mult", "mcHitMultiplicity(b_l1_mcidx)")

    # MC-unique count: number of unique MC particles with >=1 hit in L1 inner (per event)
    df = df.Define("b_l1_nmc_unique", "mcUniqueCount(b_l1_mcidx)")

    # multiplicity per MC (raw) and 25 um merged versions

    df = df.Define(
        "b_l1_mult_dedup",
        f"mcHitMultiplicityMergedConsecutiveByMC(b_l1_mcidx, "
        f"b_x[b_l1inner_mask], b_y[b_l1inner_mask], b_z[b_l1inner_mask], {dedup_radius_mm}f)")

    df = df.Define(
        "b_l1_nhits_dedup",
        f"mergedHitCountConsecutiveByMC(b_l1_mcidx, "
        f"b_x[b_l1inner_mask], b_y[b_l1inner_mask], b_z[b_l1inner_mask], {dedup_radius_mm}f)")

    df = df.Define("b_l1_dz_seq", "deltaSeqSameMC(b_l1_mcidx, b_l1_z)")
    df = df.Define("b_l1_dphi_seq", "deltaPhiSeqSameMC(b_l1_mcidx, b_l1_phi)")
    df = df.Define("b_l1_abs_dz_seq", "vecAbs(b_l1_dz_seq)")
    df = df.Define("b_l1_abs_dphi_seq", "vecAbs(b_l1_dphi_seq)")

    # vtx endcap hits
    df = df.Define("ec_x", "getSimHitPosition_x(VertexEndcapCollection)")
    df = df.Define("ec_y", "getSimHitPosition_y(VertexEndcapCollection)")
    df = df.Define("ec_z", "getSimHitPosition_z(VertexEndcapCollection)")
    df = df.Define("ec_absz", "vecAbs(ec_z)")

    df = df.Define("ec_phi_deg", "wrapVecDeg360(getSimHitPosition_phi(VertexEndcapCollection, true))")
    df = df.Define("ec_theta_deg", "getSimHitPosition_theta(VertexEndcapCollection, true)")
    df = df.Define("ec_theta_rad", "getSimHitPosition_theta(VertexEndcapCollection, false)")
    df = df.Define("ec_costh", "cosFromRadians(ec_theta_rad)")

    df = df.Define("ec_edep_gev", "getEnergyDeposition(VertexEndcapCollection)")
    df = df.Define("ec_edep_mev", "ec_edep_gev * 1000.0f")

    df = df.Define("ec_is_secondary", "isProducedBySecondary(VertexEndcapCollection)")
    df = df.Define("ec_is_primary_hit", "ec_is_secondary == false")

    df = df.Define("ec_mc", "getMCParticle(VertexEndcapCollection, MCParticles, _VertexEndcapCollection_particle.index)")
    df = df.Define("ec_mc_gen", "getMCGeneratorStatus(ec_mc)")
    df = df.Define("ec_mc_pdg_abs", "getMCPDGID(ec_mc, true)")
    df = df.Define("ec_is_elec", "ec_mc_pdg_abs==11")
    df = df.Define("ec_is_mc_primary", "ec_mc_gen==1")

    df = df.Define("ec_l1_inner_geom", f"(ec_absz>={endcap_l1_z_min_mm}f) & (ec_absz<{endcap_l1_z_inner_hi_mm}f)")
    df = df.Define("ec_l1_inner_mask", "ec_l1_inner_geom & ec_is_primary_hit & ec_is_elec & ec_is_mc_primary")

    df = df.Define("ec_l1_phi", "ec_phi_deg[ec_l1_inner_mask]")
    df = df.Define("ec_l1_theta", "ec_theta_deg[ec_l1_inner_mask]")
    df = df.Define("ec_l1_costh", "ec_costh[ec_l1_inner_mask]")
    df = df.Define("ec_l1_edep", "ec_edep_mev[ec_l1_inner_mask]")

    df = df.Define("ec_l1_x_cm", "ec_x[ec_l1_inner_mask] * 0.1f")
    df = df.Define("ec_l1_y_cm", "ec_y[ec_l1_inner_mask] * 0.1f")

    df = df.Define("ec_mcidx_all", "_VertexEndcapCollection_particle.index")
    df = df.Define("ec_l1_mcidx", "ec_mcidx_all[ec_l1_inner_mask]")

    df = df.Define("ec_l1_nmc_unique", "mcUniqueCount(ec_l1_mcidx)")

    df = df.Define("n_hits_ec_l1", "static_cast<int>(ec_l1_x_cm.size())")
    df = df.Define(
        "ec_l1_nhits_dedup",
        f"mergedHitCountConsecutiveByMC(ec_l1_mcidx, "
        f"ec_x[ec_l1_inner_mask], ec_y[ec_l1_inner_mask], ec_z[ec_l1_inner_mask], {dedup_radius_mm}f)")

    # occupancy (barrel inner l1 only)
    df = df.Define(
        "b_digis",
        f"DigitizerSimHitBarrel(BarrelGeometry({barrel_l1_r_mm}f, {barrel_l1_z_extent_mm}f, {pix_pitch_x_um}f, {pix_pitch_y_um}f), "
        f"b_x[b_l1inner_mask], b_y[b_l1inner_mask], b_z[b_l1inner_mask], b_edep_gev[b_l1inner_mask], 0.0f)")
    df = df.Define("global_tot_hits", "static_cast<int>(b_digis.size())")

    sensor_a = sensor_area_mm2(barrel_l1_r_mm, barrel_l1_z_extent_mm)
    pix_a = pixel_area_mm2(pix_pitch_x_um, pix_pitch_y_um)

    df = df.Define("global_occ", f"((global_tot_hits)/{sensor_a}f)*{pix_a}f*{cluster_size}f*{safety_factor}f")
    df = df.Define("global_occ_u6", "global_occ * 1.0e6f")

    for cfg_name, n_r, n_c, use_pixels, _, _ in occ_configs:
        df = df.Define(
            f"occ_{cfg_name}",
            f"getOccupancyBarrel(BarrelGeometry({barrel_l1_r_mm}f, {barrel_l1_z_extent_mm}f, {pix_pitch_x_um}f, {pix_pitch_y_um}f), "
            f"b_digis, {n_r}, {n_c}, {str(use_pixels).lower()})")
        df = df.Define(f"occ_{cfg_name}_max_hits", f"occ_{cfg_name}[0]")

        if use_pixels:
            window_a = pix_a * n_r * n_c
        else:
            window_a = sensor_a / float(n_r * n_c)

        df = df.Define(f"occ_{cfg_name}_local_max_occ", f"((occ_{cfg_name}_max_hits)/{window_a}f)*{pix_a}f*{cluster_size}f*{safety_factor}f",)
        df = df.Define(f"occ_{cfg_name}_local_max_occ_u6", f"occ_{cfg_name}_local_max_occ * 1.0e6f")

    # mc-with-hit-on-l1-inner selection and fraction
    df = df.Define("mc_sel_has_l1hit_mask", "maskMCSelectedHasHit(mc_sel_idx, b_l1_mcidx)")
    df = df.Define("mc_sel_with_l1hit", "mc_sel[mc_sel_has_l1hit_mask]")

    df = df.Define("n_mc_sel", "static_cast<int>(mc_sel.size())")
    df = df.Define("n_mc_sel_with_l1hit", "static_cast<int>(ROOT::VecOps::Sum(mc_sel_has_l1hit_mask))")

    df = df.Define("mc_theta_rad_hit", "getMCMomentum_theta(mc_sel_with_l1hit, false)")
    df = df.Define("mc_pt_gev_hit", "getMCMomentum_pt(mc_sel_with_l1hit)")
    df = df.Define("mc_theta_log10_hit", "log10Vec(mc_theta_rad_hit)")
    df = df.Define("mc_pt_log10_hit", "log10Vec(mc_pt_gev_hit)")

    # hist booking
    h1 = {}
    h2 = {}

    def h1_book(name, bins, col):
        h1[name] = df.Histo1D((name, "", *bins), col)

    def h2_book(name, xbins, ybins, xcol, ycol):
        h2[name] = df.Histo2D((name, "", *xbins, *ybins), xcol, ycol)

    # 1d angles
    h1_book("phi_barrel_l1", bins_phi_360, "b_l1_phi")
    h1_book("theta_barrel_l1", bins_theta_deg, "b_l1_theta")
    h1_book("costheta_barrel_l1", bins_costheta, "b_l1_costh")

    h1_book("phi_endcap_l1", bins_phi_360, "ec_l1_phi")
    h1_book("theta_endcap_l1", bins_theta_deg, "ec_l1_theta")
    h1_book("costheta_endcap_l1", bins_costheta, "ec_l1_costh")

    h1_book("z_barrel_l1", bins_z_barrel, "b_l1_z")

    # 2d hitmaps
    h2_book("z_phi_barrel_l1", bins_z_barrel, bins_phi_360, "b_l1_z", "b_l1_phi")
    h2_book("theta_phi_barrel_l1", bins_theta_deg, bins_phi_360, "b_l1_theta", "b_l1_phi")
    h2_book("costheta_phi_barrel_l1", bins_costheta, bins_phi_360, "b_l1_costh", "b_l1_phi")

    # energy deposition
    h1_book("edep_barrel_l1_mev", bins_edep_mev, "b_l1_edep")
    h1_book("edep_endcap_l1_mev", bins_edep_mev, "ec_l1_edep")

    # hits per event
    h1_book("simhits_barrel_l1_per_event", bins_hits_per_evt, "n_hits_b_l1")
    h1_book("dedup25um_simhits_barrel_l1_per_event", bins_hits_per_evt, "b_l1_nhits_dedup")

    h1_book("unique_mc_barrel_l1_per_event", bins_hits_per_evt, "b_l1_nmc_unique")

    h1_book("simhits_endcap_l1_per_event", bins_hits_per_evt, "n_hits_ec_l1")
    h1_book("dedup25um_simhits_endcap_l1_per_event", bins_hits_per_evt, "ec_l1_nhits_dedup")
    h1_book("unique_mc_endcap_l1_per_event", bins_hits_per_evt, "ec_l1_nmc_unique")

    # r zoom
    h1_book("r_zoom_barrel_13to15", bins_r_zoom, "b_r_zoom")

    # multiplicity
    h1_book("mult_per_mc_barrel_l1", bins_mult, "b_l1_mult")
    h1_book("mult_per_mc_barrel_l1_dedup25um", bins_mult, "b_l1_mult_dedup")

    # abs deltas
    h1_book("abs_dz_seq_mm", bins_abs_dz_mm, "b_l1_abs_dz_seq")
    h1_book("abs_dphi_seq_deg", bins_abs_dphi_deg, "b_l1_abs_dphi_seq")

    # abs momentum in mev
    h1_book("hit_px_mev_abs_barrel_l1", bins_hit_pcomp_mev_abs, "b_l1_px_mev_abs")
    h1_book("hit_py_mev_abs_barrel_l1", bins_hit_pcomp_mev_abs, "b_l1_py_mev_abs")
    h1_book("hit_pz_mev_abs_barrel_l1", bins_hit_pcomp_mev_abs, "b_l1_pz_mev_abs")
    h1_book("hit_pt_mev_abs_barrel_l1", bins_hit_pt_mev_abs, "b_l1_pt_mev_abs")
    h1_book("hit_p_mev_abs_barrel_l1", bins_hit_p_mev_abs, "b_l1_p_mev_abs")

    # endcap xy
    h2_book("endcap_l1_xy_cm", bins_ec_xy_cm, bins_ec_xy_cm, "ec_l1_x_cm", "ec_l1_y_cm")

    # incidence
    h1_book("incidence_angle_deg_barrel_l1", bins_incidence_deg, "b_l1_inc_deg")

    # hits-per-layer
    h1_book("hits_per_layer_barrel_3layers", (3, -0.5, 2.5), "b_layer3_sel_valid")

    # global occupancy
    h1_book("global_hits_total", bins_global_digis, "global_tot_hits")
    h1_book("global_occ_u6", bins_global_occ_u6, "global_occ_u6")

    # local occupancy
    for cfg_name, n_r, n_c, _, bins_hits, bins_occ_u6 in occ_configs:
        h1_book(f"local_max_hits_{cfg_name}", bins_hits, f"occ_{cfg_name}_max_hits")
        h1_book(f"local_max_occ_{cfg_name}_u6", bins_occ_u6, f"occ_{cfg_name}_local_max_occ_u6")

    # mc plots
    h1_book("mc_x_mm", bins_mc_x_mm, "mc_x")
    h1_book("mc_y_mm", bins_mc_y_mm, "mc_y")
    h1_book("mc_z_mm", bins_mc_z_mm, "mc_z")
    h1_book("mc_r_mm", bins_mc_r_mm, "mc_r")

    h2_book("mc_xy_mm", bins_mc_x_mm, bins_mc_y_mm, "mc_x", "mc_y")
    h2_book("mc_xz_mm", bins_mc_x_mm, bins_mc_z_mm, "mc_x", "mc_z")
    h2_book("mc_yz_mm", bins_mc_y_mm, bins_mc_z_mm, "mc_y", "mc_z")
    h2_book("mc_zr_mm", bins_mc_z_mm, bins_mc_r_mm, "mc_z", "mc_r")

    h1_book("mc_p_mev_abs", bins_mc_p_mev_abs, "mc_p_mev_abs")
    h1_book("mc_px_mev_abs", bins_mc_pcomp_mev_abs, "mc_px_mev_abs")
    h1_book("mc_py_mev_abs", bins_mc_pcomp_mev_abs, "mc_py_mev_abs")
    h1_book("mc_pz_mev_abs", bins_mc_pcomp_mev_abs, "mc_pz_mev_abs")
    h1_book("mc_pt_mev_abs", bins_mc_p_mev_abs, "mc_pt_mev_abs")

    h1_book("mc_e_mev", bins_mc_e_mev, "mc_e_mev")
    h1_book("mc_e_gev", bins_mc_e_gev, "mc_e_gev")

    h1_book("mc_beta", bins_mc_beta, "mc_beta")
    h1_book("mc_bx", bins_mc_beta_comp, "mc_bx")
    h1_book("mc_by", bins_mc_beta_comp, "mc_by")
    h1_book("mc_bz", bins_mc_beta_comp, "mc_bz")

    h2_book("mc_logtheta_logpt", bins_log_theta, bins_log_pt, "mc_theta_log10", "mc_pt_log10")
    h1_book("mc_logtheta", bins_log_theta, "mc_theta_log10")
    h1_book("mc_logpt", bins_log_pt, "mc_pt_log10")

    h2_book("mc_logtheta_logpt_l1hit", bins_log_theta, bins_log_pt, "mc_theta_log10_hit", "mc_pt_log10_hit")
    h1_book("mc_logtheta_l1hit", bins_log_theta, "mc_theta_log10_hit")
    h1_book("mc_logpt_l1hit", bins_log_pt, "mc_pt_log10_hit")

    # summary actions (also executed via RunGraphs)
    mean_simhits_barrel = df.Mean("n_hits_b_l1")
    mean_unique_mc_barrel = df.Mean("b_l1_nmc_unique")
    mean_dedup25_barrel = df.Mean("b_l1_nhits_dedup")

    mean_simhits_endcap = df.Mean("n_hits_ec_l1")
    mean_unique_mc_endcap = df.Mean("ec_l1_nmc_unique")
    mean_dedup25_endcap = df.Mean("ec_l1_nhits_dedup")

    sum_mc = df.Sum("n_mc_sel")
    sum_mc_hit = df.Sum("n_mc_sel_with_l1hit")
    sum_mc_outside_bp = df.Sum("n_mc_sel_outside_bp")

    # run everything once
    actions = []
    actions.append(count_events)
    actions.extend(list(h1.values()))
    actions.extend(list(h2.values()))
    actions.extend([
        mean_simhits_barrel,
        mean_unique_mc_barrel,
        mean_dedup25_barrel,
        mean_simhits_endcap,
        mean_unique_mc_endcap,
        mean_dedup25_endcap,
        sum_mc,
        sum_mc_hit,
        sum_mc_outside_bp,
    ])
    t_run0 = time.perf_counter()
    ROOT.RDF.RunGraphs(actions)
    t_run1 = time.perf_counter()
    log.info("TIMING: RunGraphs (event loop): %.3f s", t_run1 - t_run0)

    n_events = int(count_events.GetValue())
    if n_events <= 0:
        n_events = 1
    scale_per_event = 1.0 / float(n_events)

    log.info("materialized: %s %s events=%d", gp_tag, stream, n_events)

    # write root output to vtx dir, then copy to gen dir
    out_root_vtx = os.path.join(vtx_dir, f"{gp_tag}_ddsim_analysis.root")
    out_root_gen = os.path.join(gen_dir, f"{gp_tag}_ddsim_analysis.root")

    tf = ROOT.TFile(out_root_vtx, "RECREATE")
    tf.cd()

    for name, hh in h1.items():
        hh.GetPtr().Write(name)
    for name, hh in h2.items():
        hh.GetPtr().Write(name)

    # numeric summaries
    numbers_out = {}
    numbers_out["n_events"] = int(n_events)

    numbers_out["avg_simhits_barrel_l1"] = float(mean_simhits_barrel.GetValue())
    numbers_out["avg_unique_mc_barrel_l1"] = float(mean_unique_mc_barrel.GetValue())
    numbers_out["avg_dedup25um_simhits_barrel_l1"] = float(mean_dedup25_barrel.GetValue())

    numbers_out["avg_simhits_endcap_l1"] = float(mean_simhits_endcap.GetValue())
    numbers_out["avg_unique_mc_endcap_l1"] = float(mean_unique_mc_endcap.GetValue())
    numbers_out["avg_dedup25um_simhits_endcap_l1"] = float(mean_dedup25_endcap.GetValue())

    sum_mc_v = float(sum_mc.GetValue())
    sum_mc_hit_v = float(sum_mc_hit.GetValue())
    frac_mc_with_l1hit = (sum_mc_hit_v / sum_mc_v) if sum_mc_v > 0.0 else 0.0
    numbers_out["frac_mc_sel_with_l1hit"] = float(frac_mc_with_l1hit)

    sum_mc_outside_bp_v = float(sum_mc_outside_bp.GetValue())
    frac_mc_outside_bp = (sum_mc_outside_bp_v / sum_mc_v) if sum_mc_v > 0.0 else 0.0
    numbers_out["frac_mc_sel_outside_bp_r_gt_10mm"] = float(frac_mc_outside_bp)

    # occupancy summaries from histograms (histograms are u6; we also write raw to summary.txt)
    h_occ_u6 = h1["global_occ_u6"].GetPtr()
    occ_u6_mean, occ_u6_p95 = hist_mean_and_p95(h_occ_u6)
    numbers_out["avg_global_occupancy_u6"] = float(occ_u6_mean)
    numbers_out["p95_global_occupancy_u6"] = float(occ_u6_p95)

    numbers_out["avg_global_occupancy"] = float(occ_u6_mean) * 1.0e-6
    numbers_out["p95_global_occupancy"] = float(occ_u6_p95) * 1.0e-6

    # store key params as TParameters (kept u6 in ROOT file; text summary includes raw too)
    params = [
        tparam_i("n_events", numbers_out["n_events"]),
        tparam_d("avg_simhits_barrel_l1", numbers_out["avg_simhits_barrel_l1"]),
        tparam_d("avg_unique_mc_barrel_l1", numbers_out["avg_unique_mc_barrel_l1"]),
        tparam_d("avg_dedup25um_simhits_barrel_l1", numbers_out["avg_dedup25um_simhits_barrel_l1"]),
        tparam_d("avg_simhits_endcap_l1", numbers_out["avg_simhits_endcap_l1"]),
        tparam_d("avg_unique_mc_endcap_l1", numbers_out["avg_unique_mc_endcap_l1"]),
        tparam_d("avg_dedup25um_simhits_endcap_l1", numbers_out["avg_dedup25um_simhits_endcap_l1"]),
        tparam_d("avg_global_occupancy_u6", numbers_out["avg_global_occupancy_u6"]),
        tparam_d("p95_global_occupancy_u6", numbers_out["p95_global_occupancy_u6"]),
        tparam_d("frac_mc_sel_with_l1hit", numbers_out["frac_mc_sel_with_l1hit"]),
        tparam_d("frac_mc_sel_outside_bp_r_gt_10mm", numbers_out["frac_mc_sel_outside_bp_r_gt_10mm"]),
    ]
    write_params(tf, params)

    for cfg_name, _, _, _, _, _ in occ_configs:
        h_loc_u6 = h1[f"local_max_occ_{cfg_name}_u6"].GetPtr()
        loc_u6_mean, loc_u6_p95 = hist_mean_and_p95(h_loc_u6)

        numbers_out[f"avg_local_max_occupancy_{cfg_name}_u6"] = float(loc_u6_mean)
        numbers_out[f"p95_local_max_occupancy_{cfg_name}_u6"] = float(loc_u6_p95)

        numbers_out[f"avg_local_max_occupancy_{cfg_name}"] = float(loc_u6_mean) * 1.0e-6
        numbers_out[f"p95_local_max_occupancy_{cfg_name}"] = float(loc_u6_p95) * 1.0e-6

        tparam_d(f"avg_local_max_occupancy_{cfg_name}_u6", float(loc_u6_mean)).Write()
        tparam_d(f"p95_local_max_occupancy_{cfg_name}_u6", float(loc_u6_p95)).Write()

    tf.Close()

    # copy root to gen dir
    try:
        shutil.copy2(out_root_vtx, out_root_gen)
    except Exception:
        pass

    xlabel_1d, xylabel_2d, plot_title, occ_dr = make_plot_labels(occ_configs, pix_pitch_x_um, pix_pitch_y_um)

    # save 1D plots
    for nm, hh in h1.items():
        out_dir = gen_dir if nm.startswith("mc_") else vtx_dir
        xlab = xlabel_1d.get(nm, nm)
        ttl = plot_title(nm)

        # hits per layer: relabel bins as 1,2,3
        if nm == "hits_per_layer_barrel_3layers":
            htmp = hh.GetPtr().Clone()
            htmp.SetDirectory(0)
            htmp.GetXaxis().SetBinLabel(1, "1")
            htmp.GetXaxis().SetBinLabel(2, "2")
            htmp.GetXaxis().SetBinLabel(3, "3")
            save_hist1d(
                htmp,
                os.path.join(out_dir, nm + ".png"),
                ttl,
                xlab,
                "entries per event",
                scale=scale_per_event,
                xbuffer=1.0,
                logy=False,
            )
            continue

        # mc_e_gev is saved separately as logy
        if nm == "mc_e_gev":
            continue

        save_hist1d(hh.GetPtr(), os.path.join(out_dir, nm + ".png"), ttl, xlab, "entries per event", scale=scale_per_event, xbuffer=1.0, logy=False,)

    # Save mc energy gev as logy
    save_hist1d(h1["mc_e_gev"].GetPtr(), os.path.join(gen_dir, "mc_e_gev_logy.png"), "mc e gev logy", "mc energy (gev)", "entries per event", scale=scale_per_event, xbuffer=1.0, logy=True,)

    # save 2D plots
    for nm, hh in h2.items():
        out_dir = gen_dir if nm.startswith("mc_") else vtx_dir
        xlab, ylab = xylabel_2d.get(nm, ("", ""))
        ttl = plot_title(nm)
        save_hist2d(hh.GetPtr(), os.path.join(out_dir, nm + ".png"), ttl, xlab, ylab, scale=scale_per_event, xbuffer=1.0, ybuffer=1.0,)

    # box plots
    pt_thr = cyclotron_pt_threshold_gev(bfield_t, box_r_mm)
    log_pt_thr = math.log10(max(1e-12, pt_thr))

    theta_min_rad = theta_cut_from_z_extent_rad(box_r_mm, barrel_l1_z_extent_mm)
    log_theta_cut = math.log10(max(1e-12, theta_min_rad))

    save_hist2d_with_box(h2["mc_logtheta_logpt"].GetPtr(), os.path.join(gen_dir, "mc_logtheta_logpt_with_box.png"),"mc logtheta logpt with box", "log10(theta_rad)", "log10(pt/gev)", log_theta_cut, log_pt_thr, log_theta_max, log_pt_max, scale=scale_per_event,)
    save_hist2d_with_box(h2["mc_logtheta_logpt_l1hit"].GetPtr(), os.path.join(gen_dir, "mc_logtheta_logpt_l1hit_with_box.png"), "mc logtheta logpt l1hit with box", "log10(theta_rad)", "log10(pt/gev)", log_theta_cut, log_pt_thr, log_theta_max, log_pt_max, scale=scale_per_event,)

    # summary of vtx stuff stored in the html in this text file
    summary_lines = [
        f"input_dir: {input_dir}",
        f"tag: {gp_tag}",
        f"stream: {stream}",
        f"output_root_vtx: {out_root_vtx}",
        f"output_root_gen: {out_root_gen}",
        f"dedup_radius_um: {dedup_radius_um}",
        f"cluster_size: {cluster_size}",
        f"safety_factor: {safety_factor}",
        f"pix_pitch_um: {pix_pitch_x_um}x{pix_pitch_y_um}",
        "",
        "notes:",
        "  dedup25um_* means consecutive-hit merge within same MC using 25 um distance",
        "  unique_mc_* means count of unique MC particle IDs producing >=1 hit in the region",
        "",
        "occupancy_configs:",
    ]
    for cfg_name, n_r, n_c, use_pixels, _, _ in occ_configs:
        summary_lines.append(f"  {cfg_name}: n_r={n_r} n_c={n_c} use_pixels={use_pixels} window={occ_dr[cfg_name]}")
    summary_lines.append("")

    # Print all numbers (raw + u6 are both present for occupancy keys)
    for k in sorted(numbers_out.keys()):
        summary_lines.append(f"{k}: {numbers_out[k]}")

    write_text(os.path.join(vtx_dir, "summary.txt"), summary_lines)

    ensure_index_php(vtx_dir, index_php_source)
    ensure_index_php(gen_dir, index_php_source)

    log.info("wrote vtx: %s", vtx_dir)
    log.info("wrote gen: %s", gen_dir)

    t1 = time.perf_counter()
    log.info("done: %s %s seconds=%.3f", gp_tag, stream, (t1 - t0))


def run_input_dir(input_dir, threads, maxfiles, bfield_t, noPairs0=False, plot=False):
    if plot:
        # Only replot from already-written root outputs
        plot_one_stream(input_dir, "pairs", bfield_t)
        if not noPairs0:
            plot_one_stream(input_dir, "pairs0", bfield_t)
        else:
            log.info("flag --noPairs0 set -> skip pairs0 for %s", input_dir)
        return

    # normal analysis mode

    t_ls0 = time.perf_counter()
    pairs, pairs0 = collect_files(input_dir, maxfiles)
    t_ls1 = time.perf_counter()
    log.info("TIMING: collect_files (glob/listdir): %.3f s", t_ls1 - t_ls0)

    if pairs:
        run_one_stream(input_dir, pairs, "pairs", threads, bfield_t)
    else:
        log.warning("no ddsim_*.root in %s -> skip pairs", input_dir)

    if noPairs0:
        log.info("flag --noPairs0 set -> skip pairs0 for %s", input_dir)
        return

    if pairs0:
        run_one_stream(input_dir, pairs0, "pairs0", threads, bfield_t)
    else:
        log.info("no ddsim0_*.root in %s -> skip pairs0", input_dir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--maxfiles", type=int, default=0)
    ap.add_argument("--input-dir", default="", help="optional single input dir; if empty uses default list")
    ap.add_argument("--noPairs0", action="store_true", help="skip the pairs0 stream entirely")
    ap.add_argument("--plot", action="store_true", help="Plot only: read existing *_ddsim_analysis.root and regenerate PNGs (no event loop)")
    args = ap.parse_args()

    if args.threads and args.threads > 0:
        try:
            ROOT.ROOT.EnableImplicitMT(int(args.threads))
            log.info("Enabled ROOT implicit MT with %d threads", args.threads)
        except Exception as e:
            log.warning("Could not enable ROOT implicit MT: %s", e)

    if args.input_dir:
        datasets = [(args.input_dir, 2.0)]  # default for ad-hoc single run
    else:
        datasets = list(default_input_dirs)

    if not datasets:
        raise SystemExit("no input dirs provided")

    for d, bfield_t in datasets:
        run_input_dir(d, args.threads, args.maxfiles, bfield_t,
                  noPairs0=args.noPairs0,
                  plot=args.plot)

if __name__ == "__main__":
    main()
