import os
import csv
from datetime import datetime
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

csv_path = "/ceph/submit/data/group/fcc/kudela/gp_timing_studies/timings.csv"
out_dir  = "/home/submit/kudela/public_html/fccee/beam_background/gp_timing"

acc_dat_local = "acc.dat"

cut_x_multi = 150.0
cut_y_multi = 60.0
cut_z_multi = 2.0

# Volume axis scaling for plots
# Here we use cm^3 to keep numbers reasonable:
# 1 cm^3 = 1000 mm^3
vol_scale = 1.0 / 1000.0
vol_unit  = "cm^{3}"


def parse_ts(ts: str):
    return datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S")


def read_sigmas_from_acc_dat(path):
    if not os.path.exists(path):
        raise RuntimeError(f"Could not find {path} (needed to read sigma_x/y/z).")

    sigma_x = None
    sigma_y = None
    sigma_z = None

    with open(path, "r") as f:
        for line in f:
            s = line.strip().replace(" ", "")
            if s.startswith("sigma_x="):
                sigma_x = float(s.split("=", 1)[1].rstrip(";"))
            elif s.startswith("sigma_y="):
                sigma_y = float(s.split("=", 1)[1].rstrip(";"))
            elif s.startswith("sigma_z="):
                sigma_z = float(s.split("=", 1)[1].rstrip(";"))

    if sigma_x is None or sigma_y is None or sigma_z is None:
        raise RuntimeError("Failed to parse sigma_x, sigma_y, sigma_z from acc.dat")

    return sigma_x, sigma_y, sigma_z


def compute_outer_grid_half_extents_mm(grids, n_x, n_y, n_z, sigma_x_nm, sigma_y_nm, sigma_z_um):
    """
    Returns (X_mm, Y_mm, Z_mm) half extents for the outer tracking grid.

    IMPORTANT UNITS (based on your observed geometry):
      sigma_x, sigma_y are in nm
      sigma_z is in um

    This makes cut_z = 2*sigma_z ~ 30 mm for sigma_z=15200 um, matching z ~ [-30, +30] mm.
    """

    # Convert sigmas into mm
    sigma_x_mm = sigma_x_nm * 1.0e-6   # nm -> mm
    sigma_y_mm = sigma_y_nm * 1.0e-6   # nm -> mm
    sigma_z_mm = sigma_z_um * 1.0e-3   # um -> mm

    # Cuts in mm
    cut_x = cut_x_multi * sigma_x_mm
    cut_y = cut_y_multi * sigma_y_mm
    cut_z = cut_z_multi * sigma_z_mm

    # Primary grid half extents (mm)
    X0 = ((n_x - 2.0) / n_x) * cut_x
    Y0 = ((n_y - 2.0) / n_y) * cut_y
    Z0 = cut_z  # z does not change with grids

    # Secondary grids (half extents in mm)
    sec = {}

    x0 = X0
    y0 = Y0
    sec[1] = (x0, y0)

    x1 = x0 * 2.0
    y1 = y0 * 2.0
    sec[2] = (x1, y1)

    ratio = x0 / y0
    x2 = x1
    y2 = y1 * (ratio ** (1.0 / 3.0))
    sec[3] = (x2, y2)

    x3 = x1
    y3 = y1 * (ratio ** (2.0 / 3.0))
    sec[4] = (x3, y3)

    x4 = x1
    y4 = x1
    sec[5] = (x4, y4)

    x5 = x1 * 2.0
    y5 = x1 * 2.0
    sec[6] = (x5, y5)

    x6 = x1 * 3.0
    y6 = x1 * 3.0
    sec[7] = (x6, y6)

    x7 = x1 * 6.0
    y7 = x1 * 6.0
    sec[8] = (x7, y7)

    if grids <= 0:
        X, Y = sec[1]
    else:
        if grids not in sec:
            raise RuntimeError(f"grids={grids} not in [0..8] for this parameterization")
        X, Y = sec[grids]

    return X, Y, Z0


def make_graph(points, name, xtitle, ytitle):
    gr = ROOT.TGraph(len(points))
    gr.SetName(name)
    for i, (x, y) in enumerate(points):
        gr.SetPoint(i, float(x), float(y))

    gr.GetXaxis().SetTitle(xtitle)
    gr.GetYaxis().SetTitle(ytitle)
    gr.GetXaxis().CenterTitle(True)
    gr.GetYaxis().CenterTitle(True)

    gr.SetLineWidth(2)
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.1)
    return gr


def set_yrange_from_points(graph, points, ymin_pad_frac=0.05, ymax_pad_frac=0.15):
    if not points:
        return
    ys = [p[1] for p in points]
    y_min = min(ys)
    y_max = max(ys)
    span = (y_max - y_min) if (y_max > y_min) else (y_max if y_max > 0 else 1.0)
    graph.GetYaxis().SetRangeUser(max(0.0, y_min - ymin_pad_frac * span), y_max + ymax_pad_frac * span)


def draw_single(graph, points, outpath):
    c = ROOT.TCanvas("c", "c", 900, 700)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    set_yrange_from_points(graph, points)
    graph.Draw("ALP")
    c.SaveAs(outpath)
    del c


def draw_overlay(g_no, pts_no, g_tr, pts_tr, outpath, leg_x1=0.60, leg_y1=0.2, leg_x2=0.90, leg_y2=0.36):
    c = ROOT.TCanvas("c_ov", "c_ov", 900, 700)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    g_no.SetMarkerStyle(20)
    g_no.SetLineStyle(1)

    g_tr.SetMarkerStyle(21)
    g_tr.SetLineStyle(2)

    all_pts = pts_no + pts_tr
    set_yrange_from_points(g_no, all_pts)

    g_no.Draw("ALP")
    g_tr.Draw("LP SAME")

    leg = ROOT.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(g_tr, "Tracking (track_pairs = 1)", "lp")
    leg.AddEntry(g_no, "No tracking (track_pairs = 0)", "lp")
    leg.Draw()

    c.SaveAs(outpath)
    del c


# -------------------------
# Read timings.csv (latest per (grids, track_pairs))
# -------------------------
if not os.path.exists(csv_path):
    raise RuntimeError(f"CSV not found: {csv_path}")

latest = {}
all_real_s = []

with open(csv_path, "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        grids = int(row["grids"])
        tp    = int(row["track_pairs"])
        ts    = parse_ts(row["timestamp"])
        real_s = float(row["real_s"])
        all_real_s.append(real_s)

        key = (grids, tp)
        if key not in latest or ts > latest[key][0]:
            latest[key] = (ts, real_s)

if not all_real_s:
    raise RuntimeError("No timing rows found in CSV.")

max_real_s = max(all_real_s)
if max_real_s >= 3600.0:
    tscale = 1.0 / 3600.0
    unit_label = "hours"
else:
    tscale = 1.0 / 60.0
    unit_label = "minutes"

ytitle = f"Runtime per BX ({unit_label})"

# Compute geometry meta per grid
sigma_x, sigma_y, sigma_z = read_sigmas_from_acc_dat(acc_dat_local)

n_x = 128
n_y = 128
n_z = 128

grid_meta = {}
for g in range(0, 9):
    X_mm, Y_mm, Z_mm = compute_outer_grid_half_extents_mm(g, n_x, n_y, n_z, sigma_x, sigma_y, sigma_z)

    vol_mm3 = (2.0 * X_mm) * (2.0 * Y_mm) * (2.0 * Z_mm)

    dx = (2.0 * X_mm) / n_x
    dy = (2.0 * Y_mm) / n_y
    dz = (2.0 * Z_mm) / n_z
    voxel_mm3 = dx * dy * dz

    grid_meta[g] = {
        "X_mm": X_mm,
        "Y_mm": Y_mm,
        "Z_mm": Z_mm,
        "vol_mm3": vol_mm3,
        "voxel_mm3": voxel_mm3,
    }

print(f"[CHECK] Z half-extent is constant: Z=±{grid_meta[1]['Z_mm']:.3f} mm (full thickness {2*grid_meta[1]['Z_mm']:.3f} mm)")


# Helpers: build timing point lists
def build_points_runtime_vs_grids():
    pts_no = []
    pts_tr = []
    for (grids, tp), (_ts, real_s) in latest.items():
        y = real_s * tscale
        x = grids
        if tp == 0:
            pts_no.append((x, y))
        elif tp == 1:
            pts_tr.append((x, y))
    pts_no.sort(key=lambda p: p[0])
    pts_tr.sort(key=lambda p: p[0])
    return pts_no, pts_tr


def build_points_runtime_vs_volume():
    pts_no = []
    pts_tr = []
    for (grids, tp), (_ts, real_s) in latest.items():
        y = real_s * tscale
        x = grid_meta[grids]["vol_mm3"] * vol_scale
        if tp == 0:
            pts_no.append((x, y))
        elif tp == 1:
            pts_tr.append((x, y))
    pts_no.sort(key=lambda p: p[0])
    pts_tr.sort(key=lambda p: p[0])
    return pts_no, pts_tr


# plot shell volume and voxel volume vs grid number (assumes the )
def build_points_shell_volume_vs_grids():
    pts = []
    for g in range(1, 9):
        v_shell_mm3 = grid_meta[g]["vol_mm3"] - grid_meta[g - 1]["vol_mm3"]
        v_shell = v_shell_mm3 * vol_scale  # to cm^3
        pts.append((g, v_shell))
    return pts


def build_points_voxel_volume_vs_grids():
    pts = []
    for g in range(0, 9):
        vvox_mm3 = grid_meta[g]["voxel_mm3"]
        pts.append((g, vvox_mm3))
    return pts


def draw_simple(points, name, xtitle, ytitle_local, outpath):
    gr = make_graph(points, name, xtitle, ytitle_local)
    draw_single(gr, points, outpath)


def draw_incremental_runtime_per_grid(outpath):
    # compute deltaT(g) = T(g)-T(g-1) for track_pairs=0 and 1 separately
    deltas_no = []
    deltas_tr = []

    # build map grids -> runtime minutes
    t_no = {}
    t_tr = {}
    for (g, tp), (_ts, real_s) in latest.items():
        if tp == 0:
            t_no[g] = real_s * tscale
        elif tp == 1:
            t_tr[g] = real_s * tscale

    for g in range(1, 9):
        if g in t_no and (g - 1) in t_no:
            deltas_no.append((g, t_no[g] - t_no[g - 1]))
        if g in t_tr and (g - 1) in t_tr:
            deltas_tr.append((g, t_tr[g] - t_tr[g - 1]))

    g_no = make_graph(deltas_no, "g_deltaT_no", "Number of grids", f"Incremental runtime per added grid ({unit_label})")
    g_tr = make_graph(deltas_tr, "g_deltaT_tr", "Number of grids", f"Incremental runtime per added grid ({unit_label})")

    c = ROOT.TCanvas("c_dT", "c_dT", 900, 700)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    # style
    g_no.SetMarkerStyle(20)
    g_no.SetLineStyle(1)
    g_tr.SetMarkerStyle(21)
    g_tr.SetLineStyle(2)

    all_pts = deltas_no + deltas_tr
    set_yrange_from_points(g_no, all_pts)

    g_no.Draw("ALP")
    g_tr.Draw("LP SAME")

    leg = ROOT.TLegend(0.55, 0.20, 0.88, 0.36)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(g_tr, "Tracking (track_pairs = 1)", "lp")
    leg.AddEntry(g_no, "No tracking (track_pairs = 0)", "lp")
    leg.Draw()

    c.SaveAs(outpath)
    del c


def draw_incremental_runtime_per_shell_volume(outpath):
    # compute (T(g)-T(g-1)) / (V(g)-V(g-1)) for g>=1
    pts_no = []
    pts_tr = []

    t_no = {}
    t_tr = {}
    for (g, tp), (_ts, real_s) in latest.items():
        if tp == 0:
            t_no[g] = real_s * tscale
        elif tp == 1:
            t_tr[g] = real_s * tscale

    for g in range(1, 9):
        dv_mm3 = grid_meta[g]["vol_mm3"] - grid_meta[g - 1]["vol_mm3"]
        dv = dv_mm3 * vol_scale  # cm^3
        if dv <= 0:
            continue

        if g in t_no and (g - 1) in t_no:
            dt = t_no[g] - t_no[g - 1]
            pts_no.append((g, dt / dv))

        if g in t_tr and (g - 1) in t_tr:
            dt = t_tr[g] - t_tr[g - 1]
            pts_tr.append((g, dt / dv))

    g_no = make_graph(pts_no, "g_dt_dv_no", "Number of grids", f"Incremental runtime per incremental volume ({unit_label} / {vol_unit})")
    g_tr = make_graph(pts_tr, "g_dt_dv_tr", "Number of grids", f"Incremental runtime per incremental volume ({unit_label} / {vol_unit})")

    c = ROOT.TCanvas("c_dTdV", "c_dTdV", 900, 700)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.06)

    # style
    g_no.SetMarkerStyle(20)
    g_no.SetLineStyle(1)
    g_tr.SetMarkerStyle(21)
    g_tr.SetLineStyle(2)

    all_pts = pts_no + pts_tr
    set_yrange_from_points(g_no, all_pts)

    g_no.Draw("ALP")
    g_tr.Draw("LP SAME")

    leg = ROOT.TLegend(0.55, 0.20, 0.88, 0.36)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(g_tr, "Tracking (track_pairs = 1)", "lp")
    leg.AddEntry(g_no, "No tracking (track_pairs = 0)", "lp")
    leg.Draw()

    c.SaveAs(outpath)
    del c


# make output directory and write plots
os.makedirs(out_dir, exist_ok=True)

# existing style plots
pts_no, pts_tr = build_points_runtime_vs_grids()
g_no = make_graph(pts_no, "g_no_vs_grids", "Number of grids", ytitle)
g_tr = make_graph(pts_tr, "g_tr_vs_grids", "Number of grids", ytitle)
draw_overlay(g_no, pts_no, g_tr, pts_tr, os.path.join(out_dir, "timing_overlay_vs_grids.png"))

pts_no, pts_tr = build_points_runtime_vs_volume()
g_no = make_graph(pts_no, "g_no_vs_vol", f"Outer tracking region volume ({vol_unit})", ytitle)
g_tr = make_graph(pts_tr, "g_tr_vs_vol", f"Outer tracking region volume ({vol_unit})", ytitle)
draw_overlay(g_no, pts_no, g_tr, pts_tr, os.path.join(out_dir, "timing_overlay_vs_volume.png"))

# incremental plots
draw_incremental_runtime_per_shell_volume(os.path.join(out_dir, "incremental_runtime_per_incremental_volume.png"))
draw_incremental_runtime_per_grid(os.path.join(out_dir, "incremental_runtime_per_added_grid.png"))

shell_pts = build_points_shell_volume_vs_grids()
draw_simple(shell_pts, "g_shell_volume", "Number of grids", f"Shell volume added: V(g)-V(g-1) ({vol_unit})", os.path.join(out_dir, "shell_volume_added_vs_grids.png"))

voxel_pts = build_points_voxel_volume_vs_grids()
draw_simple(voxel_pts, "g_voxel_volume", "Number of grids", "Outer grid voxel volume (mm^{3})", os.path.join(out_dir, "outer_grid_voxel_volume_vs_grids.png"))

print(f"[OK] Wrote plots to: {out_dir}/")
print(f"[OK] Using runtime units: {unit_label} (max real_s={max_real_s:.2f} s)")
print("[OK] New geometry plots created:")
print("  - shell_volume_added_vs_grids.png")
print("  - outer_grid_voxel_volume_vs_grids.png")
