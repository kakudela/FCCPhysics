import os
import shutil
import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "--clean", action="store_true", help="Delete compare output directory before plotting")
args, _unknown = parser.parse_known_args()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)

index_php_src = "/home/submit/kudela/public_html/fccee/beam_background/index.php"
# Helpers
def ensure_dir(d):
    os.makedirs(d, exist_ok=True)
    try_copy_index_php(d)


def clear_dir_contents(d):
    if not d:
        return
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        return

    for name in os.listdir(d):
        p = os.path.join(d, name)
        try:
            if os.path.isdir(p) and not os.path.islink(p):
                shutil.rmtree(p)
            else:
                os.unlink(p)
        except Exception as e:
            print(f"[warn] could not remove '{p}': {e}")



def try_copy_index_php(d):
    if not os.path.isfile(index_php_src):
        return
    dst = os.path.join(d, "index.php")
    if os.path.isfile(dst):
        return
    try:
        shutil.copy(index_php_src, dst)
    except Exception:
        pass


def open_root(path):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def fetch_hist(f, hist_name):
    obj = f.Get(hist_name)
    if not obj:
        raise RuntimeError(f"hist not found: '{hist_name}' in file: {f.GetName()}")
    if not isinstance(obj, ROOT.TH1):
        raise RuntimeError(f"object '{hist_name}' is not a TH1 in file: {f.GetName()}")

    hc = obj.Clone(f"{hist_name}__clone__{abs(hash((f.GetName(), hist_name)))}")
    hc.SetDirectory(0)
    hc.SetTitle("")  # we draw titles with TLatex
    return hc


# analysis.py writes TParameter<int>("n_events") into the ROOT file, fall back to 1 if missing
def read_n_events(f):
    obj = f.Get("n_events")
    if obj and hasattr(obj, "GetVal"):
        try:
            v = int(obj.GetVal())
            return v if v > 0 else 1
        except Exception:
            pass
    return 1


def integral_with_uof(h):
    # include under/overflow
    return float(h.Integral(0, h.GetNbinsX() + 1))


def apply_norm(h, norm_mode, n_events): # ROOT hists are NOT per-event normalized by analysis.py.
    if n_events <= 0:
        n_events = 1

    # per-event scaling (always for both modes)
    h.Scale(1.0 / float(n_events))

    if norm_mode == "per_event":
        return "per_event"

    if norm_mode == "unit_area":
        integ = integral_with_uof(h)
        if integ > 0.0:
            h.Scale(1.0 / integ)
        return "unit_area_after_per_event"

    raise RuntimeError(f"unknown norm_mode: {norm_mode}")


def pick_colors(n):
    base = [
        ROOT.kBlue + 1,
        ROOT.kRed + 1,
        ROOT.kOrange + 7,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kBlack,
        ROOT.kCyan + 1,
        ROOT.kViolet + 1,
        ROOT.kTeal + 2,
    ]
    return [base[i % len(base)] for i in range(n)]


def apply_style(h, color):
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetMarkerSize(0)


def make_legend(x1, y1, x2, y2, n_entries):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    if n_entries > 6:
        leg.SetNColumns(2)
    return leg


def draw_title(canvas, title):
    canvas.cd()
    t = ROOT.TLatex()
    t.SetNDC(True)
    t.SetTextFont(42)
    t.SetTextSize(0.04)
    t.DrawLatex(0.12, 0.96, title)
    return t


def _resolve_legend_pos(plot, group_legend_pos):
    if "legend_pos" in plot and plot["legend_pos"] is not None:
        return plot["legend_pos"]
    if group_legend_pos is not None:
        return group_legend_pos
    return "tr"


def _resolve_legend_style(plot, group_legend_style):
    style = plot.get("legend_style", None)
    if style is None:
        style = group_legend_style if group_legend_style is not None else "l"
    style = str(style).strip()
    if style == "":
        style = "l"

    allowed = set("lpf")
    if any(ch not in allowed for ch in style):
        raise RuntimeError(f"invalid legend_style='{style}' (allowed chars: l, p, f)")
    return style


def overlay_1d(plot, samples, out_png, show_legend=True, group_legend_pos=None, group_legend_style=None):
    hist_name = plot["hist"]
    norm_mode = plot.get("norm", "per_event")
    title = plot["title"]
    x_title = plot["x_title"]
    y_title = plot["y_title"]

    legend_pos_resolved = _resolve_legend_pos(plot, group_legend_pos)
    legend_style_resolved = _resolve_legend_style(plot, group_legend_style)

    if legend_pos_resolved == "none":
        show_legend = False

    x_range = plot.get("x_range", None)  # (xmin,xmax) or None
    y_range = plot.get("y_range", None)  # (ymin,ymax) or None

    file_cache = {}
    hists = []  # (hist, label, path)
    any_nonzero = False

    for s in samples:
        path = s["path"]
        label = s["label"]

        if path not in file_cache:
            file_cache[path] = open_root(path)
        f = file_cache[path]

        try:
            hc = fetch_hist(f, hist_name)
        except Exception as e:
            print(f"[warn] {label}: {e}")
            continue

        n_events = read_n_events(f)
        integ_before = integral_with_uof(hc)

        norm_tag = apply_norm(hc, norm_mode, n_events)
        integ_after = integral_with_uof(hc)

        print(
            f"[info] hist='{hist_name}' sample='{label}' "
            f"n_events={n_events} entries={hc.GetEntries():.0f} "
            f"integ={integ_before:.6g} -> integ_after={integ_after:.6g} norm={norm_tag}"
        )

        if integ_after <= 0.0:
            print(f"[warn] hist '{hist_name}' has ZERO integral after norm for sample '{label}' -> will still draw (flat).")
        else:
            any_nonzero = True

        hists.append((hc, label, path))

    if not hists:
        print(f"[warn] No histograms found for '{hist_name}' in any sample -> skipping output '{out_png}'")
        for f in file_cache.values():
            f.Close()
        return

    colors = pick_colors(len(hists))
    for i, (h, _, _) in enumerate(hists):
        apply_style(h, colors[i])

    # y-range (after normalization) with optional x-range trimming
    ymax = 0.0
    for h, _, _ in hists:
        if x_range is not None:
            h.GetXaxis().SetRangeUser(float(x_range[0]), float(x_range[1]))
        m = float(h.GetMaximum())
        if m > ymax:
            ymax = m
        if x_range is not None:
            h.GetXaxis().UnZoom()

    if ymax <= 0.0:
        ymax = 1.0

    c = ROOT.TCanvas(f"c_{abs(hash((hist_name, out_png)))}", "", 900, 750)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)
    c.SetLogy(False)

    h0 = hists[0][0]

    # apply x trim (view only)
    if x_range is not None:
        h0.GetXaxis().SetRangeUser(float(x_range[0]), float(x_range[1]))

    # force y-axis to start at 0, with optional override
    if y_range is not None:
        h0.SetMinimum(float(y_range[0]))
        h0.SetMaximum(float(y_range[1]))
    else:
        h0.SetMinimum(0.0)
        h0.SetMaximum(1.25 * ymax)

    # hard-typed axis titles
    h0.GetXaxis().SetTitle(x_title)
    h0.GetYaxis().SetTitle(y_title)

    h0.Draw("hist")

    for h, _, _ in hists[1:]:
        if x_range is not None:
            h.GetXaxis().SetRangeUser(float(x_range[0]), float(x_range[1]))
        h.Draw("hist same")

    draw_title(c, title)

    if show_legend and len(hists) > 1:
        if legend_pos_resolved == "tr":
            leg = make_legend(0.62, 0.67, 0.92, 0.87, len(hists))
        elif legend_pos_resolved == "tl":
            leg = make_legend(0.14, 0.67, 0.44, 0.87, len(hists))
        elif legend_pos_resolved == "br":
            leg = make_legend(0.70, 0.16, 0.92, 0.36, len(hists))
        elif legend_pos_resolved == "bl":
            leg = make_legend(0.14, 0.14, 0.44, 0.36, len(hists))
        elif legend_pos_resolved == "tc":
            leg = make_legend(0.45, 0.70, 0.65, 0.90, len(hists))
        elif legend_pos_resolved == "cc":
            leg = make_legend(0.50, 0.43, 0.65, 0.63, len(hists))
        else:
            leg = make_legend(0.62, 0.67, 0.92, 0.87, len(hists))

        for h, label, _ in hists:
            leg.AddEntry(h, label, legend_style_resolved)
        leg.Draw()

    out_dir = os.path.dirname(out_png)
    ensure_dir(out_dir)
    c.SaveAs(out_png)
    c.Close()

    for f in file_cache.values():
        f.Close()

    if not any_nonzero:
        print(f"[warn] All samples had zero integral for '{hist_name}' -> output is flat/empty: {out_png}")
    else:
        print(f"[ok] wrote {out_png}")


def run_group(group):
    outdir = group["outdir"]
    samples = group["samples"]
    show_legend = group.get("show_legend", True)
    group_legend_pos = group.get("legend_pos", "tr")
    group_legend_style = group.get("legend_style", "l")
    plots = group["plots"]

    ensure_dir(outdir)

    for p in plots:
        out_png = os.path.join(outdir, p["out"])
        overlay_1d(
            plot=p,
            samples=samples,
            out_png=out_png,
            show_legend=show_legend,
            group_legend_pos=group_legend_pos,
            group_legend_style=group_legend_style,
        )

def main():
    outdir = "/home/submit/kudela/public_html/fccee/beam_background/ddsim/compare"
    ensure_dir(outdir)
    if args.clean:
        clear_dir_contents(outdir)
        ensure_dir(outdir)  # restore index.php

    group_mult = {
        "outdir": outdir,
        "show_legend": True,
        "legend_pos": "cc",
        "legend_style": "l",
        "samples": [
            {"label": "0T", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_ddsim_analysis.root"},
            {"label": "1T", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_1T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_1T_grids8_ddsim_analysis.root"},
            {"label": "2T", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
            {"label": "3T", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_3T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_3T_grids8_ddsim_analysis.root"},
            {"label": "4T", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_4T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_4T_grids8_ddsim_analysis.root"},
        ],
        "plots": [
            {
                "hist": "mult_per_mc_barrel_l1",
                "out": "mult_per_mc_barrel_l1_overlay.png",
                "norm": "per_event",
                "title": "",
                "x_title": "layer 1 simhits per MC",
                "y_title": "entries per event",
                "x_range": (0.0, 22.0),
                "legend_pos": "cc",
            },
            {
                "hist": "mult_per_mc_barrel_l1_dedup25um",
                "out": "mult_per_mc_barrel_l1_dedup25um_overlay.png",
                "norm": "per_event",
                "title": "",
                "x_title": "layer 1 simhits per MC (25 um merge)",
                "y_title": "entries per event",
                "x_range": (0.0, 22.0),
                "legend_pos": "cc",
            },
            {
                "hist": "mult_per_mc_barrel_l1",
                "out": "mult_per_mc_barrel_l1_overlay_unit_area.png",
                "norm": "unit_area",
                "title": "",
                "x_title": "layer 1 simhits per MC",
                "y_title": "unit area",
                "x_range": (0.0, 22.0),
                "legend_pos": "cc",
            },
            {
                "hist": "mult_per_mc_barrel_l1_dedup25um",
                "out": "mult_per_mc_barrel_l1_dedup25um_overlay_unit_area.png",
                "norm": "unit_area",
                "title": "",
                "x_title": "layer 1 simhits per MC (25 um merge)",
                "y_title": "unit area",
                "x_range": (0.0, 22.0),
                "legend_pos": "cc",
            },
        ],
    }

    group_fixes = {
        "outdir": outdir,
        "show_legend": True,
        "legend_pos": "tc",
        "legend_style": "l",
        "samples": [
            {"label": "ORIGINAL", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_noBoundary/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_noBoundary_ddsim_analysis.root"},
            {"label": "0T GRIDS8","path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_2T_ddsim/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_2T_ddsim_ddsim_analysis.root"},
            {"label": "PAIRS0",   "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs0/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
            {"label": "VTX000",   "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_vtx000/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_vtx000_ddsim_analysis.root"},
            {"label": "2T GRIDS8","path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
        ],
        "plots": [
            # barrel L1 (per-event)
            {"hist": "phi_barrel_l1",     "out": "fixes_phi_barrel_l1.png",     "norm": "per_event", "title": "", "x_title": "barrel L1 hit phi (degrees)",     "y_title": "entries per event", "legend_pos": "tc"},
            {"hist": "theta_barrel_l1",   "out": "fixes_theta_barrel_l1.png",   "norm": "per_event", "title": "", "x_title": "barrel L1 hit theta (degrees)",   "y_title": "entries per event", "legend_pos": "cc"},
            {"hist": "costheta_barrel_l1","out": "fixes_costheta_barrel_l1.png","norm": "per_event", "title": "", "x_title": "barrel L1 cos(theta)",            "y_title": "entries per event", "legend_pos": "cc"},
            {"hist": "z_barrel_l1",       "out": "fixes_z_barrel_l1.png",       "norm": "per_event", "title": "", "x_title": "barrel L1 hit z (mm)",            "y_title": "entries per event", "legend_pos": "tc"},

            # endcap L1 (per-event)
            {"hist": "phi_endcap_l1",     "out": "fixes_phi_endcap_l1.png",     "norm": "per_event", "title": "", "x_title": "endcap L1 hit phi (degrees)",     "y_title": "entries per event", "legend_pos": "tc"},
            {"hist": "theta_endcap_l1",   "out": "fixes_theta_endcap_l1.png",   "norm": "per_event", "title": "", "x_title": "endcap L1 hit theta (degrees)",   "y_title": "entries per event", "legend_pos": "cc"},
            {"hist": "costheta_endcap_l1","out": "fixes_costheta_endcap_l1.png","norm": "per_event", "title": "", "x_title": "endcap L1 cos(theta)",            "y_title": "entries per event", "legend_pos": "cc"},

            # barrel L1 (unit-area)
            {"hist": "phi_barrel_l1",     "out": "fixes_phi_barrel_l1_unit_area.png",      "norm": "unit_area", "title": "", "x_title": "barrel L1 hit phi (degrees)",   "y_title": "unit area", "legend_pos": "tc"},
            {"hist": "theta_barrel_l1",   "out": "fixes_theta_barrel_l1_unit_area.png",    "norm": "unit_area", "title": "", "x_title": "barrel L1 hit theta (degrees)", "y_title": "unit area", "legend_pos": "cc"},
            {"hist": "costheta_barrel_l1","out": "fixes_costheta_barrel_l1_unit_area.png", "norm": "unit_area", "title": "", "x_title": "barrel L1 cos(theta)",          "y_title": "unit area", "legend_pos": "cc"},
            {"hist": "z_barrel_l1",       "out": "fixes_z_barrel_l1_unit_area.png",        "norm": "unit_area", "title": "", "x_title": "barrel L1 hit z (mm)",          "y_title": "unit area", "legend_pos": "br"},

            # endcap L1 (unit-area)
            {"hist": "phi_endcap_l1",     "out": "fixes_phi_endcap_l1_unit_area.png",      "norm": "unit_area", "title": "", "x_title": "endcap L1 hit phi (degrees)",   "y_title": "unit area", "legend_pos": "br"},
            {"hist": "theta_endcap_l1",   "out": "fixes_theta_endcap_l1_unit_area.png",    "norm": "unit_area", "title": "", "x_title": "endcap L1 hit theta (degrees)", "y_title": "unit area", "legend_pos": "cc"},
            {"hist": "costheta_endcap_l1","out": "fixes_costheta_endcap_l1_unit_area.png", "norm": "unit_area", "title": "", "x_title": "endcap L1 cos(theta)",          "y_title": "unit area", "legend_pos": "cc"},
        ],
    }

    group_incidence = {
        "outdir": outdir,
        "show_legend": True,
        "legend_pos": "tr",
        "legend_style": "l",
        "samples": [
            {"label": "2T beam fields", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
            {"label": "2T no beam fields", "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_nbf_ddsim_analysis.root"},
        ],
        "plots": [
            {
                "hist": "incidence_angle_deg_barrel_l1",
                "out": "incidence_angle_deg_barrel_l1.png",
                "norm": "per_event",
                "title": "",
                "x_title": "incidence angle wrt radial normal (degrees)",
                "y_title": "entries per event",
                "legend_pos": "tr",
            },
            {
                "hist": "incidence_angle_deg_barrel_l1",
                "out": "incidence_angle_deg_barrel_l1_unit_area.png",
                "norm": "unit_area",
                "title": "",
                "x_title": "incidence angle wrt radial normal (degrees)",
                "y_title": "unit area",
                "legend_pos": "tr",
            },
            {
                "hist": "incidence_angle_phi_deg_barrel_l1",
                "out": "incidence_angle_phi_deg_barrel_l1.png",
                "norm": "per_event",
                "title": "",
                "x_title": "incidence angle (phi component) wrt radial normal (degrees)",
                "y_title": "entries per event",
                "legend_pos": "tr",
            },
            {
                "hist": "incidence_angle_phi_deg_barrel_l1",
                "out": "incidence_angle_phi_deg_barrel_l1_unit_area.png",
                "norm": "unit_area",
                "title": "",
                "x_title": "incidence angle (phi component) wrt radial normal (degrees)",
                "y_title": "unit area",
                "legend_pos": "tr",
            },
            {
                "hist": "incidence_angle_theta_deg_barrel_l1",
                "out": "incidence_angle_theta_deg_barrel_l1.png",
                "norm": "per_event",
                "title": "",
                "x_title": "incidence angle (theta component) wrt radial normal (degrees)",
                "y_title": "entries per event",
                "legend_pos": "tr",
            },
            {
                "hist": "incidence_angle_theta_deg_barrel_l1",
                "out": "incidence_angle_theta_deg_barrel_l1_unit_area.png",
                "norm": "unit_area",
                "title": "",
                "x_title": "incidence angle (theta component) wrt radial normal (degrees)",
                "y_title": "unit area",
                "legend_pos": "tr",
            },

        ]
    }

    run_group(group_mult)
    run_group(group_fixes)
    run_group(group_incidence)


if __name__ == "__main__":
    main()
