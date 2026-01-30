import os
import shutil
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

index_php_src = "/home/submit/kudela/public_html/fccee/beam_background/index.php"

# helpers:

def ensure_dir(d: str):
    os.makedirs(d, exist_ok=True)
    try_copy_index_php(d)


def try_copy_index_php(d: str):
    if not os.path.isfile(index_php_src):
        return
    dst = os.path.join(d, "index.php")
    if os.path.isfile(dst):
        return
    try:
        shutil.copy(index_php_src, dst)
    except Exception:
        pass


def open_root(path: str) -> ROOT.TFile:
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def read_scalar(f: ROOT.TFile, key: str) -> float:
    obj = f.Get(key)
    if not obj:
        raise RuntimeError(f"parameter not found: '{key}' in file: {f.GetName()}")

    # TParameter<T>
    if hasattr(obj, "GetVal"):
        try:
            return float(obj.GetVal())
        except Exception:
            pass

    # fallback: number stored in TNamed title/name
    if isinstance(obj, ROOT.TNamed):
        try:
            return float(obj.GetTitle())
        except Exception:
            pass
        try:
            return float(obj.GetName())
        except Exception:
            pass

    raise RuntimeError(f"parameter '{key}' not readable as scalar (type={obj.ClassName()}) in file: {f.GetName()}")


def draw_title(canvas: ROOT.TCanvas, title: str):
    if not title:
        return None
    canvas.cd()
    t = ROOT.TLatex()
    t.SetNDC(True)
    t.SetTextFont(42)
    t.SetTextSize(0.04)
    t.DrawLatex(0.12, 0.96, title)
    return t


def _auto_range(vals, pad_frac=0.15, min_top_frac=0.20, abs_min_pad=0.0):
    if not vals:
        return (0.0, 1.0)

    vmin = min(vals)
    vmax = max(vals)

    # Handle completely empty/degenerate cases safely
    if math.isfinite(vmax) is False or math.isfinite(vmin) is False:
        return (0.0, 1.0)

    span = vmax - vmin

    pad_from_span = span * float(pad_frac)
    pad_from_level = abs(vmax) * float(min_top_frac)
    pad = max(pad_from_span, pad_from_level, float(abs_min_pad))

    # If everything is exactly 0 and pad is 0, force something sane
    if vmax == 0.0 and pad == 0.0:
        pad = 1.0

    # ymin behavior: start at 0
    if vmin >= 0.0:
        ymin = 0.0
        ymax = vmax + pad
        if ymax <= 0.0:
            ymax = 1.0
        return (ymin, ymax)

    # If negatives ever appear, keep them visible with padding
    ymin = vmin - pad
    ymax = vmax + pad
    if math.isclose(ymin, ymax):
        ymax = ymin + 1.0
    return (ymin, ymax)


def _auto_xrange(xs, pad_frac=0.05):
    xmin = min(xs)
    xmax = max(xs)
    if math.isclose(xmin, xmax):
        return (xmin - 1.0, xmax + 1.0)
    pad = pad_frac * (xmax - xmin)
    return (xmin - pad, xmax + pad)


def plot_param_vs_x(
    *,
    param: str,
    y_title: str,
    out_png: str,
    title: str,
    x_title: str,
    points: list,   # list of dicts: {"x":..., "path":...}
    x_range=None,
    y_range=None,
    connect=False,
    marker_style=20,
    marker_size=1.2,
):
    # Read points
    xs, ys = [], []
    file_cache = {}

    for pt in points:
        x = float(pt["x"])
        path = pt["path"]

        if path not in file_cache:
            try:
                file_cache[path] = open_root(path)
            except Exception as e:
                print(f"[warn] could not open '{path}': {e}")
                continue

        f = file_cache[path]
        try:
            y = read_scalar(f, param)
        except Exception as e:
            print(f"[warn] x={x:g} param='{param}' in '{path}': {e}")
            continue

        xs.append(x)
        ys.append(float(y))
        print(f"[info] param='{param}' x={x:g} y={y:g} file='{path}'")

    for f in file_cache.values():
        f.Close()

    if not xs:
        print(f"[warn] No valid points for param '{param}' -> skipping '{out_png}'")
        return

    # Sort by x
    pts = sorted(zip(xs, ys), key=lambda t: t[0])
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]

    gr = ROOT.TGraph(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        gr.SetPoint(i, x, y)

    gr.SetMarkerStyle(int(marker_style))
    gr.SetMarkerSize(float(marker_size))
    gr.SetLineWidth(2)

    # Ranges
    if x_range is None:
        xmin, xmax = _auto_xrange(xs)
    else:
        xmin, xmax = float(x_range[0]), float(x_range[1])

    if y_range is None:
        # Always ymin=0 when possible, always give headroom
        ymin, ymax = _auto_range(ys, pad_frac=0.15, min_top_frac=0.08, abs_min_pad=0.0)
    else:
        ymin, ymax = float(y_range[0]), float(y_range[1])

    c = ROOT.TCanvas(f"c_{abs(hash((param, out_png)))}", "", 900, 750)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)
    c.SetLogy(False)

    frame = c.DrawFrame(xmin, ymin, xmax, ymax)
    frame.GetXaxis().SetTitle(x_title)
    frame.GetYaxis().SetTitle(y_title)

    gr.Draw("PL" if connect else "P")
    draw_title(c, title)

    ensure_dir(os.path.dirname(out_png))
    c.SaveAs(out_png)
    c.Close()

    print(f"[ok] wrote {out_png}")


def run_study(study: dict, param_defs: list):
    outdir = study["outdir"]
    x_title = study["x_title"]
    points = study["points"]

    title_prefix = study.get("title_prefix", "")
    title_suffix = study.get("title_suffix", "")
    x_range = study.get("x_range", None)

    connect = bool(study.get("connect", False))
    marker_style = study.get("marker_style", 20)
    marker_size = study.get("marker_size", 1.2)

    y_overrides = study.get("y_overrides", {})  # per-param overrides

    ensure_dir(outdir)

    for pd in param_defs:
        param = pd["param"]
        out_name = pd["out"]
        y_title = pd["y_title"]
        base_title = pd.get("title", "")

        # Build plot title for this study (optional, but useful)
        title_parts = []
        if title_prefix:
            title_parts.append(title_prefix)
        if base_title:
            title_parts.append(base_title)
        if title_suffix:
            title_parts.append(title_suffix)
        title = "  ".join([t for t in title_parts if t])

        y_range = None
        if param in y_overrides and "y_range" in y_overrides[param]:
            y_range = y_overrides[param]["y_range"]
        elif "y_range" in pd:
            y_range = pd["y_range"]

        out_png = os.path.join(outdir, out_name)
        plot_param_vs_x(
            param=param,
            y_title=y_title,
            out_png=out_png,
            title=title,
            x_title=x_title,
            points=points,
            x_range=x_range,
            y_range=y_range,
            connect=connect,
            marker_style=marker_style,
            marker_size=marker_size,
        )


def main():
    # Global parameter definitions (shared across all studies)
    # This is what ensures y-axis titles are identical across studies
    PARAM_DEFS = [
        {"param": "avg_simhits_barrel_l1",     "out": "avg_simhits_barrel_l1.png",     "y_title": "avg simhits (barrel L1)"},
        {"param": "avg_unique_mc_barrel_l1",  "out": "avg_unique_mc_barrel_l1.png",  "y_title": "avg unique MC (barrel L1)"},
        {"param": "avg_simhits_endcap_l1",    "out": "avg_simhits_endcap_l1.png",    "y_title": "avg simhits (endcap L1)"},
        {"param": "avg_unique_mc_endcap_l1",  "out": "avg_unique_mc_endcap_l1.png",  "y_title": "avg unique MC (endcap L1)"},
        {"param": "avg_global_occupancy_u6",  "out": "avg_global_occupancy_u6.png",  "y_title": "avg global occupancy (x 1e-6)"},
        {"param": "p95_global_occupancy_u6",  "out": "p95_global_occupancy_u6.png",  "y_title": "p95 global occupancy (x 1e-6)"},
        {
            "param": "frac_mc_sel_with_l1hit",
            "out": "frac_mc_sel_with_l1hit.png",
            "y_title": "fraction of selected MC w/ >/= L1 hit",
            "y_range": (0.0, 0.25),
        },
    ]

    studies = [
        {
            "name": "granularity_scan_2T",
            "outdir": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/occupancy_studies/grid_granularity/2T",
            "x_title": "number of grid divisions per axis (n)",
            "connect": True,
            "points": [
                {"x": 32, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z32_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z32_2T_grids8_ddsim_analysis.root"},
                {"x": 64, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z64_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z64_2T_grids8_ddsim_analysis.root"},
                {"x": 128, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
                {"x": 256, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z256_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z256_2T_grids8_ddsim_analysis.root"},
                {"x": 320, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z320_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z320_2T_grids8_ddsim_analysis.root"},
                {"x": 512, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z512_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z512_2T_grids8_ddsim_analysis.root"},
            ],
        },

        {
            "name": "bfield_scan",
            "outdir": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/occupancy_studies/detector_magnetic_field",
            "x_title": "CLD_o2_v07 solenoid B field (T)",
            "connect": True,
            "points": [
                {"x": 0.0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_ddsim_analysis.root"},
                {"x": 0.1, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0p1T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0p1T_grids8_ddsim_analysis.root"},
                {"x": 0.2, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0p2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0p2T_grids8_ddsim_analysis.root"},
                {"x": 0.3, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0p3T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0p3T_grids8_ddsim_analysis.root"},
                {"x": 0.4, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0p4T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0p4T_grids8_ddsim_analysis.root"},
                {"x": 0.5, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0p5T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0p5T_grids8_ddsim_analysis.root"},
                {"x": 0.75, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0p75T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0p75T_grids8_ddsim_analysis.root"},
                {"x": 1.0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_1T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_1T_grids8_ddsim_analysis.root"},
                {"x": 2.0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
                {"x": 3.0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_3T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_3T_grids8_ddsim_analysis.root"},
                {"x": 4.0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_4T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_4T_grids8_ddsim_analysis.root"},
            ],
        },

        {
            "name": "grids_scan_2T",
            "outdir": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/occupancy_studies/grids/2T",
            "x_title": "grids",
            "connect": True,
            "points": [
                {"x": 0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1/pairs0/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_ddsim_analysis.root"},
                {"x": 1, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_ddsim_analysis.root"},
                {"x": 2, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids2/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids2_ddsim_analysis.root"},
                {"x": 3, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids3/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids3_ddsim_analysis.root"},
                {"x": 4, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids4/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids4_ddsim_analysis.root"},
                {"x": 5, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids5/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids5_ddsim_analysis.root"},
                {"x": 6, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids6/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids6_ddsim_analysis.root"},
                {"x": 7, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids7/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids7_ddsim_analysis.root"},
                {"x": 8, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_ddsim_analysis.root"},
            ],
        },

        {
            "name": "grids_scan_2T_nbf",
            "outdir": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/occupancy_studies/grids/2T_nbf",
            "x_title": "grids",
            "connect": True,
            "points": [
                {"x": 0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_nbf/pairs0/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_nbf_ddsim_analysis.root"},
                {"x": 1, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids1_nbf_ddsim_analysis.root"},
                {"x": 2, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids2_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids2_nbf_ddsim_analysis.root"},
                {"x": 3, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids3_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids3_nbf_ddsim_analysis.root"},
                {"x": 4, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids4_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids4_nbf_ddsim_analysis.root"},
                {"x": 5, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids5_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids5_nbf_ddsim_analysis.root"},
                {"x": 6, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids6_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids6_nbf_ddsim_analysis.root"},
                {"x": 7, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids7_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids7_nbf_ddsim_analysis.root"},
                {"x": 8, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_2T_grids8_nbf_ddsim_analysis.root"},
            ],
        },

        {
            "name": "grids_scan_0T",
            "outdir": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/occupancy_studies/grids/0T",
            "x_title": "grids",
            "connect": True,
            "points": [
                {"x": 0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1/pairs0/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_ddsim_analysis.root"},
                {"x": 1, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_ddsim_analysis.root"},
                {"x": 2, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids2/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids2_ddsim_analysis.root"},
                {"x": 3, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids3/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids3_ddsim_analysis.root"},
                {"x": 4, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids4/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids4_ddsim_analysis.root"},
                {"x": 5, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids5/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids5_ddsim_analysis.root"},
                {"x": 6, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids6/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids6_ddsim_analysis.root"},
                {"x": 7, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids7/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids7_ddsim_analysis.root"},
                {"x": 8, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_ddsim_analysis.root"},
            ],
        },

        {
            "name": "grids_scan_0T_nbf",
            "outdir": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/occupancy_studies/grids/0T_nbf",
            "x_title": "grids",
            "connect": True,
            "points": [
                {"x": 0, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_nbf/pairs0/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_nbf_ddsim_analysis.root"},
                {"x": 1, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids1_nbf_ddsim_analysis.root"},
                {"x": 2, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids2_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids2_nbf_ddsim_analysis.root"},
                {"x": 3, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids3_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids3_nbf_ddsim_analysis.root"},
                {"x": 4, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids4_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids4_nbf_ddsim_analysis.root"},
                {"x": 5, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids5_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids5_nbf_ddsim_analysis.root"},
                {"x": 6, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids6_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids6_nbf_ddsim_analysis.root"},
                {"x": 7, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids7_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids7_nbf_ddsim_analysis.root"},
                {"x": 8, "path": "/home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_nbf/pairs/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_nbf_ddsim_analysis.root"},
            ],
        },
    ]

    for st in studies:
        run_study(study=st, param_defs=PARAM_DEFS)


if __name__ == "__main__":
    main()
