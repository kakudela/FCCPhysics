import os
import math
import shutil
from array import array
import ROOT
import re
import threading


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def ensure_index_php(dst_dir, index_php_source):
    if not index_php_source:
        return
    ensure_dir(dst_dir)
    dst = os.path.join(dst_dir, "index.php")
    if os.path.exists(dst):
        return
    try:
        shutil.copy2(index_php_source, dst)
    except Exception:
        pass


def ensure_public_dir(public_ddsim_base, kind, gp_tag, stream, index_php_source):
    p0 = os.path.join(public_ddsim_base, kind)
    p1 = os.path.join(p0, gp_tag)
    p2 = os.path.join(p1, stream)

    ensure_dir(p0)
    ensure_index_php(p0, index_php_source)

    ensure_dir(p1)
    ensure_index_php(p1, index_php_source)

    ensure_dir(p2)
    ensure_index_php(p2, index_php_source)

    return p2

def apply_buffer(a, b, buf, hard_min, hard_max):
    if buf <= 1.0:
        return a, b
    w = b - a
    extra = 0.5 * (buf - 1.0) * w
    lo = a - extra
    hi = b + extra
    if hard_min is not None:
        lo = max(lo, hard_min)
    if hard_max is not None:
        hi = min(hi, hard_max)
    return lo, hi


def _buffer_axes_1d(h, xbuffer):
    if xbuffer <= 1.0:
        return
    xax = h.GetXaxis()
    lo, hi = xax.GetXmin(), xax.GetXmax()
    lo2, hi2 = apply_buffer(lo, hi, xbuffer, lo, hi)
    xax.SetLimits(lo2, hi2)


def _buffer_axes_2d(h, xbuffer, ybuffer):
    xax = h.GetXaxis()
    yax = h.GetYaxis()

    xlo, xhi = xax.GetXmin(), xax.GetXmax()
    ylo, yhi = yax.GetXmin(), yax.GetXmax()

    if xbuffer > 1.0:
        xlo2, xhi2 = apply_buffer(xlo, xhi, xbuffer, xlo, xhi)
        xax.SetLimits(xlo2, xhi2)
    if ybuffer > 1.0:
        ylo2, yhi2 = apply_buffer(ylo, yhi, ybuffer, ylo, yhi)
        yax.SetLimits(ylo2, yhi2)


def save_hist1d(h, out_png, title, xlabel, ylabel, scale=1.0, xbuffer=1.0, logy=False):
    if h is None:
        return
    hc = h.Clone()
    hc.SetDirectory(0)
    if scale != 1.0:
        hc.Scale(scale)

    c = ROOT.TCanvas("c1", "", 800, 600)
    hc.SetTitle(title)
    hc.GetXaxis().SetTitle(xlabel)
    hc.GetYaxis().SetTitle(ylabel)

    _buffer_axes_1d(hc, xbuffer)

    # Always start at 0 for linear plots
    if logy:
        c.SetLogy(True)
        hc.SetMinimum(1e-12)   # log scale cannot start at 0
    else:
        hc.SetMinimum(0.0)

    hc.SetLineWidth(2)
    hc.Draw("HIST")
    c.SaveAs(out_png)
    c.Close()


def save_hist2d(h, out_png, title, xlabel, ylabel, scale=1.0, xbuffer=1.0, ybuffer=1.0):
    if h is None:
        return
    hc = h.Clone()
    hc.SetDirectory(0)
    if scale != 1.0:
        hc.Scale(scale)

    c = ROOT.TCanvas("c2", "", 900, 750)
    c.SetRightMargin(0.15)

    hc.SetTitle(title)
    hc.GetXaxis().SetTitle(xlabel)
    hc.GetYaxis().SetTitle(ylabel)

    # Remove 2D stats box
    old_optstat = ROOT.gStyle.GetOptStat()
    ROOT.gStyle.SetOptStat(0)
    hc.SetStats(0)

    _buffer_axes_2d(hc, xbuffer, ybuffer)

    hc.Draw("COLZ")
    c.SaveAs(out_png)
    c.Close()

    ROOT.gStyle.SetOptStat(old_optstat)

def save_hist2d_with_box(h, out_png, title, xlabel, ylabel, lxmin, lymin, lxmax, lymax, scale=1.0):
    if h is None:
        return

    hc = h.Clone()
    hc.SetDirectory(0)
    if scale != 1.0:
        hc.Scale(scale)

    # Axis ranges (in user coords)
    xax = hc.GetXaxis()
    yax = hc.GetYaxis()
    xmin = float(xax.GetXmin())
    xmax = float(xax.GetXmax())
    ymin = float(yax.GetXmin())
    ymax = float(yax.GetXmax())

    def _finite(x):
        try:
            return math.isfinite(float(x))
        except Exception:
            return False

    def _clamp(v, lo, hi):
        return max(lo, min(hi, v))

    # Convert + clamp requested box coords to visible ranges
    # If non-finite, fall back to axis edge.
    x1 = float(lxmin) if _finite(lxmin) else xmin
    y1 = float(lymin) if _finite(lymin) else ymin
    x2 = float(lxmax) if _finite(lxmax) else xmax
    y2 = float(lymax) if _finite(lymax) else ymax

    # Clamp
    x1 = _clamp(x1, xmin, xmax)
    x2 = _clamp(x2, xmin, xmax)
    y1 = _clamp(y1, ymin, ymax)
    y2 = _clamp(y2, ymin, ymax)

    # Ensure proper ordering (TBox expects x1<x2, y1<y2)
    if x2 < x1:
        x1, x2 = x2, x1
    if y2 < y1:
        y1, y2 = y2, y1

    c = ROOT.TCanvas("c3", "", 900, 750)
    c.SetRightMargin(0.15)

    hc.SetTitle(title)
    hc.GetXaxis().SetTitle(xlabel)
    hc.GetYaxis().SetTitle(ylabel)

    old_optstat = ROOT.gStyle.GetOptStat()
    ROOT.gStyle.SetOptStat(0)
    hc.SetStats(0)

    hc.Draw("COLZ")

    box = ROOT.TBox(x1, y1, x2, y2)
    box.SetFillStyle(0)
    box.SetLineColor(ROOT.kRed)
    box.SetLineWidth(3)
    box.Draw("same")

    c.SaveAs(out_png)
    c.Close()

    ROOT.gStyle.SetOptStat(old_optstat)

def hist_mean_and_p95(h):
    if h is None:
        return (0.0, 0.0)
    mean = h.GetMean()
    probs = array("d", [0.95])
    q = array("d", [0.0])
    h.GetQuantiles(1, q, probs)
    return (float(mean), float(q[0]))


def sensor_area_mm2(R_mm, z_half_mm):
    return (2.0 * math.pi * R_mm) * (2.0 * z_half_mm)


def pixel_area_mm2(pX_um, pY_um):
    return (pX_um * pY_um) * 1e-6


def write_text(path, lines):
    ensure_dir(os.path.dirname(path))
    with open(path, "w") as f:
        for line in lines:
            f.write(str(line) + "\n")


def compute_occ_from_hist(h_tot, sensorA, pixA, cls, safety_factor):
    avg_hits, p95_hits = hist_mean_and_p95(h_tot)
    O_avg = (avg_hits / sensorA) * pixA * cls * safety_factor
    O_p95 = (p95_hits / sensorA) * pixA * cls * safety_factor
    return avg_hits, p95_hits, O_avg, O_p95


def compute_global_occ_from_hits(hits, sensorA, pixA, cls, safety_factor):
    return (hits / sensorA) * pixA * cls * safety_factor

_badfile_re = re.compile(r'opening file "([^"]+)"')
_cache_lock = threading.Lock()

def _extract_bad_file_from_exc(exc):
    m = _badfile_re.search(str(exc))
    return m.group(1) if m else None

def filter_good_files(files, treename="events", max_report=10):
    good = []
    n_zombie = 0
    n_notree = 0
    n_exc = 0

    for f in files:
        tf = None
        try:
            tf = ROOT.TFile.Open(f)
            if not tf or tf.IsZombie():
                n_zombie += 1
                if n_zombie <= max_report:
                    print(f"[prune] zombie/truncated: {f}")
                continue

            t = tf.Get(treename)
            if not t:
                n_notree += 1
                if n_notree <= max_report:
                    print(f"[prune] missing tree '{treename}': {f}")
                continue

            good.append(f)

        except Exception as e:
            n_exc += 1
            if n_exc <= max_report:
                print(f"[prune] exception: {f} -> {e}")
            continue

        finally:
            try:
                if tf:
                    tf.Close()
            except Exception:
                pass

    print(f"[prune] kept={len(good)} zombie={n_zombie} notree={n_notree} exc={n_exc}")
    return good

def weed_files_on_open(files, treename="events", max_loops=1):
    files = list(files)
    loops = 0
    while True:
        if not files:
            return files
        loops += 1
        if loops > max_loops:
            return files
        try:
            probe = ROOT.RDataFrame(treename, files)
            _ = int(probe.Count().GetValue())
            return files
        except Exception as exc:
            bad = _extract_bad_file_from_exc(exc)
            if bad and bad in files:
                files = [x for x in files if x != bad]
                continue
            files = files[:-1]
            continue
