import pathlib, sys, os, math

IN_DIR   = pathlib.Path("/ceph/submit/data/group/fcc/ee/detector/guineapig/FCCee_Z_4IP_04may23_FCCee_Z128").resolve()
OUT_BASE = pathlib.Path("/ceph/submit/data/user/k/kudela/beam_backgrounds/altered_gp_pairs").resolve()

HX = 1.32555e6     # nm
HY = 0.00186e6
HZ = 25.4e6
TOL = 1e-6         # nm

def project_to_box(x, y, z, px, py, pz):
    if abs(x) <= HX+TOL and abs(y) <= HY+TOL and abs(z) <= HZ+TOL:
        return x, y, z, False, 0.0, 0.0

    dx, dy, dz = -px, -py, -pz
    t_enter, t_exit = -float('inf'), float('inf')

    for coord, d, half in ((x, dx, HX), (y, dy, HY), (z, dz, HZ)):
        if d == 0.0:            # parallel
            continue
        t1 = (+half - coord) / d
        t2 = (-half - coord) / d
        if t1 > t2: t1, t2 = t2, t1
        t_enter = max(t_enter, t1)
        t_exit  = min(t_exit,  t2)

    miss = t_enter > t_exit
    if not miss and t_exit >= 0:
        x += dx * t_enter
        y += dy * t_enter
        z += dz * t_enter
    else:                       # clamp fallback
        x = max(min(x, HX), -HX)
        y = max(min(y, HY), -HY)
        z = max(min(z, HZ), -HZ)
    return x, y, z, miss, t_enter, t_exit

def main():
    if not IN_DIR.is_dir(): sys.exit(f"Input directory not found: {IN_DIR}")
    files = sorted(IN_DIR.rglob("*.pairs"))
    if not files: sys.exit("No *.pairs files under input dir")

    OUT_ROOT = OUT_BASE / IN_DIR.name

    total=outside_orig=outside_new=0
    int_miss=b_tot=bx=by=bz=bxy=bxz=byz=bxyz=0
    kept=inside_kept=0

    int_raw = []; int_fmt = []
    back_raw= []; back_fmt= []

    for src in files:
        dst = OUT_ROOT / src.relative_to(IN_DIR)
        dst.parent.mkdir(parents=True, exist_ok=True)

        with src.open() as fin, dst.open("w") as fout:
            for ln in fin:
                if not ln.strip() or ln.lstrip().startswith("#"):
                    fout.write(ln); continue
                cols = ln.split()
                if len(cols) < 7:  fout.write(ln); continue

                E, px, py, pz = map(float, cols[:4])
                x, y, z       = map(float, cols[4:7])

                total += 1
                if abs(x)>HX or abs(y)>HY or abs(z)>HZ: outside_orig += 1

                bx_f = (x>HX and px<0) or (x<-HX and px>0)
                by_f = (y>HY and py<0) or (y<-HY and py>0)
                bz_f = (z>HZ and pz<0) or (z<-HZ and pz>0)
                back_flag = bx_f or by_f or bz_f
                if back_flag:
                    b_tot += 1
                    if   bx_f and by_f and bz_f: bxyz += 1
                    elif bx_f and by_f:          bxy  += 1
                    elif bx_f and bz_f:          bxz  += 1
                    elif by_f and bz_f:          byz  += 1
                    elif bx_f:                   bx   += 1
                    elif by_f:                   by   += 1
                    else:                        bz   += 1

                x_new,y_new,z_new,miss_flag,t_en,t_ex = project_to_box(x,y,z,px,py,pz)

                if miss_flag: int_miss += 1
                if abs(x_new)>HX+TOL or abs(y_new)>HY+TOL or abs(z_new)>HZ+TOL:
                    outside_new += 1

                # collect example lines
                if miss_flag and len(int_raw)<15:
                    int_raw.append(ln.rstrip())
                    int_fmt.append(
                        f"E={E:.6f}  x={x:.6f}  y={y:.6f}  z={z:.6f}  "
                        f"px={px:.6f}  py={py:.6f}  pz={pz:.6f}  "
                        f"t_enter={t_en:.3f}  t_exit={t_ex:.3f}")
                if back_flag and len(back_raw)<15:
                    back_raw.append(ln.rstrip())
                    back_fmt.append(
                        f"E={E:.6f}  x={x:.6f}  y={y:.6f}  z={z:.6f}  "
                        f"px={px:.6f}  py={py:.6f}  pz={pz:.6f}  "
                        f"t_enter={t_en:.3f}  t_exit={t_ex:.3f}")

                # write only tracks that are neither pathological case
                if miss_flag or back_flag:
                    continue

                kept += 1
                if abs(x_new)<=HX+TOL and abs(y_new)<=HY+TOL and abs(z_new)<=HZ+TOL:
                    inside_kept += 1

                cols[4:7] = [f"{x_new:.6f}",f"{y_new:.6f}",f"{z_new:.6f}"]
                fout.write(" ".join(cols)+"\n")

    print(f"\nOriginal outside: {outside_orig}/{total} ({outside_orig/total:.2%})")
    print(f"Altered  outside: {outside_new}/{total} ({outside_new/total:.2%})")
    print(f"Interval-miss    : {int_miss}/{total} ({int_miss/total:.2%})")
    print(f"Backward-miss    : {b_tot}/{total} ({b_tot/total:.2%})")
    if b_tot:
        f = lambda n: f"{n} ({n/b_tot:.2%})"
        print("   only x :", f(bx))
        print("   only y :", f(by))
        print("   only z :", f(bz))
        print("   x & y  :", f(bxy))
        print("   x & z  :", f(bxz))
        print("   y & z  :", f(byz))
        print("   x y z  :", f(bxyz))
    if kept:
        print(f"\nPost-filter inside: {inside_kept}/{kept} ({inside_kept/kept:.2%})")

    # ---- print samples -------------------------------------------------
    def dump(title, raw, fmt):
        print(f"\n======= {title} (showing {len(raw)} examples) =======")
        for r,f in zip(raw,fmt):
            print(r)
        print("--- formatted ---")
        for r,f in zip(raw,fmt):
            print(f)

    dump("Interval-miss originals", int_raw, int_fmt)
    dump("Backward-momentum originals", back_raw, back_fmt)

    print("\nAltered files saved under:", OUT_ROOT)

if __name__=="__main__":
    main()
