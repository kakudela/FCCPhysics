import numpy as np
import matplotlib.pyplot as plt
import os

# Input data
b_field_T = np.array([0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3], dtype=float)
max_occs   = np.array([228.0, 262.6, 309.0, 307.4, 294.8, 233.1, 194.6, 163.5, 149.5], dtype=float)
avg_occs   = np.array([172.0, 215.2, 268.8, 277.3, 262.1, 203.9, 163.5, 142.6, 125.9], dtype=float)

def plot_b_field_occupancy(b_field, y_data, label, outdir="b_field_occupancy_analysis"):
    os.makedirs(outdir, exist_ok=True)
    plt.figure(figsize=(10, 7))
    plt.plot(b_field, y_data, 'o-', ms=8, lw=2, color='blue', label='Data Points')
    plt.title(f"{label} vs. Detector B-Field", fontsize=16)
    plt.xlabel("Detector B-Field (T)", fontsize=14)
    plt.ylabel(f"{label} (hits/cell)", fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.legend(fontsize=12)
    plt.xticks(b_field)
    plt.tight_layout()
    outpath = os.path.join(outdir, f"{label.lower().replace(' ', '_')}_vs_b_field.png")
    plt.savefig(outpath)
    plt.close()
    print(f"[{label}] Plot saved to {outpath}")

if __name__ == "__main__":
    plot_b_field_occupancy(b_field_T, max_occs, "Max Occupancy")
    plot_b_field_occupancy(b_field_T, avg_occs, "Average Occupancy")
