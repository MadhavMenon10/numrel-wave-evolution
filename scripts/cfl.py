import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wave_tools import DATA_DIR, FIGURE_DIR, load_run, check_run

DELTA_X = 0.1
FINAL_TIME = 100.0
COURANT_NUMBERS = np.round(np.arange(0.20, 2.001, 0.05), 4)
BLOWUP_THRESHOLD = 10.0
PREDICTED_LIMITS = {2: np.sqrt(2.0), 4: np.sqrt(1.5)}


def largest_field_value(stencil_order, courant_number):
    delta_t = courant_number * DELTA_X
    num_steps = int(round(FINAL_TIME / delta_t))
    t_final = num_steps * delta_t
    csv = DATA_DIR / f"cfl_order{stencil_order}.csv"
    check_run(DELTA_X, delta_t, stencil_order, t_final, t_final, csv)
    t, x, phi = load_run(csv)
    final_field = phi[-1]
    if not np.all(np.isfinite(final_field)):
        return np.inf
    return float(np.max(np.abs(final_field)))


def report_stability_limit(stencil_order, largest_values):
    is_stable = largest_values < BLOWUP_THRESHOLD
    predicted = PREDICTED_LIMITS[stencil_order]
    if np.all(is_stable):
        print(f"Order {stencil_order} stencil: stable for every Courant number tested, up to {COURANT_NUMBERS[-1]}.")
        return
    first_unstable = int(np.argmax(~is_stable))
    print(
        f"Order {stencil_order} stencil: stable up to a Courant number of "
        f"{COURANT_NUMBERS[first_unstable - 1]:.2f}, unstable from {COURANT_NUMBERS[first_unstable]:.2f} "
        f"(predicted limit: {predicted:.4f})."
    )


def main():
    FIGURE_DIR.mkdir(exist_ok = True)
    fig, ax = plt.subplots(figsize = (7, 4.5))
    for stencil_order in (2, 4):
        largest_values = np.array([largest_field_value(stencil_order, c) for c in COURANT_NUMBERS])
        report_stability_limit(stencil_order, largest_values)
        drawable = np.minimum(np.where(np.isfinite(largest_values), largest_values, 1e12), 1e12)
        ax.semilogy(COURANT_NUMBERS, drawable, "o-", ms = 3, label = f"order {stencil_order}")
        ax.axvline(PREDICTED_LIMITS[stencil_order], ls = ":", lw = 1)
    ax.axhline(BLOWUP_THRESHOLD, color = "gray", lw = 1)
    ax.set_xlabel(r"Courant number $\lambda = \Delta t / \Delta x$")
    ax.set_ylabel(r"largest $|\Phi|$ at $t \approx 100$")
    ax.set_title("Stability scan (dotted lines: predicted limits)")
    ax.grid(alpha = 0.3)
    ax.legend()
    fig.tight_layout()
    out = FIGURE_DIR / "cfl_scan.png"
    fig.savefig(out, dpi = 150)
    print(f"Saved the stability plot to {out}")


if __name__ == "__main__":
    main()
