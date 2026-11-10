import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wave_tools import DATA_DIR, FIGURE_DIR, L2_norm, check_run, load_run

DELTA_X_COARSE = 0.1
DELTA_X_MEDIUM = 0.05
DELTA_X_HIGH = 0.025
COURANT_NUMBER = 0.5
FINAL_TIME = 20.0
SAVE_INTERVAL = 0.05
SNAPSHOT_TIMES = [2.5, 5.0, 7.5]


def grid_stride(delta_x_fine):
    ratio = DELTA_X_COARSE / delta_x_fine
    stride = int(round(ratio))
    if abs(ratio - stride) > 1e-9:
        raise ValueError("Error in diff between ratio and stride")
    return stride


def run_resolution(stencil_order, delta_x):
    csv = DATA_DIR / f"self_convergence_order{stencil_order}_dx{delta_x}.csv"
    check_run(delta_x, COURANT_NUMBER * delta_x, stencil_order, FINAL_TIME, SAVE_INTERVAL, csv)
    return load_run(csv)


def compute_differences(stencil_order):
    t, x, phi_coarse = run_resolution(stencil_order, DELTA_X_COARSE)
    t_medium, x_medium, phi_medium = run_resolution(stencil_order, DELTA_X_MEDIUM)
    t_high, x_high, phi_high = run_resolution(stencil_order, DELTA_X_HIGH)
    if not (np.allclose(t, t_medium) and np.allclose(t, t_high)):
        raise RuntimeError("Runs do not share the same output times")
    phi_medium = phi_medium[:, ::grid_stride(DELTA_X_MEDIUM)]
    phi_high = phi_high[:, ::grid_stride(DELTA_X_HIGH)]
    coarse_minus_medium = phi_coarse - phi_medium
    medium_minus_high = phi_medium - phi_high
    return t, x, coarse_minus_medium, medium_minus_high


def convergence_factor(coarse_minus_medium, medium_minus_high):
    numerator = np.array([L2_norm(row, DELTA_X_COARSE) for row in coarse_minus_medium])
    denominator = np.array([L2_norm(row, DELTA_X_COARSE) for row in medium_minus_high])
    with np.errstate(divide = "ignore", invalid = "ignore"):
        return numerator / denominator


def convergence_factor_q(n):
    return (DELTA_X_COARSE ** n - DELTA_X_MEDIUM ** n) / (DELTA_X_MEDIUM ** n - DELTA_X_HIGH ** n)

def report(stencil_order, t, factors):
    expected = convergence_factor_q(stencil_order)
    print(f"\nStencil order {stencil_order}, resolutions dx = {DELTA_X_COARSE}, {DELTA_X_MEDIUM}, {DELTA_X_HIGH}")
    print(f"  Expected convergence factor Q_{stencil_order} = {expected:.3f}")
    for snapshot_time in SNAPSHOT_TIMES + [FINAL_TIME]:
        index = int(np.argmin(np.abs(t - snapshot_time)))
        factor = factors[index]
        print(f"  t = {t[index]:5.2f}: factor = {factor:8.3f}")


def main():
    FIGURE_DIR.mkdir(exist_ok = True)
    fig_points, axes_points = plt.subplots(2, len(SNAPSHOT_TIMES), figsize = (4.5 * len(SNAPSHOT_TIMES), 7), sharex = True)
    fig_factor, axes_factor = plt.subplots(1, 2, figsize = (11, 4))
    for row, stencil_order in enumerate((2, 4)):
        t, x, coarse_minus_medium, medium_minus_high = compute_differences(stencil_order)
        scale = convergence_factor_q(stencil_order)
        factors = convergence_factor(coarse_minus_medium, medium_minus_high)
        report(stencil_order, t, factors)
        for column, snapshot_time in enumerate(SNAPSHOT_TIMES):
            index = int(np.argmin(np.abs(t - snapshot_time)))
            ax = axes_points[row, column]
            ax.plot(x, coarse_minus_medium[index], lw = 1.8, label = r"$\Phi_c - \Phi_m$")
            ax.plot(x, scale * medium_minus_high[index], "--", lw = 1.8, label = rf"$Q_{stencil_order} \, (\Phi_m - \Phi_h)$, $Q_{stencil_order}$ = {scale:.2f}")
            ax.set_title(f"order {stencil_order}, t = {t[index]:.2f}")
            ax.grid(alpha = 0.3)
            if row == 1:
                ax.set_xlabel("x")
            if column == 0:
                ax.set_ylabel("difference")
                ax.legend(loc = "upper right", fontsize = 8)
        ax = axes_factor[row]
        ax.plot(t[1:], factors[1:], lw = 1.5)
        ax.axhline(scale, color = "k", ls = ":", lw = 1)
        ax.set_yscale("log")
        ax.set_xlabel("t")
        ax.set_ylabel(r"$\| \Phi_c - \Phi_m \|_2 \, / \, \| \Phi_m - \Phi_h \|_2$")
        ax.set_title(f"order {stencil_order}")
        ax.grid(alpha = 0.3, which = "both")
    fig_points.suptitle("Pointwise self-convergence: coarse-minus-medium against scaled medium-minus-high")
    fig_points.tight_layout()
    out_points = FIGURE_DIR / "self_convergence_pointwise.png"
    fig_points.savefig(out_points, dpi = 150)
    fig_factor.tight_layout()
    out_factor = FIGURE_DIR / "self_convergence_factor.png"
    fig_factor.savefig(out_factor, dpi = 150)

if __name__ == "__main__":
    main()
