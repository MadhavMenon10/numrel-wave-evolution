import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wave_tools import DATA_DIR, FIGURE_DIR, exact_solution, load_run, check_run


TIMES = list(range(21))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stencil_order", type = int, default = 2, choices = (2, 4))
    parser.add_argument("--delta_x", type = float, default = 0.1)
    parser.add_argument("--delta_t", type = float, default = 0.05)
    parser.add_argument("--t_final", type = float, default = 20.0)
    args = parser.parse_args()
    csv = DATA_DIR / f"anim_order{args.stencil_order}.csv"
    check_run(args.delta_x, args.delta_t, args.stencil_order, args.t_final, args.delta_t, csv)
    t, x, phi = load_run(csv)
    FIGURE_DIR.mkdir(exist_ok = True)
    fig, axes = plt.subplots(3, 7, figsize = (22, 8), sharex = True, sharey = True)
    for ax, t0 in zip(axes.ravel(), TIMES):
        r = int(np.argmin(np.abs(t - t0)))
        ax.plot(x, exact_solution(t[r], x), "k--", lw = 1, label = "exact")
        ax.plot(x, phi[r], lw = 1.5, label = "numerical")
        ax.set_title(f"t = {t[r]:.2f}", fontsize = 9)
        ax.grid(alpha = 0.3)
    for ax in axes[-1]:
        ax.set_xlabel("x")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\Phi$")
    axes[0, 0].legend(loc = "upper right", fontsize = 8)
    fig.suptitle(
        f"Stencil order {args.stencil_order} "
        f"(delta_x = {args.delta_x}, delta_t = {args.delta_t})"
    )
    fig.tight_layout()
    out = FIGURE_DIR / f"snapshots_order{args.stencil_order}.png"
    fig.savefig(out, dpi = 150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
