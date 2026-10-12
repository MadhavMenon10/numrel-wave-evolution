import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from wave_tools import DATA_DIR, FIGURE_DIR, exact_solution, load_run, check_run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stencil_order", type = int, default = 2, choices = (2, 4))
    parser.add_argument("--delta_x", type = float, default = 0.1)
    parser.add_argument("--delta_t", type = float, default = 0.05)
    parser.add_argument("--t_final", type = float, default = 20.0)
    parser.add_argument("--stride", type = int, default = 2)
    parser.add_argument("--fps", type = int, default = 20)
    args = parser.parse_args()
    csv = DATA_DIR / f"anim_order{args.stencil_order}.csv"
    check_run(args.delta_x, args.delta_t, args.stencil_order, args.t_final, args.delta_t, csv)
    t, x, phi = load_run(csv)

    frames = range(0, len(t), args.stride)
    fig, ax = plt.subplots(figsize = (7, 4))
    (exact_line,) = ax.plot(x, exact_solution(t[0], x), "k--", lw = 1, label = "exact")
    (line,) = ax.plot(x, phi[0], lw = 2, label = "numerical")
    ax.set_xlim(x[0], x[-1] + (x[1] - x[0]))
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\Phi$")
    ax.grid(alpha = 0.3)
    ax.legend(loc = "upper right")
    title = ax.set_title("")

    def update(r):
        line.set_ydata(phi[r])
        exact_line.set_ydata(exact_solution(t[r], x))
        title.set_text(f"order {args.stencil_order} stencil, t = {t[r]:.2f}")
        return line, exact_line, title

    anim = FuncAnimation(fig, update, frames = frames, blit = False)
    FIGURE_DIR.mkdir(exist_ok = True)
    out = FIGURE_DIR / f"wave_order{args.stencil_order}.gif"
    anim.save(out, writer = PillowWriter(fps = args.fps))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
