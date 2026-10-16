import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wave_tools import DATA_DIR, FIGURE_DIR, L2_error, L2_norm, check_run, exact_solution, load_run

DELTA_XS = [0.2, 0.1, 0.05, 0.025]
COURANT_NUMBER = 0.5
FINAL_TIME = 5.0


def measure_error(stencil_order, delta_x):
    delta_t = COURANT_NUMBER * delta_x
    csv = DATA_DIR / f"convergence_order{stencil_order}_dx{delta_x}.csv"
    check_run(delta_x, delta_t, stencil_order, FINAL_TIME, FINAL_TIME, csv)
    t, x, phi = load_run(csv)
    exact = exact_solution(t[-1], x)
    absolute_error = L2_error(phi[-1], exact, delta_x)
    relative_error = absolute_error / L2_norm(exact, delta_x)
    return absolute_error, relative_error


def study_stencil(stencil_order):
    errors = [measure_error(stencil_order, delta_x) for delta_x in DELTA_XS]
    absolute_errors = np.array([error[0] for error in errors])
    relative_errors = np.array([error[1] for error in errors])
    return absolute_errors, relative_errors


def report(stencil_order, absolute_errors, relative_errors):
    delta_xs = np.array(DELTA_XS)
    local_orders = np.diff(np.log(absolute_errors)) / np.diff(np.log(delta_xs))
    fitted_order = np.polyfit(np.log(delta_xs), np.log(absolute_errors), 1)[0]
    print(f"Stencil order {stencil_order}, Courant number {COURANT_NUMBER}, evolved to t = {FINAL_TIME}")
    print(f"{'delta_x':>8}  {'L2 error':>10}  {'relative error':>14}  {'local order':>11}")
    for index, delta_x in enumerate(DELTA_XS):
        order_text = "" if index == 0 else f"{local_orders[index - 1]:11.2f}"
        print(f"  {delta_x:8.4f}  {absolute_errors[index]:10.3e}  {relative_errors[index]:14.3e}  {order_text}")
    print(f"Fitted order of convergence: {fitted_order:.2f}")
    return fitted_order


def main():
    FIGURE_DIR.mkdir(exist_ok = True)
    delta_xs = np.array(DELTA_XS)
    fig, ax = plt.subplots(figsize = (6, 4.5))
    for stencil_order, marker in ((2, "o"), (4, "s")):
        absolute_errors, relative_errors = study_stencil(stencil_order)
        fitted_order = report(stencil_order, absolute_errors, relative_errors)
        ax.loglog(delta_xs, absolute_errors, marker + "-", label = f"Stencil order {stencil_order} (fitted slope {fitted_order:.2f})")
        ideal = absolute_errors[0] * (delta_xs / delta_xs[0]) ** stencil_order
        ax.loglog(delta_xs, ideal, "k:", lw = 0.8)
    ax.set_xlabel(r"$\Delta x$")
    ax.set_ylabel(rf"$L^2$ error at $t = {FINAL_TIME}$")
    ax.set_title(r"Convergence at fixed $\Delta t / \Delta x = 0.5$ (dotted: ideal slopes)")
    ax.grid(alpha = 0.3, which = "both")
    ax.legend()
    fig.tight_layout()
    out = FIGURE_DIR / "convergence.png"
    fig.savefig(out, dpi = 150)
    print(f"\nSaved the convergence plot to {out}")


if __name__ == "__main__":
    main()
