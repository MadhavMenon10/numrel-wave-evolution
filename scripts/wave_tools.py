import subprocess
from pathlib import Path
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
EXECUTABLE_DIR = REPO_ROOT / "build" / "wave_solver"
DATA_DIR = REPO_ROOT / "data"
FIGURE_DIR = REPO_ROOT / "figures"

L = 20.0

def load_run(path):
    data = np.loadtxt(path, delimiter=",", ndmin=2)
    t = data[:, 0]
    phi = data[:, 1:]
    n = phi.shape[1]
    x = np.arange(n) * (L / n)
    return t, x, phi

def run_solver(delta_x, delta_t, stencil_order, t_final, save_interval, out_path):
    DATA_DIR.mkdir(exist_ok = True)
    cmd = [str(EXECUTABLE_DIR), str(delta_x), str(delta_t), str(stencil_order), str(t_final), str(save_interval), str(out_path)]
    return subprocess.run(cmd, capture_output=True, text=True)

def check_run(delta_x, delta_t, stencil_order, t_final, save_interval, out_path):
    result = run_solver(delta_x, delta_t, stencil_order, t_final, save_interval, out_path)
    if result.returncode != 0:
        raise RuntimeError(f"Solver failed")
    return out_path

def exact_solution(t, x):
    def f(s):
        return np.exp(-((s - 10.0) ** 2))
    total = np.zeros_like(x, dtype=float)
    for i in range(-10, 11):
       total += f(x - t + L * i) + f(x + t + L * i)
    return 0.5 * total


# To see if solver converges
def L2_norm(phi, delta_x):
    return float(np.sqrt(delta_x * np.sum(phi ** 2)))

def L2_error(phi_approx, phi_exact, delta_x):
    return float(np.sqrt(delta_x * np.sum((phi_approx - phi_exact) ** 2)))

