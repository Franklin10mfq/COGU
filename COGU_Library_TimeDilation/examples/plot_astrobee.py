"""
Astrobee trajectory plots -- the same six figures `plot_trajectory.py` drew in
COGU_Library_Enhanced, plus the 3D view with the keep-out spheres.

With --reference it overlays a saved run so two libraries (or Python vs the
generated C) can be compared on one axis instead of by eye across two windows.

    python examples/plot_astrobee.py --save out/
    python examples/plot_astrobee.py --reference ref_enhanced.npz --save out/
    python examples/plot_astrobee.py --save out/ --html out/astrobee_3d.html

A reference file is an .npz holding `x` (nx, T+1) and `u` (nu, T).  Running
this script with --dump writes one in that format.
"""

import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
import numpy as np

from astrobee_model import OBSTACLES, astrobee_problem
from scvx_toolbox import solve_trajectory

# rows plotted per figure: (title, ylabel, state/control, indices, labels, step?)
PANELS = [
    ('quaternion', 'Quaternion sequence', 'x', [6, 7, 8, 9],
     ['q0', 'q1', 'q2', 'q3'], False),
    ('position', 'Position sequence [m]', 'x', [0, 1, 2],
     ['x', 'y', 'z'], False),
    ('velocity', 'Velocity sequence [m/s]', 'x', [3, 4, 5],
     ['vx', 'vy', 'vz'], False),
    ('angular_velocity', 'Angular velocity sequence [rad/s]', 'x', [10, 11, 12],
     ['wx', 'wy', 'wz'], False),
    ('acceleration', 'Control acceleration sequence [m/s^2]', 'u', [0, 1, 2],
     ['ax', 'ay', 'az'], True),
    ('torque', 'Control torque sequence [Nm]', 'u', [3, 4, 5],
     ['taux', 'tauy', 'tauz'], True),
]


def draw(x, u, tf, save_dir, reference=None, ref_label='reference',
         label='scvx_toolbox'):
    import matplotlib.pyplot as plt

    T = x.shape[1] - 1
    t_x = np.linspace(0.0, tf, T + 1)
    t_u = np.linspace(0.0, tf, u.shape[1])
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for name, ylabel, var, idx, labels, is_step in PANELS:
        arr = x if var == 'x' else u
        t = t_x if var == 'x' else t_u
        fig, ax = plt.subplots(figsize=(7, 4))
        plot = ax.step if is_step else ax.plot
        for n, (i, lab) in enumerate(zip(idx, labels)):
            plot(t, arr[i, :], label=lab, color=colors[n % len(colors)])
        if name == 'quaternion':
            ax.plot(t, np.linalg.norm(x[6:10, :], axis=0), '--',
                    color='0.4', label='||q||')
        if reference is not None:
            ref = reference['x'] if var == 'x' else reference['u']
            t_ref = (np.linspace(0.0, tf, ref.shape[1]) if var == 'u' else t_x)
            for n, i in enumerate(idx):
                (ax.step if is_step else ax.plot)(
                    t_ref, ref[i, :], ':', linewidth=2.2,
                    color=colors[n % len(colors)], alpha=0.75,
                    label=None if n else ref_label)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(ylabel)
        ax.grid(True)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        if save_dir:
            path = os.path.join(save_dir, f'astrobee_{name}.png')
            fig.savefig(path, dpi=130)
            print(f"  wrote {path}")

    # 3D trajectory with the keep-out spheres
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x[0], x[1], x[2], color=colors[0], linewidth=2, label=label)
    if reference is not None:
        r = reference['x']
        ax.plot(r[0], r[1], r[2], ':', color='0.3', linewidth=2.2,
                label=ref_label)
    ug, vg = np.mgrid[0:2 * np.pi:24j, 0:np.pi:12j]
    for centre, radius in OBSTACLES:
        ax.plot_surface(radius * np.cos(ug) * np.sin(vg) + centre[0],
                        radius * np.sin(ug) * np.sin(vg) + centre[1],
                        radius * np.cos(vg) + centre[2],
                        color='tab:red', alpha=0.22, linewidth=0)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    # Equal span on all three axes, or the keep-out spheres are drawn as
    # ellipsoids and the clearance is impossible to read off the picture.
    pad = 0.5
    span = max(x[i].max() - x[i].min() for i in range(3)) / 2 + pad
    for setter, i in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
        mid = (x[i].max() + x[i].min()) / 2
        setter(mid - span, mid + span)
    ax.set_box_aspect((1, 1, 1))
    ax.legend(fontsize=8)
    ax.set_title('Astrobee guidance (SCvx, ECOS)')
    fig.tight_layout()
    if save_dir:
        path = os.path.join(save_dir, 'astrobee_3d.png')
        fig.savefig(path, dpi=130)
        print(f"  wrote {path}")
    return plt


def clearance(x):
    """Minimum distance to each keep-out sphere surface, over the trajectory."""
    out = []
    for centre, radius in OBSTACLES:
        d = np.linalg.norm(x[0:3, :] - centre.reshape(3, 1), axis=0)
        out.append(d.min() - radius)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--time-mode', default='ZOH')
    ap.add_argument('--T', type=int, default=30)
    ap.add_argument('--tf', type=float, default=200.0)
    ap.add_argument('--legacy', action='store_true',
                    help='reproduce COGU_Library_Enhanced exactly')
    ap.add_argument('--reference', help='.npz with x and u to overlay')
    ap.add_argument('--ref-label', default='reference')
    ap.add_argument('--save', help='directory to write PNGs into')
    ap.add_argument('--dump', help='.npz to write this run into')
    ap.add_argument('--show', action='store_true')
    args = ap.parse_args()

    if not args.show:
        matplotlib.use('Agg')

    kwargs = astrobee_problem(T=args.T, tf=args.tf, time_mode=args.time_mode,
                              legacy_quaternion=args.legacy)
    result = solve_trajectory(size_N=20, max_iter=25,
                              pin_aware_trust_region=not args.legacy,
                              defect_penalty_space=('unscaled' if args.legacy
                                                    else 'scaled'),
                              verbose=True, **kwargs)

    x = result['x_opt_unscaled']
    u = result['u_opt_unscaled']

    print(f"\nstatus {result['status']}, {result['iterations']} iterations, "
          f"tf = {result['tf']:.4f} s")
    for i, c in enumerate(clearance(x), start=1):
        print(f"  clearance to obstacle {i}: {c:+.4f} m")

    if args.dump:
        np.savez(args.dump, x=x, u=u)
        print(f"  wrote {args.dump}")

    ref = np.load(args.reference) if args.reference else None
    if args.save:
        os.makedirs(args.save, exist_ok=True)
    plt = draw(x, u, result['tf'], args.save, ref, args.ref_label)
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
