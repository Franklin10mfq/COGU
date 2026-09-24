"""
Subprocess worker that runs CVXPYgen with a clean heap.

Invoked by cgen.generate_solver_subprocess:

    python _cpg_worker.py <payload.pkl>

The pickled payload carries:
    package_root : path to prepend to sys.path so the package imports
    build_kwargs : arguments to build_problem (slices and numpy arrays
                   pickle natively, unlike JSON)
    solver_dir   : destination directory for the generated C
    solver       : the embedded solver to target, e.g. 'ECOS'
    extra_path   : optional directories to prepend to PATH (for cmake)

Why a separate process: DPP canonicalization needs a large contiguous
allocation.  In the parent, SymPy and the SCvx solve have already fragmented
the heap, and the allocation fails.  A fresh child has the whole block.
"""

import os
import pickle
import sys


def main():
    with open(sys.argv[1], 'rb') as fp:
        payload = pickle.load(fp)

    # CVXPYgen also tries to compile an optional Python wrapper, which needs a
    # toolchain on PATH.  If it is missing the wrapper fails but the C code --
    # the part that matters here -- has already been written.
    extra_path = payload.get('extra_path')
    if extra_path:
        os.environ['PATH'] = os.pathsep.join(list(extra_path) + [os.environ.get('PATH', '')])

    sys.path.insert(0, payload['package_root'])
    from scvx_toolbox.problem import build_problem

    # cvxpygen imports every solver interface eagerly, and the pdaqp one pulls
    # in juliacall, which tries to download a Julia toolchain on first import.
    # We never use pdaqp, so stub it out rather than require network access.
    if 'pdaqp' not in sys.modules:
        import types
        _stub = types.ModuleType('pdaqp')
        _stub.MPQP = object
        sys.modules['pdaqp'] = _stub

    from cvxpygen import cpg

    # cp.Problem is not picklable, but its inputs are: rebuild it identically.
    prob_dict = build_problem(**payload['build_kwargs'])
    solver_dir = payload['solver_dir']

    try:
        cpg.generate_code(prob_dict['problem'], code_dir=solver_dir,
                          solver=payload.get('solver', 'ECOS'))
    except Exception as exc:
        c_src = os.path.join(solver_dir, 'c', 'src', 'cpg_solve.c')
        if os.path.isfile(c_src):
            print(f"[worker] C generated; optional Python wrapper skipped ({exc})")
        else:
            raise


if __name__ == '__main__':
    main()
