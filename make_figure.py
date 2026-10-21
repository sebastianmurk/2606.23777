#!/usr/bin/env python3
"""Generate the standalone three-panel Euclidean-flow / KM comparison."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from figure_config import STYLE, NUMERICS
from likelihood import PauliLikelihood
from km_iteration import km_trace
from plotting import render

BASE = Path(__file__).resolve().parent


def model_for(case):
    return PauliLikelihood(case["center"], case["relative_weights_xz"],
                           case["weight_y"], case["axis_angle_degrees"])


def fingerprint(scenarios):
    content = json.dumps({"scenarios": scenarios, "numerics": NUMERICS},
                         sort_keys=True).encode()
    for name in ("likelihood.py", "km_iteration.py"):
        content += (BASE / name).read_bytes()
    return hashlib.sha256(content).hexdigest()


def compute(scenarios, cache):
    cache.mkdir(parents=True, exist_ok=True)
    records = {}
    for case in scenarios:
        key, model = case["id"], model_for(case)
        flow = model.integrate(case["initial"],
                               samples=NUMERICS["flow_samples"],
                               rtol=NUMERICS["flow_rtol"], atol=NUMERICS["flow_atol"],
                               tolerance=NUMERICS["flow_endpoint_tolerance"])
        km = km_trace(model, case["initial"],
                      precision=NUMERICS["km_decimal_digits"],
                      step_tolerance=NUMERICS["km_step_tolerance"],
                      max_iterations=NUMERICS["km_max_iterations"])
        values = model.value(flow.points)
        if np.min(np.diff(values)) < -2e-11:
            raise ArithmeticError("Euclidean likelihood decreases beyond solver tolerance")
        if len(flow.crossing_points) and not np.allclose(
                np.linalg.norm(flow.crossing_points, axis=1), 1, atol=2e-10, rtol=0):
            raise ArithmeticError("Euclidean boundary event residual is too large")
        if not np.all(km.boundary_gaps > 0):
            raise ArithmeticError("A computed KM iterate is not full rank")
        mle = model.physical_maximum()
        np.savez_compressed(cache / f"scenario_{key}.npz",
                            flow_times=flow.times, flow_points=flow.points,
                            crossing_times=flow.crossing_times,
                            crossing_points=flow.crossing_points,
                            km_points=km.points, km_boundary_gaps=km.boundary_gaps,
                            mle=mle)
        records[key] = {
            "model": model.information(),
            "initial_contour_level": float(model.value(case["initial"])),
            "euclidean": {
                "boundary_crossing_times": flow.crossing_times.tolist(),
                "boundary_crossing_points": flow.crossing_points.tolist(),
                "endpoint_error": float(np.linalg.norm(flow.points[-1]-model.center)),
                "minimum_sampled_likelihood_increment": float(np.diff(values).min()),
                "maximum_sampled_radius": float(np.linalg.norm(flow.points, axis=1).max()),
                "solver_steps": flow.solver_steps,
            },
            "km": {
                "updates": km.iterations, "classification": km.classification,
                "endpoint": km.endpoint.tolist(), "last_step": km.last_step,
                "boundary_gap": km.boundary_gap,
                "likelihood_gap_to_mle": km.likelihood_gap_to_mle,
                "minimum_likelihood_increment": km.minimum_likelihood_increment,
                "maximum_determinant_identity_error": km.purity_identity_error,
            },
        }
        print(f"({key}) Euclidean crossings: {len(flow.crossing_times)}; "
              f"KM: {km.iterations} updates, {km.classification}", flush=True)
    (cache / "diagnostics.json").write_text(json.dumps({
        "fingerprint": fingerprint(scenarios), "numerics": NUMERICS,
        "interpretation": "Solid black inside the physical disk; dashed black outside. "
                          "The purple KM line joins all computed discrete iterates; "
                          "joining segments guide the eye and are not a continuous flow.",
        "panels": records,
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true", help="Reuse the included numerical cache")
    parser.add_argument("--dpi", type=int, default=STYLE["dpi"], help="PNG DPI; PDF/SVG remain vector")
    parser.add_argument("--output-dir", type=Path, default=BASE / "outputs")
    args = parser.parse_args()
    if args.dpi < 50:
        parser.error("Require dpi >= 50")
    scenarios = json.loads((BASE / "scenarios.json").read_text())
    cache = BASE / "data"
    if not args.plot_only:
        compute(scenarios, cache)
    metadata = json.loads((cache / "diagnostics.json").read_text())
    if metadata["fingerprint"] != fingerprint(scenarios):
        parser.error("Physics or numerical settings changed; rerun without --plot-only")
    cases = []
    for case in scenarios:
        with np.load(cache / f"scenario_{case['id']}.npz", allow_pickle=False) as saved:
            data = {name: saved[name].copy() for name in saved.files}
        cases.append((case, model_for(case), data))
    render(cases, args.output_dir, args.dpi)


if __name__ == "__main__":
    main()
