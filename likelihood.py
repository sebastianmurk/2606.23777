"""Genuine informationally complete qubit likelihoods and Euclidean flows.

The drawing plane is r=(x,z), corresponding to rho=(I+x X+z Z)/2.
Two rotated orthogonal Pauli axes lie in this plane. A third Y setting,
with balanced outcomes, makes the full three-dimensional model IC.

The unconstrained mathematical extension is valid whenever BOTH rotated
coordinates are in (-1,1); positivity of rho additionally requires |r|<=1.
The continuous flows below are explicitly EUCLIDEAN, not Petz/KM updates.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize, minimize_scalar


@dataclass
class FlowResult:
    """Sampled curve plus the ODE's dense output and true disk crossings."""

    times: np.ndarray
    points: np.ndarray
    crossing_times: np.ndarray
    crossing_points: np.ndarray
    tmax: float
    evaluate: Callable
    solver_steps: int


class PauliLikelihood:
    """Weighted binary likelihood for three orthogonal Pauli measurements.

    ``weights_xz`` are the RELATIVE weights of the two in-plane settings;
    they are multiplied by ``1-y_weight``. ``center`` is the unrestricted
    likelihood maximum in drawing coordinates, which may lie outside the
    physical disk, but its rotated coordinates must remain in (-1,1).

    The omitted setting-probability constants never affect contours,
    gradients, likelihood comparisons, or optimizers.
    """

    def __init__(self, center, weights_xz=(0.5, 0.5), y_weight=0.1,
                 angle_degrees=0.0):
        self.center = np.asarray(center, dtype=float)
        self.weights_xz = np.asarray(weights_xz, dtype=float)
        self.y_weight = float(y_weight)
        self.angle_degrees = float(angle_degrees)
        if self.center.shape != (2,) or self.weights_xz.shape != (2,):
            raise ValueError("center and weights_xz must each have two entries")
        if (not np.all(np.isfinite(self.center)) or
                not np.all(self.weights_xz > 0) or
                not np.isclose(self.weights_xz.sum(), 1.0) or
                not 0 < self.y_weight < 1):
            raise ValueError("Use positive normalized weights_xz and 0<y_weight<1")
        theta = np.deg2rad(self.angle_degrees)
        self.rotation = np.array([[np.cos(theta), -np.sin(theta)],
                                  [np.sin(theta), np.cos(theta)]])
        self.means = self.to_measurement(self.center)
        if np.any(np.abs(self.means) >= 1):
            raise ValueError("center predicts invalid outcome probabilities")
        self.weights = (1 - self.y_weight) * self.weights_xz
        self._physical_maximum = None
        self._boundary_cache = {}

    def to_measurement(self, points):
        """Drawing (x,z) -> coordinates along the two rotated Pauli axes."""
        return np.asarray(points, dtype=float) @ self.rotation

    def from_measurement(self, points):
        """Rotated Pauli coordinates -> drawing (x,z)."""
        return np.asarray(points, dtype=float) @ self.rotation.T

    def value(self, points):
        """L for an array of shape (...,2); NaN outside its extension domain."""
        u = self.to_measurement(points)
        valid = np.all(np.abs(u) < 1, axis=-1)
        with np.errstate(invalid="ignore", divide="ignore"):
            terms = ((1 + self.means) / 2 * np.log((1 + u) / 2) +
                     (1 - self.means) / 2 * np.log((1 - u) / 2))
            out = terms @ self.weights - self.y_weight * np.log(2)
        return np.where(valid, out, np.nan)

    def gradient(self, points):
        """Ordinary Euclidean gradient, in drawing coordinates."""
        u = self.to_measurement(points)
        if np.any(np.abs(u) >= 1):
            raise ValueError("gradient requested outside the likelihood domain")
        gu = self.weights * (self.means - u) / (1 - u * u)
        return self.from_measurement(gu)

    def physical_maximum(self):
        """Unique physical MLE in this IC model, returned as (x,z).

        The balanced Y setting forces the optimum to y=0. For an external
        center, solve the positive Lagrange-multiplier KKT equations using
        nested bracketed scalar roots, rather than projecting the center.
        """
        if self._physical_maximum is not None:
            return self._physical_maximum.copy()
        if np.linalg.norm(self.means) <= 1:
            result = self.center.copy()
        else:
            def at_multiplier(lam):
                u = np.zeros(2)
                for j, (mj, wj) in enumerate(zip(self.means, self.weights)):
                    if mj == 0:
                        continue
                    sign, b = np.sign(mj), abs(mj)
                    uj = brentq(lambda v: wj * (b-v)/(1-v*v)-lam*v,
                                0.0, b, xtol=1e-14, rtol=1e-14)
                    u[j] = sign * uj
                return u

            high = 1.0
            while np.linalg.norm(at_multiplier(high)) > 1:
                high *= 2
            lam = brentq(lambda q: np.dot(at_multiplier(q), at_multiplier(q))-1,
                         0.0, high, xtol=1e-13, rtol=1e-13)
            result = self.from_measurement(at_multiplier(lam))
        self._physical_maximum = result.copy()
        return result

    def integrate(self, initial, tolerance=1e-10, samples=2400,
                  rtol=2e-10, atol=2e-12):
        """Integrate unconstrained Euclidean ascent with exact event location.

        In measurement coordinates, each coordinate obeys
        u'_j = w_j (m_j-u_j)/(1-u_j^2), remains between u_j(0) and m_j,
        and approaches m_j. Radau handles the widely separated time scales
        in panel (c) efficiently. No state projection or boundary stopping
        is applied: these are the mathematical unconstrained paths.
        """
        initial = np.asarray(initial, dtype=float)
        if initial.shape != (2,) or np.linalg.norm(initial) >= 1:
            raise ValueError("The initial state must be strictly inside the disk")
        u0 = self.to_measurement(initial)
        distance = np.linalg.norm(u0-self.means)
        if tolerance <= 0 or samples < 10:
            raise ValueError("Use positive tolerance and samples>=10")
        if distance <= tolerance:
            def constant(ts):
                ts = np.asarray(ts)
                return np.broadcast_to(initial, ts.shape + (2,)).copy()
            return FlowResult(np.array([0.]), initial[None, :], np.array([]),
                              np.empty((0, 2)), 0., constant, 0)

        # |u(t)-m| <= |u(0)-m| exp(-min(w)*t) gives a safe finite horizon.
        horizon = float(np.log(distance/tolerance) / np.min(self.weights) * 1.02)

        def rhs(t, u):
            return self.weights * (self.means-u) / (1-u*u)

        def jac(t, u):
            return np.diag(self.weights * (-1-u*u+2*self.means*u)/(1-u*u)**2)

        def disk_event(t, u):
            return np.dot(u, u)-1

        def converged(t, u):
            return np.linalg.norm(u-self.means)-tolerance

        converged.terminal = True
        converged.direction = -1
        sol = solve_ivp(rhs, (0., horizon), u0, method="Radau", jac=jac,
                        rtol=rtol, atol=atol, dense_output=True,
                        events=(disk_event, converged))
        if not sol.success:
            raise RuntimeError("Gradient-flow integration failed: " + sol.message)

        def evaluate(ts):
            ts = np.asarray(ts)
            if np.any(ts < 0) or np.any(ts > sol.t[-1]*(1+1e-13)):
                raise ValueError("Requested times outside the computed flow")
            return self.from_measurement(np.moveaxis(sol.sol(ts), 0, -1))

        tmax = float(sol.t[-1])
        # Log time resolves both the quick initial bend and the slow tail;
        # event and solver nodes are also retained, so no crossing is missed.
        times = np.unique(np.concatenate((
            [0.], np.geomspace(max(1e-10, tmax*1e-10), tmax, samples),
            sol.t, sol.t_events[0])))
        points = evaluate(times)
        if np.any(np.abs(self.to_measurement(points)) >= 1):
            raise RuntimeError("The integrated path left the probability domain")
        crossing_times = np.asarray(sol.t_events[0], dtype=float)
        crossing_points = evaluate(crossing_times) if crossing_times.size else np.empty((0, 2))
        return FlowResult(times, points, crossing_times, crossing_points,
                          tmax, evaluate, len(sol.t)-1)

    def boundary_maximum(self, full_bloch=True, starts=90):
        """Numerical maximum on the physical boundary, with a KKT residual.

        Full-sphere multistart is needed to test an upper level set in the
        full Bloch ball: checking only the drawing circle is insufficient.
        This is a numerical diagnostic, not an interval-certified theorem.
        ``point`` has drawing coordinates (x,z,y) if full_bloch=True.
        """
        key = (bool(full_bloch), int(starts))
        if key in self._boundary_cache:
            return self._boundary_cache[key].copy()
        if not full_bloch:
            theta = np.linspace(0., 2*np.pi, max(360, starts*8), endpoint=False)
            points = np.stack((np.cos(theta), np.sin(theta)), axis=-1)
            values = self.value(points)
            j = np.nanargmax(values)
            step = 2*np.pi/len(theta)
            opt = minimize_scalar(lambda t: -float(self.value([np.cos(t), np.sin(t)])),
                                  bounds=(theta[j]-step, theta[j]+step),
                                  method="bounded", options={"xatol": 1e-13})
            point = np.array([np.cos(opt.x), np.sin(opt.x)])
            g = self.gradient(point)
            result = {"point": point, "value": float(self.value(point)),
                      "residual": float(np.linalg.norm(g-np.dot(g, point)*point))}
        else:
            weights = np.r_[self.weights, self.y_weight]
            means = np.r_[self.means, 0.]

            def value(u):
                return np.sum(weights*((1+means)/2*np.log((1+u)/2) +
                                       (1-means)/2*np.log((1-u)/2)))

            def grad(u):
                return weights*(means-u)/(1-u*u)

            best = None
            for j in range(starts):
                y = 1-2*(j+0.5)/starts
                phi = j*np.pi*(3-np.sqrt(5))
                u0 = np.array([np.sqrt(1-y*y)*np.cos(phi),
                               np.sqrt(1-y*y)*np.sin(phi), y])
                opt = minimize(lambda u: -value(u), u0, jac=lambda u: -grad(u),
                               constraints={"type": "eq", "fun": lambda u: u@u-1,
                                            "jac": lambda u: 2*u},
                               bounds=[(-1+1e-12, 1-1e-12)]*3, method="SLSQP",
                               options={"ftol": 2e-13, "maxiter": 400})
                if (opt.success and abs(opt.x@opt.x-1) < 1e-8 and
                        (best is None or opt.fun < best.fun)):
                    best = opt
            if best is None:
                raise RuntimeError("Full-sphere likelihood maximization failed")
            u = best.x
            g = grad(u)
            point = np.r_[self.from_measurement(u[:2]), u[2]]
            result = {"point": point, "value": float(-best.fun),
                      "residual": float(np.linalg.norm(g-(g@u)*u)),
                      "starts": starts}
        self._boundary_cache[key] = result
        return result.copy()

    def information(self):
        """JSON-friendly model specification, including actual POVM weights."""
        return {"center": self.center.tolist(),
                "measurement_means_uvy": [*self.means.tolist(), 0.],
                "measurement_weights_uvy": [*self.weights.tolist(), self.y_weight],
                "angle_degrees": self.angle_degrees,
                "physical_maximum": self.physical_maximum().tolist()}
