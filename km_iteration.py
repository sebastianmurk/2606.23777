"""Physical KM iteration for the same Pauli measurement models.

This module deliberately does not extend the update outside the Bloch ball.
All iteration arithmetic uses Decimal; float conversion occurs only for
graphics and summary diagnostics. No clipping, projection, or ridge is used.
"""

from dataclasses import dataclass
from decimal import Decimal, localcontext

import numpy as np


_PI = Decimal(
    "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481"
)


@dataclass
class KMTrace:
    points: np.ndarray
    iterations: int
    endpoint: np.ndarray
    last_step: float
    boundary_gap: float
    likelihood_change: float
    likelihood_gap_to_mle: float
    classification: str
    purity_identity_error: float
    minimum_likelihood_increment: float
    boundary_gaps: np.ndarray
    precision: int


def _number(value):
    return Decimal(str(value))


def _rotation(angle_degrees, precision):
    """Decimal sine/cosine, so the high-precision map is truly orthogonal."""
    theta = _number(angle_degrees) * _PI / Decimal(180)
    theta = (theta + _PI) % (2*_PI) - _PI
    sine, cosine = theta, Decimal(1)
    s_term, c_term = theta, Decimal(1)
    threshold = Decimal(10) ** (-(precision+5))
    n = 0
    while max(abs(s_term), abs(c_term)) > threshold:
        n += 1
        s_term *= -theta*theta / Decimal((2*n)*(2*n+1))
        c_term *= -theta*theta / Decimal((2*n-1)*(2*n))
        sine += s_term
        cosine += c_term
    return cosine, sine


def km_trace(model, initial, precision=80, step_tolerance="1e-24",
               max_iterations=20000):
    """Return actual physical KM iterates for ``PauliLikelihood``.

    For in-plane coordinates u=(u1,u2), g=dL/du, q=sqrt(s), s=1-|u|^2,

        du = q*g - q/(1+q) * (u.g) * u,
        s_next = s * [(1-u.g)^2 - |g|^2].

    The second identity is det(rho_next)=det(rho)*det(R), used to track
    tiny positive boundary gaps without subtracting two nearly equal
    numbers. Its consistency with 1-|u_next|^2 is checked at every step.

    Boundary labels describe finite-precision evidence of an asymptotic
    approach; a computed iterate is never replaced by a pure state.
    """
    if not 40 <= int(precision) <= 120:
        raise ValueError("Use 40 to 120 Decimal digits")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be positive")
    initial = np.asarray(initial, dtype=float)
    if initial.shape != (2,) or not np.all(np.isfinite(initial)):
        raise ValueError("initial must contain two finite coordinates")

    # Guard digits protect trigonometric construction and the purity check.
    with localcontext() as ctx:
        ctx.prec = int(precision) + 10
        one = Decimal(1)
        tolerance = _number(step_tolerance)
        if tolerance <= 0 or tolerance < Decimal(10)**(-(precision//2)):
            raise ValueError("Use a positive step tolerance no smaller than 10^(-precision/2)")
        cosine, sine = _rotation(model.angle_degrees, precision)

        def to_measurement(point):
            x, z = map(_number, point)
            return [cosine*x+sine*z, -sine*x+cosine*z]

        def to_drawing(point):
            u, v = point
            return [float(cosine*u-sine*v), float(sine*u+cosine*v)]

        means = to_measurement(model.center)
        relative = list(map(_number, model.weights_xz))
        if sum(relative) != one:
            raise ValueError("Decimal relative measurement weights must sum exactly to one")
        y_weight = _number(model.y_weight)
        weights = [(one-y_weight)*w for w in relative]
        u = to_measurement(initial)
        s = one-sum(x*x for x in u)
        if s <= 0 or any(abs(x) >= one for x in means):
            raise ValueError("Use a strictly physical initial state and valid outcome means")

        def likelihood(point):
            return sum(w*((one+m)/2*((one+x)/2).ln() +
                          (one-m)/2*((one-x)/2).ln())
                       for w, m, x in zip(weights, means, point)) - y_weight*Decimal(2).ln()

        initial_likelihood = likelihood(u)
        previous_likelihood = initial_likelihood
        points = [to_drawing(u)]
        gaps = [float(s)]
        minimum_increment = Decimal("Infinity")
        identity_error = Decimal(0)
        numerical_slack = Decimal(10)**(-(precision-8))
        last_step = Decimal("Infinity")
        stopped = False
        for iteration in range(1, int(max_iterations)+1):
            g = [w*(m-x)/(one-x*x) for w, m, x in zip(weights, means, u)]
            dot = sum(x*y for x, y in zip(u, g))
            det_r = (one-dot)**2-sum(x*x for x in g)
            if det_r <= 0:
                raise ArithmeticError("The KM update lost positive definiteness")
            q = s.sqrt()
            delta = [q*y-q/(one+q)*dot*x for x, y in zip(u, g)]
            updated = [x+y for x, y in zip(u, delta)]
            updated_s = s*det_r
            identity_error = max(identity_error, abs(one-sum(x*x for x in updated)-updated_s))
            if identity_error > numerical_slack:
                raise ArithmeticError("Precision is insufficient for the purity recurrence")
            if any(abs(x) >= one for x in updated):
                raise ArithmeticError("An iterate left the probability domain")
            current_likelihood = likelihood(updated)
            increment = current_likelihood-previous_likelihood
            if increment < -numerical_slack:
                raise ArithmeticError("Likelihood decreased beyond the arithmetic tolerance")
            minimum_increment = min(minimum_increment, increment)
            u, s = updated, updated_s
            previous_likelihood = current_likelihood
            last_step = sum(x*x for x in delta).sqrt()
            points.append(to_drawing(u))
            gaps.append(float(s))
            if last_step <= tolerance:
                stopped = True
                break

        endpoint = np.asarray(points[-1])
        mle = model.physical_maximum()
        mle_gap = float(model.value(mle)-model.value(endpoint))
        near_mle = np.linalg.norm(endpoint-mle) < 1e-8
        if not stopped:
            classification = "iteration limit reached"
        elif near_mle:
            classification = "boundary MLE" if s < Decimal("1e-12") else "interior MLE"
        elif s < Decimal("1e-18"):
            classification = "nonoptimal boundary approach" if mle_gap > 1e-10 else "boundary approach"
        else:
            classification = "step tolerance reached"

        return KMTrace(np.asarray(points), iteration, endpoint,
                         float(last_step), float(s),
                         float(previous_likelihood-initial_likelihood), mle_gap,
                         classification, float(identity_error),
                         float(minimum_increment), np.asarray(gaps), int(precision))
