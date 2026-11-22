"""Finite witness certificates and independent checks for paper_witness/main.tex.

Run python3 verify_witness.py for checks, or --certificate example.json
to evaluate a supplied finite representation. This is not a paper parser.
Exact finite checks use Fraction; numerical checks are explicitly separate.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


def projection(obs, target, weights):
    groups = {}
    for i, label in enumerate(obs):
        groups.setdefault(label, []).append(i)
    means = {
        label: sum(weights[i] * target[i] for i in indices)
        / sum(weights[i] for i in indices)
        for label, indices in groups.items()
    }
    return tuple(means[label] for label in obs)


def norm2(vector, weights):
    return sum(w * x * x for w, x in zip(weights, vector))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def refines(fine, coarse):
    factor = {}
    for a, b in zip(fine, coarse):
        if a in factor and factor[a] != b:
            return False
        factor[a] = b
    return True


def partitions(n):
    if n == 0:
        yield ()
        return

    def extend(prefix):
        if len(prefix) == n:
            yield prefix
        else:
            for label in range(max(prefix) + 2):
                yield from extend(prefix + (label,))

    yield from extend((0,))


def certificate(data):
    states = data["states"]
    n = len(states)
    if not n or len(set(states)) != n:
        raise ValueError("States must be nonempty and unique.")
    obs = tuple(data["observation"])
    diagnostic = tuple(data["diagnostic"])
    target = tuple(F(str(v)) for v in data["target"])
    weights = tuple(F(str(v)) for v in data["weights"])
    if any(len(x) != n for x in (obs, diagnostic, target, weights)):
        raise ValueError("Every vector must have one entry per state.")
    if any(w <= 0 for w in weights) or sum(weights) != 1:
        raise ValueError("Weights must be positive and sum exactly to one.")
    provenance = data.get("diagnostic_provenance")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("The diagnostic must have a stated provenance.")
    augmented = tuple(zip(obs, diagnostic))
    old_fit = projection(obs, target, weights)
    new_fit = projection(augmented, target, weights)
    old_loss = norm2(subtract(target, old_fit), weights)
    new_loss = norm2(subtract(target, new_fit), weights)
    gain = norm2(subtract(new_fit, old_fit), weights)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if obs[i] == obs[j] and target[i] != target[j]:
                bound = (target[i] - target[j]) ** 2 / (
                    1 / weights[i] + 1 / weights[j]
                )
                pairs.append({"states": [states[i], states[j]],
                              "lower_bound": str(bound)})
    if old_loss - new_loss != gain:
        raise AssertionError("Refinement identity failed.")
    return {
        "domain_size": n,
        "diagnostic_provenance": provenance,
        "preserves_original_observation": refines(augmented, obs),
        "strict_refinement": not refines(obs, augmented),
        "diagnostic_previously_expressible": refines(obs, diagnostic),
        "old_observation_loss": str(old_loss),
        "new_observation_loss": str(new_loss),
        "projection_gain": str(gain),
        "obstruction_pairs": pairs,
        "scope": "Exact finite instance; supplied target and diagnostic required.",
    }


def check_finite():
    n = 4
    ps = list(partitions(n))
    assert len(ps) == 15
    weights = (F(1, 10), F(2, 10), F(3, 10), F(4, 10))
    nested = [(fine, coarse) for fine in ps for coarse in ps
              if refines(fine, coarse)]
    loss_cases = 0
    gain_cases = 0
    universal_cases = 0
    for h in product((-1, 0, 1), repeat=n):
        projections = {obs: projection(obs, h, weights) for obs in ps}
        for obs in ps:
            fit = projections[obs]
            loss = norm2(subtract(h, fit), weights)
            assert (loss == 0) == refines(obs, h)
            # Independent Pythagorean decomposition against the constant model.
            mean = sum(w * v for w, v in zip(weights, h))
            constant = (mean,) * n
            assert norm2(subtract(h, constant), weights) == (
                loss + norm2(subtract(fit, constant), weights)
            )
            for i in range(n):
                for j in range(i + 1, n):
                    if obs[i] == obs[j]:
                        lower = F((h[i] - h[j]) ** 2) / (
                            1 / weights[i] + 1 / weights[j])
                        assert loss >= lower
            loss_cases += 1
        for fine, coarse in nested:
            old, new = projections[coarse], projections[fine]
            assert (
                norm2(subtract(h, old), weights)
                - norm2(subtract(h, new), weights)
                == norm2(subtract(new, old), weights)
            )
            gain_cases += 1
    for coarse in ps:
        for diagnostic in product((0, 1), repeat=n):
            augmented = tuple(zip(coarse, diagnostic))
            assert refines(augmented, coarse)
            assert refines(augmented, diagnostic)
            assert refines(coarse, augmented) == refines(coarse, diagnostic)
            for candidate in ps:
                if refines(candidate, coarse) and refines(candidate, diagnostic):
                    assert refines(candidate, augmented)
                    universal_cases += 1
    example = json.loads(Path(__file__).with_name("example.json").read_text())
    result = certificate(example)
    assert result["old_observation_loss"] == "1"
    assert result["new_observation_loss"] == "0"
    assert result["projection_gain"] == "1"
    for field, value in (("weights", [0, 0, 0, 1]),
                         ("target", [1, 2]),
                         ("diagnostic_provenance", "")):
        invalid = dict(example, **{field: value})
        try:
            certificate(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid certificate input was accepted.")
    return {"loss_cases": loss_cases, "refinement_cases": gain_cases,
            "universal_property_cases": universal_cases}


def check_distributor():
    # Enumerate finite models on a line independently of Hilbert projection.
    points = (-2, -1, 0, 1, 2)
    models = [tuple(p for i, p in enumerate(points) if mask & (1 << i))
              for mask in range(1, 1 << len(points))]
    def distance(x, model):
        return min(abs(x - m) for m in model)
    def excess(a, b):
        return max(distance(x, b) for x in a)
    cases = 0
    for a in models:
        for b in models:
            ab = excess(a, b)
            for f in points:
                for g in points:
                    assert distance(g, b) <= abs(g - f) + distance(f, a) + ab
                    cases += 1
            for c in models:
                assert excess(a, c) <= ab + excess(b, c)
    # Regression: squared distance cannot replace distance in these actions.
    assert abs(0 - 2) ** 2 > abs(0 - 1) ** 2 + abs(1 - 2) ** 2
    return {"distributor_cases": cases, "excess_triangles": len(models) ** 3}


def check_affine_symbolic():
    import sympy as s
    t, a, b = s.symbols("t a b", real=True)
    value = s.integrate((t ** 3 - a - b * t) ** 2, (t, 0, 1))
    derivative = s.integrate((3 * t ** 2 - b) ** 2, (t, 0, 1))
    objective = value + derivative
    l2_fit = s.solve([s.diff(value, a), s.diff(value, b)], (a, b))
    weighted_fit = s.solve([s.diff(objective, a), s.diff(objective, b)], (a, b))
    assert l2_fit[b] == s.Rational(9, 10)
    assert weighted_fit[b] == s.Rational(129, 130)
    gap = s.simplify(objective.subs(l2_fit) - objective.subs(weighted_fit))
    assert gap == s.Rational(3, 325)
    p = l2_fit[a] + l2_fit[b] * t
    q = weighted_fit[a] + weighted_fit[b] * t
    projected_gap = s.integrate((p - q) ** 2 + s.diff(p - q, t) ** 2, (t, 0, 1))
    assert s.simplify(projected_gap - gap) == 0
    return {"l2_slope": str(l2_fit[b]), "weighted_slope": str(weighted_fit[b]),
            "gap": str(gap)}


def check_spectral_numeric():
    import numpy as np
    from scipy.integrate import quad

    def integrate(fn, lo, hi):
        return (quad(lambda t: float(np.real(fn(t))), lo, hi, epsabs=1e-11)[0]
                + 1j * quad(lambda t: float(np.imag(fn(t))), lo, hi, epsabs=1e-11)[0])

    def energy(fn, df, lo, width, lam):
        hi = lo + width
        gram = np.array([[width, 0], [0, width ** 3 / 12 + lam * width]])
        rhs = np.array([integrate(fn, lo, hi),
                        integrate(lambda t: (t - lo - width / 2) * fn(t) + lam * df(t),
                                  lo, hi)])
        a, b = np.linalg.solve(gram, rhs)
        residual = lambda t: (
            abs(fn(t) - a - b * (t - lo - width / 2)) ** 2
            + lam * abs(df(t) - b) ** 2)
        return quad(residual, lo, hi, epsabs=1e-11)[0]

    z, c = 0.3 + 1.7j, 0.8 - 0.4j
    worst = 0.0
    for width in (0.2, 0.7, 1.3):
        for lam in (0.02, 1.0):
            base = energy(lambda t: np.exp(z*t), lambda t: z*np.exp(z*t),
                          0, width, lam)
            scaled = width * energy(
                lambda t: np.exp(z*width*t),
                lambda t: z*width*np.exp(z*width*t), 0, 1, lam/width**2)
            assert np.isclose(base, scaled, rtol=1e-9, atol=1e-12)
            for endpoint in (2.0, 5.0):
                direct = energy(lambda t: c*np.exp(z*t),
                                lambda t: c*z*np.exp(z*t),
                                endpoint-width, width, lam)
                predicted = abs(c)**2*np.exp(2*z.real*(endpoint-width))*base
                relative = abs(direct-predicted)/predicted
                worst = max(worst, relative)
                assert relative < 1e-8
    fn = lambda t: np.exp(t)*np.cos(2*t)
    dfn = lambda t: np.exp(t)*(np.cos(2*t)-2*np.sin(2*t))
    c_star = energy(fn, dfn, 0, 1, 0)
    assert c_star > 0
    for width in (0.5, 0.1, 0.02):
        direct = energy(lambda t: fn(t/width), lambda t: dfn(t/width)/width,
                        0, width, 0)
        assert np.isclose(direct, width*c_star, rtol=1e-9)
    # Equal scalar phase penalties with different Z/3 holonomies.
    phase_distance = lambda a, b: min((a-b) % 3, (b-a) % 3)
    assert phase_distance(0, 1) == phase_distance(0, 2)
    assert phase_distance(1, 2) > 0
    return {"single_mode_max_relative_error": worst, "affine_scaling_c_star": c_star}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path)
    args = parser.parse_args()
    if args.certificate:
        print(json.dumps(certificate(json.loads(args.certificate.read_text())), indent=2))
    else:
        result = {"finite_exact": check_finite(),
                  "enriched_exact": check_distributor(),
                  "projection_symbolic": check_affine_symbolic(),
                  "spectral_numerical": check_spectral_numeric()}
        print(json.dumps(result, indent=2))
        print("PASS: exact finite and symbolic checks; numerical sanity checks.")
        print("These checks do not replace proofs or establish the infinite-spectrum SPTB claim.")


if __name__ == "__main__":
    main()
