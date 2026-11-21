import Mathlib

/-!
A two-dimensional operator-metric rigidity lemma for Paper 3.
A = [[eta, -gamma], [gamma, eta]], G = [[a,c], [c,b]].
The diagonal entries of Aᵀ G + G A are encoded below.
This establishes generic linear algebra, not an arithmetic metric or RH.
-/

namespace StructuralSpectralRigidity

def defect11 (eta gamma a c : ℝ) : ℝ := 2 * eta * a + 2 * gamma * c
def defect22 (eta gamma b c : ℝ) : ℝ := 2 * eta * b - 2 * gamma * c

theorem defect_trace (eta gamma a b c : ℝ) :
    defect11 eta gamma a c + defect22 eta gamma b c = 2 * eta * (a + b) := by
  dsimp [defect11, defect22]
  ring

theorem positive_metric_excludes_drift
    (eta gamma a b c : ℝ) (ha : 0 < a) (hb : 0 < b)
    (h11 : defect11 eta gamma a c = 0)
    (h22 : defect22 eta gamma b c = 0) : eta = 0 := by
  have htrace : 2 * eta * (a + b) = 0 := by
    rw [← defect_trace eta gamma a b c, h11, h22]
    ring
  have hsum : a + b ≠ 0 := ne_of_gt (add_pos ha hb)
  have hprod : eta * (a + b) = 0 := by
    nlinarith [htrace]
  rcases mul_eq_zero.mp hprod with hη | hab
  · exact hη
  · exact False.elim (hsum hab)

theorem quantitative_defect_bound
    (eta gamma a b c : ℝ) (ha : 0 < a) (hb : 0 < b) :
    |eta| ≤
      (|defect11 eta gamma a c| + |defect22 eta gamma b c|) / (2 * (a + b)) := by
  have hpositive : 0 < (2 : ℝ) * (a + b) := by
    have hs : 0 < a + b := add_pos ha hb
    positivity
  have hfactor :
      |defect11 eta gamma a c + defect22 eta gamma b c| =
        |eta| * (2 * (a + b)) := by
    rw [defect_trace]
    calc
      |2 * eta * (a + b)| = |eta * (2 * (a + b))| := by congr 1 <;> ring
      _ = |eta| * |2 * (a + b)| := abs_mul _ _
      _ = |eta| * (2 * (a + b)) := by rw [abs_of_pos hpositive]
  apply (le_div_iff₀ hpositive).2
  rw [← hfactor]
  exact abs_add_le _ _

theorem rotation_preserves_diagonal_metric (gamma : ℝ) :
    defect11 0 gamma 1 0 = 0 ∧ defect22 0 gamma 1 0 = 0 := by
  constructor <;> norm_num [defect11, defect22]

def casimirImag (beta gamma : ℝ) : ℝ := gamma * (1 - 2 * beta)

theorem critical_line_of_casimir_reality
    (beta gamma : ℝ) (hgamma : gamma ≠ 0)
    (hreal : casimirImag beta gamma = 0) : beta = (1 : ℝ) / 2 := by
  dsimp [casimirImag] at hreal
  rcases mul_eq_zero.mp hreal with hz | hb
  · exact False.elim (hgamma hz)
  · linarith


/-- Real part of the arithmetic Casimir eigenvalue rho*(1-rho). -/
def casimirReal (beta gamma : ℝ) : ℝ := beta * (1 - beta) + gamma ^ 2

/-- Conditional critical-line result: the existence of a positive
compatible metric for the phase generator of rho*(1-rho) forces
beta=1/2 when gamma is nonzero. This constructs no such metric. -/
theorem critical_line_from_compatible_phase_metric
    (beta gamma a b c : ℝ) (ha : 0 < a) (hb : 0 < b) (hgamma : gamma ≠ 0)
    (h11 : defect11 (-casimirImag beta gamma) (casimirReal beta gamma) a c = 0)
    (h22 : defect22 (-casimirImag beta gamma) (casimirReal beta gamma) b c = 0) :
    beta = (1 : ℝ) / 2 := by
  have hz := positive_metric_excludes_drift
    (-casimirImag beta gamma) (casimirReal beta gamma) a b c ha hb h11 h22
  have hreal : casimirImag beta gamma = 0 := by linarith
  exact critical_line_of_casimir_reality beta gamma hgamma hreal


/-- If a two-observation real sampling minor is nonzero, the measurements
separate every nonzero coefficient vector. This is generic finite-dimensional
linear algebra; the number-theoretic nonzero minor for p=2,3,5 is separate. -/
theorem two_observation_injective
    (a b c d x y : ℝ)
    (hdet : a * d - b * c ≠ 0)
    (h1 : a * x + b * y = 0)
    (h2 : c * x + d * y = 0) : x = 0 ∧ y = 0 := by
  have hx : (a * d - b * c) * x = 0 := by
    calc
      (a * d - b * c) * x =
          d * (a * x + b * y) - b * (c * x + d * y) := by ring
      _ = 0 := by rw [h1, h2]; ring
  have hy : (a * d - b * c) * y = 0 := by
    calc
      (a * d - b * c) * y =
          a * (c * x + d * y) - c * (a * x + b * y) := by ring
      _ = 0 := by rw [h1, h2]; ring
  constructor
  · rcases mul_eq_zero.mp hx with h | h
    · exact False.elim (hdet h)
    · exact h
  · rcases mul_eq_zero.mp hy with h | h
    · exact False.elim (hdet h)
    · exact h

/-- A nonzero observation minor yields positive squared-observation energy. -/
theorem two_observation_positive
    (a b c d x y : ℝ)
    (hdet : a * d - b * c ≠ 0)
    (hxy : x ≠ 0 ∨ y ≠ 0) :
    0 < (a * x + b * y) ^ 2 + (c * x + d * y) ^ 2 := by
  by_contra hnot
  have hle : (a * x + b * y) ^ 2 + (c * x + d * y) ^ 2 ≤ 0 :=
    le_of_not_gt hnot
  have hsq1 : (a * x + b * y) ^ 2 = 0 := by
    nlinarith [sq_nonneg (c * x + d * y)]
  have hsq2 : (c * x + d * y) ^ 2 = 0 := by
    nlinarith [sq_nonneg (a * x + b * y)]
  have h1 : a * x + b * y = 0 := by
    have hm : (a * x + b * y) * (a * x + b * y) = 0 := by
      simpa [pow_two] using hsq1
    rcases mul_eq_zero.mp hm with h | h <;> exact h
  have h2 : c * x + d * y = 0 := by
    have hm : (c * x + d * y) * (c * x + d * y) = 0 := by
      simpa [pow_two] using hsq2
    rcases mul_eq_zero.mp hm with h | h <;> exact h
  obtain ⟨hx, hy⟩ := two_observation_injective a b c d x y hdet h1 h2
  rcases hxy with h | h
  · exact h hx
  · exact h hy

/-- Additional observations cannot erase a previously strict frame bound. -/
theorem three_observation_positive
    (a b c d e f x y : ℝ)
    (hdet : a * d - b * c ≠ 0)
    (hxy : x ≠ 0 ∨ y ≠ 0) :
    0 < (a * x + b * y) ^ 2 + (c * x + d * y) ^ 2 +
      (e * x + f * y) ^ 2 := by
  have hp := two_observation_positive a b c d x y hdet hxy
  nlinarith [sq_nonneg (e * x + f * y)]

end StructuralSpectralRigidity
