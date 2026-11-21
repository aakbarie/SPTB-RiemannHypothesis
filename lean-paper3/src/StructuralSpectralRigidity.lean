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


/-- Arithmetic input to three-prime sampling: after clearing denominators,
a positive power of 5 cannot equal a product of powers of 2 and 3.
This is a factorization statement, not yet the complex-phase theorem. -/
theorem prime_power_nonalias (m n : ℕ) (hn : 0 < n) :
    (5 : ℕ) ^ n * 2 ^ m ≠ 3 ^ m * 2 ^ n := by
  intro heq
  have hp : Nat.Prime 5 := by norm_num
  have h5 : 5 ∣ (5 : ℕ) ^ n := dvd_pow_self 5 (by omega)
  have hl : 5 ∣ (5 : ℕ) ^ n * 2 ^ m := dvd_mul_of_dvd_left h5 _
  have hr : 5 ∣ (3 : ℕ) ^ m * 2 ^ n := heq ▸ hl
  rcases (hp.dvd_mul).mp hr with h3 | h2
  · have hprime : 5 ∣ (3 : ℕ) := hp.dvd_of_dvd_pow h3
    norm_num at hprime
  · have hprime : 5 ∣ (2 : ℕ) := hp.dvd_of_dvd_pow h2
    norm_num at hprime


/-- Genuine prime-specific two-mode observability for real exponents.
The 2- and 3-prime evaluation columns are linearly independent. -/
theorem real_two_prime_minor_nonzero (s t : ℝ) (hst : s ≠ t) :
    Real.exp (-s * Real.log 2) * Real.exp (-t * Real.log 3) ≠
      Real.exp (-t * Real.log 2) * Real.exp (-s * Real.log 3) := by
  intro heq
  have heqExp :
      Real.exp (-s * Real.log 2 + -t * Real.log 3) =
        Real.exp (-t * Real.log 2 + -s * Real.log 3) := by
    simpa only [Real.exp_add] using heq
  have heqArg := Real.exp_injective heqExp
  have hprod : (s - t) * (Real.log 3 - Real.log 2) = 0 := by
    nlinarith [heqArg]
  rcases mul_eq_zero.mp hprod with hst' | hlog
  · exact hst (sub_eq_zero.mp hst')
  · have hl : Real.log 2 < Real.log 3 := Real.log_lt_log (by norm_num) (by norm_num)
    linarith

/-- On the real-exponent two-mode subspace, two actual prime samples
yield strictly positive observed energy. No RH-related fact is used. -/
theorem real_two_prime_frame_positive (s t x y : ℝ)
    (hst : s ≠ t) (hxy : x ≠ 0 ∨ y ≠ 0) :
    0 < (Real.exp (-s * Real.log 2) * x + Real.exp (-t * Real.log 2) * y) ^ 2
      + (Real.exp (-s * Real.log 3) * x + Real.exp (-t * Real.log 3) * y) ^ 2 := by
  apply two_observation_positive
  · exact sub_ne_zero.mpr (real_two_prime_minor_nonzero s t hst)
  · exact hxy


/-- An arithmetic logarithmic nonaliasing identity for the actual prime ratios.
It excludes positive integer rational relations between the two log gaps,
using the kernel-checked unique-factorization lemma above. -/
theorem nat_log_prime_nonalias (m n : ℕ) (hn : 0 < n) :
    (n : ℝ) * Real.log ((5 : ℝ) / 2) ≠
      (m : ℝ) * Real.log ((3 : ℝ) / 2) := by
  intro h
  have hlogpow :
      Real.log (((5 : ℝ) / 2) ^ n) =
        Real.log (((3 : ℝ) / 2) ^ m) := by
    simpa only [Real.log_pow] using h
  have hpow : (((5 : ℝ) / 2) ^ n) = (((3 : ℝ) / 2) ^ m) := by
    have hex := congrArg Real.exp hlogpow
    rwa [Real.exp_log (by positivity), Real.exp_log (by positivity)] at hex
  have hcross : (5 : ℝ) ^ n * (2 : ℝ) ^ m =
      (3 : ℝ) ^ m * (2 : ℝ) ^ n := by
    rw [div_pow, div_pow] at hpow
    exact (div_eq_div_iff (by positivity) (by positivity)).mp hpow
  have hNat : (5 : ℕ) ^ n * 2 ^ m = 3 ^ m * 2 ^ n := by
    exact_mod_cast hcross
  exact (prime_power_nonalias m n hn) hNat

/-- The prime-log ratio cannot be any nonnegative rational with a positive
natural denominator. This is the arithmetic half of the three-prime theorem. -/
theorem nat_prime_log_ratio_nonalias (m n : ℕ) (hn : 0 < n) :
    Real.log ((5 : ℝ) / 2) / Real.log ((3 : ℝ) / 2) ≠
      (m : ℝ) / (n : ℝ) := by
  intro h
  have hd : Real.log ((3 : ℝ) / 2) ≠ 0 :=
    ne_of_gt (Real.log_pos (by norm_num))
  have hn0 : (n : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hn)
  have hcross := (div_eq_div_iff hd hn0).mp h
  have hh : (n : ℝ) * Real.log ((5 : ℝ) / 2) =
      (m : ℝ) * Real.log ((3 : ℝ) / 2) := by
    nlinarith [hcross]
  exact nat_log_prime_nonalias m n hn hh


/-- Nonzero radial scaling of both complex phase observations preserves
exact observational equality. This is one finite-dimensional input to
the information-dimension discussion; it does not formalize entropy limits. -/
theorem nonzero_radial_scaling_preserves_phase_distinction
    (a b : ℝ) (ha : a ≠ 0) (hb : b ≠ 0)
    (z₁ z₂ w₁ w₂ : ℂ) :
    ((a : ℂ) * z₁ = (a : ℂ) * w₁ ∧
     (b : ℂ) * z₂ = (b : ℂ) * w₂) ↔
    (z₁ = w₁ ∧ z₂ = w₂) := by
  have haC : (a : ℂ) ≠ 0 := by exact_mod_cast ha
  have hbC : (b : ℂ) ≠ 0 := by exact_mod_cast hb
  constructor
  · rintro ⟨h₁, h₂⟩
    constructor
    · exact mul_left_cancel₀ haC h₁
    · exact mul_left_cancel₀ hbC h₂
  · rintro ⟨h₁, h₂⟩
    constructor
    · rw [h₁]
    · rw [h₂]

end StructuralSpectralRigidity
