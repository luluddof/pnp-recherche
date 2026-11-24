import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.Transvection
import Mathlib.RingTheory.Nilpotent.Basic
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Order.Interval.Finset.Fin
import Mathlib.LinearAlgebra.Eigenspace.Triangularizable
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv
import Mathlib.LinearAlgebra.LinearIndependent.BaseChange
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure

/-!
# Lemma T (triangular swap lemma) and Theorem 1 (r = 1, 2) — G25-P6

Source: `code/gen24/g24g/notes.md` §1–2.  Conventions: indices `Fin N × Fin N`,
`Swap (a,b) (c,d) = [a = d ∧ b = c]`, `K' = ∑_l E_l ⊗ₖ C_l` (Mathlib Kronecker:
`(E ⊗ₖ C) (a,b) (c,d) = E a c * C b d`, so `E` acts on the first tensor factor).

* `card_Om_le_rank_of_upper` : core, upper-triangular `E_l`: `|{(b,j) : j < b}| ≤ rank (Swap - K')`
  (via `(Swap - K') Swap = 1 - Θ`, `Θ = K' Swap`, and `Θ` strictly triangular, hence nilpotent,
  on the lower pairs ordered by the second coordinate).
* `lemmaT` : general form (invertible `X, Y`, every `X E_l Y⁻¹` upper triangular), any field.
* `swap_sub_kron_rank` : r = 1, any field.
* `lemmaP` : two matrices over an algebraically closed field are simultaneously equivalent to upper
  triangular ones (direct induction: common direction `v` with `E₁v, E₂v ∈ K w`, then deflate;
  no Kronecker canonical form needed).
* `swap_sub_two_kron_rank` (alg. closed), `swap_sub_two_kron_rank'` : r = 2, any field
  (through `AlgebraicClosure K`, `rank_map_algebraMap_le`).
-/

open Matrix
open scoped Kronecker

namespace Triang

variable {K : Type*} [Field K]

/-- `Swap (a,b) (c,d) = [a = d ∧ b = c]`, i.e. `e_x ⊗ e_y ↦ e_y ⊗ e_x`. -/
def Swap (K : Type*) [Field K] (N : ℕ) : Matrix (Fin N × Fin N) (Fin N × Fin N) K :=
  Matrix.of fun a b => if a.1 = b.2 ∧ a.2 = b.1 then 1 else 0

variable {N : ℕ}

theorem mul_swap_apply (M : Matrix (Fin N × Fin N) (Fin N × Fin N) K) (p q : Fin N × Fin N) :
    (M * Swap K N) p q = M p (q.2, q.1) := by
  rw [Matrix.mul_apply, Finset.sum_eq_single (q.2, q.1)]
  · simp [Swap]
  · intro r _ hr
    have : ¬ (r.1 = q.2 ∧ r.2 = q.1) := by
      rintro ⟨h1, h2⟩; exact hr (Prod.ext h1 h2)
    simp [Swap, this]
  · simp

theorem swap_mul_apply (M : Matrix (Fin N × Fin N) (Fin N × Fin N) K) (p q : Fin N × Fin N) :
    (Swap K N * M) p q = M (p.2, p.1) q := by
  rw [Matrix.mul_apply, Finset.sum_eq_single (p.2, p.1)]
  · simp [Swap]
  · intro r _ hr
    have : ¬ (p.1 = r.2 ∧ p.2 = r.1) := by
      rintro ⟨h1, h2⟩; exact hr (Prod.ext h2.symm h1.symm)
    simp [Swap, this]
  · simp

theorem swap_mul_swap : Swap K N * Swap K N = 1 := by
  ext p q
  rw [mul_swap_apply]
  simp only [Swap, Matrix.of_apply, Matrix.one_apply]
  by_cases h : p = q
  · subst h; simp
  · have : ¬ (p.1 = q.1 ∧ p.2 = q.2) := fun ⟨h1, h2⟩ => h (Prod.ext h1 h2)
    simp [this, h]

theorem kron_mul_swap (A B : Matrix (Fin N) (Fin N) K) :
    (A ⊗ₖ B) * Swap K N = Swap K N * (B ⊗ₖ A) := by
  ext p q
  rw [mul_swap_apply, swap_mul_apply]
  simp [Matrix.kroneckerMap_apply, mul_comm]

/-- The "lower pairs" `{(b, j) : j < b}`. -/
abbrev Om (N : ℕ) := {p : Fin N × Fin N // p.2 < p.1}

theorem card_Om : Fintype.card (Om N) = N * (N - 1) / 2 := by
  rw [Fintype.card_subtype, Finset.card_eq_sum_ones, Finset.sum_filter, Fintype.sum_prod_type]
  have h : ∀ i : Fin N, (∑ j : Fin N, if j < i then 1 else 0) = (i : ℕ) := by
    intro i
    rw [← Finset.sum_filter, ← Finset.card_eq_sum_ones]
    have : (Finset.univ.filter fun j : Fin N => j < i) = Finset.Iio i := by
      ext j; simp
    rw [this, Fin.card_Iio]
  simp_rw [h]
  rw [Fin.sum_univ_eq_sum_range (fun i => i) N, Finset.sum_range_id]

/-- Core of Lemma T: upper-triangular `E_l`, arbitrary `C_l`. -/
theorem card_Om_le_rank_of_upper {ι : Type*} [Fintype ι]
    (E C : ι → Matrix (Fin N) (Fin N) K)
    (hE : ∀ l i j, j < i → E l i j = 0) :
    Fintype.card (Om N) ≤ (Swap K N - ∑ l, E l ⊗ₖ C l).rank := by
  classical
  set Kp := ∑ l, E l ⊗ₖ C l with hKp
  set Θ := Kp * Swap K N with hΘ
  have h1 : (Swap K N - Kp) * Swap K N = 1 - Θ := by
    rw [Matrix.sub_mul, swap_mul_swap]
  have hΘapp : ∀ a b c d : Fin N, Θ (a, b) (c, d) = ∑ l, E l a d * C l b c := by
    intro a b c d
    rw [hΘ, mul_swap_apply, hKp, Matrix.sum_apply]
    simp [Matrix.kroneckerMap_apply]
  let T : Matrix (Om N) (Om N) K := Θ.submatrix Subtype.val Subtype.val
  have hT : ∀ p q : Om N, T p q ≠ 0 → p.1.2 < q.1.2 := by
    intro p q hpq
    by_contra hlt
    rw [not_lt] at hlt
    apply hpq
    show Θ p.1 q.1 = 0
    obtain ⟨⟨a, b⟩, hab⟩ := p
    obtain ⟨⟨c, d⟩, hcd⟩ := q
    simp only at hab hcd hlt
    rw [hΘapp]
    refine Finset.sum_eq_zero fun l _ => ?_
    rw [hE l a d (lt_of_le_of_lt hlt hab), zero_mul]
  have hpow : ∀ k : ℕ, ∀ p q : Om N, (T ^ k) p q ≠ 0 → (p.1.2 : ℕ) + k ≤ q.1.2 := by
    intro k
    induction k with
    | zero =>
      intro p q h
      rw [pow_zero, Matrix.one_apply] at h
      by_cases hpq : p = q
      · subst hpq; simp
      · simp [hpq] at h
    | succ k ih =>
      intro p q h
      rw [pow_succ, Matrix.mul_apply] at h
      obtain ⟨r, _, hr⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
      have h1 := ih p r (left_ne_zero_of_mul hr)
      have h2 := hT r q (right_ne_zero_of_mul hr)
      have : (r.1.2 : ℕ) < q.1.2 := h2
      omega
  have hnil : IsNilpotent T := by
    refine ⟨N, ?_⟩
    ext p q
    by_contra h
    have := hpow N p q h
    have := q.1.2.isLt
    simp only [Matrix.zero_apply] at *
    omega
  have hsub : (1 - Θ).submatrix Subtype.val Subtype.val = 1 - T := by
    ext p q
    simp only [Matrix.submatrix_apply, Matrix.sub_apply, Matrix.one_apply, T]
    by_cases hpq : p = q
    · subst hpq; simp
    · have : p.1 ≠ q.1 := fun h => hpq (Subtype.ext h)
      simp [hpq, this]
  calc Fintype.card (Om N) = (1 - T).rank := (Matrix.rank_of_isUnit _ hnil.isUnit_one_sub).symm
    _ = ((1 - Θ).submatrix Subtype.val Subtype.val).rank := by rw [hsub]
    _ ≤ (1 - Θ).rank := Matrix.rank_submatrix_le _ _ _
    _ = ((Swap K N - Kp) * Swap K N).rank := by rw [h1]
    _ ≤ (Swap K N - Kp).rank := Matrix.rank_mul_le_left _ _


theorem isUnit_det_kron {A B : Matrix (Fin N) (Fin N) K} (hA : IsUnit A.det) (hB : IsUnit B.det) :
    IsUnit (A ⊗ₖ B).det := by
  rw [det_kronecker]
  exact (hA.pow _).mul (hB.pow _)

/-- The symmetry `M ↦ (X ⊗ Y) M (Y⁻¹ ⊗ X⁻¹)`: it fixes `Swap` and sends `E ⊗ C` to
`X E Y⁻¹ ⊗ Y C X⁻¹`. -/
theorem conj_swap_sub {ι : Type*} [Fintype ι] (E C : ι → Matrix (Fin N) (Fin N) K)
    (X Y : Matrix (Fin N) (Fin N) K) (hX : IsUnit X.det) (hY : IsUnit Y.det) :
    (X ⊗ₖ Y) * (Swap K N - ∑ l, E l ⊗ₖ C l) * (Y⁻¹ ⊗ₖ X⁻¹) =
      Swap K N - ∑ l, (X * E l * Y⁻¹) ⊗ₖ (Y * C l * X⁻¹) := by
  rw [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_sum, Matrix.sum_mul]
  congr 1
  · rw [kron_mul_swap, Matrix.mul_assoc, ← mul_kronecker_mul, Matrix.mul_nonsing_inv _ hY,
      Matrix.mul_nonsing_inv _ hX, one_kronecker_one, Matrix.mul_one]
  · refine Finset.sum_congr rfl fun l _ => ?_
    rw [← mul_kronecker_mul, ← mul_kronecker_mul]

/-- **Lemma T** (triangular swap lemma). If some invertible `X, Y` make every `X E_l Y⁻¹`
upper triangular, then for all `C_l`, `rank (Swap - ∑_l E_l ⊗ C_l) ≥ N(N-1)/2`. -/
theorem lemmaT {ι : Type*} [Fintype ι] (E C : ι → Matrix (Fin N) (Fin N) K)
    (X Y : Matrix (Fin N) (Fin N) K) (hX : IsUnit X.det) (hY : IsUnit Y.det)
    (hT : ∀ l i j, j < i → (X * E l * Y⁻¹) i j = 0) :
    N * (N - 1) / 2 ≤ (Swap K N - ∑ l, E l ⊗ₖ C l).rank := by
  have hP := isUnit_det_kron hX hY
  have hQ := isUnit_det_kron (Matrix.isUnit_nonsing_inv_det _ hY)
    (Matrix.isUnit_nonsing_inv_det _ hX)
  rw [← card_Om, ← Matrix.rank_mul_eq_right_of_isUnit_det _ _ hP,
    ← Matrix.rank_mul_eq_left_of_isUnit_det _ _ hQ, conj_swap_sub E C X Y hX hY]
  exact card_Om_le_rank_of_upper _ _ hT

/-- Every square matrix over a field is equivalent to an upper triangular (indeed diagonal) one. -/
theorem exists_equiv_upper (E : Matrix (Fin N) (Fin N) K) :
    ∃ X Y : Matrix (Fin N) (Fin N) K, IsUnit X.det ∧ IsUnit Y.det ∧
      ∀ i j, j < i → (X * E * Y⁻¹) i j = 0 := by
  obtain ⟨L, L', D, h⟩ := Matrix.Pivot.exists_list_transvec_mul_mul_list_transvec_eq_diagonal E
  set P := (L'.map TransvectionStruct.toMatrix).prod
  have hP : IsUnit P.det := by rw [TransvectionStruct.det_toMatrix_prod]; exact isUnit_one
  refine ⟨(L.map TransvectionStruct.toMatrix).prod, P⁻¹, by rw [TransvectionStruct.det_toMatrix_prod]; exact isUnit_one,
    Matrix.isUnit_nonsing_inv_det _ hP, ?_⟩
  intro i j hij
  rw [Matrix.nonsing_inv_nonsing_inv _ hP, h, Matrix.diagonal_apply_ne _ (ne_of_gt hij)]

/-- **Corollary (r = 1)**, any field: `rank (Swap - E ⊗ C) ≥ N(N-1)/2`. -/
theorem swap_sub_kron_rank (E C : Matrix (Fin N) (Fin N) K) :
    N * (N - 1) / 2 ≤ (Swap K N - E ⊗ₖ C).rank := by
  obtain ⟨X, Y, hX, hY, hT⟩ := exists_equiv_upper E
  have := lemmaT (fun _ : Unit => E) (fun _ : Unit => C) X Y hX hY (fun _ => hT)
  simpa using this

/-! ## Lemma P (pairs over an algebraically closed field) -/

/-- For `v ≠ 0` there is an invertible matrix whose column `0` is `v`. -/
theorem exists_col0 {n : ℕ} (v : Fin (n + 1) → K) (hv : v ≠ 0) :
    ∃ Z : Matrix (Fin (n + 1)) (Fin (n + 1)) K, IsUnit Z.det ∧ ∀ i, Z i 0 = v i := by
  classical
  obtain ⟨k, hk⟩ : ∃ k, v k ≠ 0 := by
    by_contra h
    push Not at h
    exact hv (funext h)
  set σ := Equiv.swap (0 : Fin (n + 1)) k
  set u : Fin (n + 1) → K := fun i => v (σ i)
  let Nm : Matrix (Fin (n + 1)) (Fin (n + 1)) K :=
    Matrix.of fun i j => if j = 0 then u i else (1 : Matrix (Fin (n + 1)) (Fin (n + 1)) K) i j
  have hlow : Nm.IsLowerTriangular := by
    intro i j hij
    have hij' : i < j := hij
    have hj : j ≠ 0 := by
      rintro rfl; exact absurd hij' (Fin.not_lt_zero _)
    simp [Nm, hj, Matrix.one_apply, ne_of_lt hij']
  have hdet : Nm.det = v k := by
    rw [det_of_isLowerTriangular _ hlow, Fin.prod_univ_succ]
    simp [Nm, u, σ, Fin.succ_ne_zero]
  refine ⟨Nm.submatrix σ id, ?_, ?_⟩
  · rw [det_permute, hdet]
    refine IsUnit.mul ?_ (isUnit_iff_ne_zero.2 hk)
    rcases Int.units_eq_one_or (Equiv.Perm.sign σ) with h | h <;> simp [h]
  · intro i
    simp [Nm, u, σ]


/-- `diag(1, A)`. -/
def ext1 {n : ℕ} (A : Matrix (Fin n) (Fin n) K) : Matrix (Fin (n + 1)) (Fin (n + 1)) K :=
  Matrix.of fun i j => Fin.cases (Fin.cases 1 (fun _ => 0) j) (fun i' => Fin.cases 0 (fun j' => A i' j') j) i

@[simp] theorem ext1_zero_zero {n : ℕ} (A : Matrix (Fin n) (Fin n) K) : ext1 A 0 0 = 1 := rfl
@[simp] theorem ext1_zero_succ {n : ℕ} (A : Matrix (Fin n) (Fin n) K) (j : Fin n) : ext1 A 0 j.succ = 0 := rfl
@[simp] theorem ext1_succ_zero {n : ℕ} (A : Matrix (Fin n) (Fin n) K) (i : Fin n) : ext1 A i.succ 0 = 0 := rfl
@[simp] theorem ext1_succ_succ {n : ℕ} (A : Matrix (Fin n) (Fin n) K) (i j : Fin n) :
    ext1 A i.succ j.succ = A i j := rfl

theorem det_ext1 {n : ℕ} (A : Matrix (Fin n) (Fin n) K) : (ext1 A).det = A.det := by
  rw [det_succ_column_zero, Fin.sum_univ_succ]
  simp only [ext1_zero_zero, ext1_succ_zero, mul_zero, zero_mul, Finset.sum_const_zero, add_zero,
    Fin.val_zero, pow_zero, one_mul, Fin.succAbove_zero]
  congr 1

theorem ext1_mul_mul_succ_zero {n : ℕ} (A B : Matrix (Fin n) (Fin n) K)
    (M : Matrix (Fin (n + 1)) (Fin (n + 1)) K) (hM : ∀ a : Fin n, M a.succ 0 = 0) (i : Fin n) :
    (ext1 A * M * ext1 B) i.succ 0 = 0 := by
  simp [Matrix.mul_apply, Fin.sum_univ_succ, hM]

theorem ext1_mul_mul_succ_succ {n : ℕ} (A B : Matrix (Fin n) (Fin n) K)
    (M : Matrix (Fin (n + 1)) (Fin (n + 1)) K) (i j : Fin n) :
    (ext1 A * M * ext1 B) i.succ j.succ = (A * M.submatrix Fin.succ Fin.succ * B) i j := by
  simp [Matrix.mul_apply, Fin.sum_univ_succ]


/-- A common direction: `v ≠ 0` with `E₁ v, E₂ v ∈ K w`, `w ≠ 0`. -/
theorem exists_common_dir [IsAlgClosed K] {n : ℕ} (E1 E2 : Matrix (Fin (n + 1)) (Fin (n + 1)) K) :
    ∃ v w : Fin (n + 1) → K, v ≠ 0 ∧ w ≠ 0 ∧ ∃ c1 c2 : K, E1 *ᵥ v = c1 • w ∧ E2 *ᵥ v = c2 • w := by
  classical
  by_cases hd : E1.det = 0
  · obtain ⟨v, hv, hEv⟩ := Matrix.exists_mulVec_eq_zero_iff.2 hd
    by_cases h2 : E2 *ᵥ v = 0
    · refine ⟨v, Pi.single 0 1, hv, ?_, 0, 0, by simp [hEv], by simp [h2]⟩
      intro h
      have := congrFun h 0
      simp at this
    · exact ⟨v, E2 *ᵥ v, hv, h2, 0, 1, by simp [hEv], by simp⟩
  · have hu : IsUnit E1.det := isUnit_iff_ne_zero.2 hd
    obtain ⟨c, hc⟩ := Module.End.exists_eigenvalue (Matrix.toLin' (E1⁻¹ * E2))
    obtain ⟨v, hv⟩ := hc.exists_hasEigenvector
    have hv0 : v ≠ 0 := hv.2
    have hvc : (E1⁻¹ * E2) *ᵥ v = c • v := by
      have := Module.End.mem_eigenspace_iff.1 hv.1
      simpa [Matrix.toLin'_apply] using this
    have hE2 : E2 *ᵥ v = c • (E1 *ᵥ v) := by
      have := congrArg (fun x => E1 *ᵥ x) hvc
      simp only [Matrix.mulVec_mulVec, ← Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hu,
        Matrix.one_mul, Matrix.mulVec_smul] at this
      exact this
    refine ⟨v, E1 *ᵥ v, hv0, ?_, 1, c, by simp, hE2⟩
    intro h
    apply hv0
    have := congrArg (fun x => E1⁻¹ *ᵥ x) h
    simpa [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hu] using this

/-- Lemma P, in the form `X E_i Z` upper triangular. -/
theorem lemmaP_aux [IsAlgClosed K] : ∀ (n : ℕ) (E1 E2 : Matrix (Fin n) (Fin n) K),
    ∃ X Z : Matrix (Fin n) (Fin n) K, IsUnit X.det ∧ IsUnit Z.det ∧
      (X * E1 * Z).IsUpperTriangular ∧ (X * E2 * Z).IsUpperTriangular
  | 0, _, _ => ⟨1, 1, by simp, by simp, fun i => i.elim0, fun i => i.elim0⟩
  | n + 1, E1, E2 => by
    classical
    obtain ⟨v, w, hv, hw, c1, c2, h1, h2⟩ := exists_common_dir E1 E2
    obtain ⟨Y0, hY0, hY0c⟩ := exists_col0 v hv
    obtain ⟨Z0, hZ0, hZ0c⟩ := exists_col0 w hw
    have hYv : Y0 *ᵥ Pi.single 0 1 = v := by
      rw [Matrix.mulVec_single_one]; funext i; exact hY0c i
    have hZw : Z0 *ᵥ Pi.single 0 1 = w := by
      rw [Matrix.mulVec_single_one]; funext i; exact hZ0c i
    have hcol : ∀ (E : Matrix (Fin (n + 1)) (Fin (n + 1)) K) (c : K), E *ᵥ v = c • w →
        ∀ a : Fin n, (Z0⁻¹ * E * Y0) a.succ 0 = 0 := by
      intro E c hE a
      have : (Z0⁻¹ * E * Y0) a.succ 0 = ((Z0⁻¹ * E * Y0) *ᵥ Pi.single 0 1) a.succ := by
        rw [Matrix.mulVec_single_one]; rfl
      rw [this, ← Matrix.mulVec_mulVec, hYv, ← Matrix.mulVec_mulVec, hE, Matrix.mulVec_smul,
        ← hZw, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hZ0, Matrix.one_mulVec]
      simp [Fin.succ_ne_zero]
    obtain ⟨X', Z', hX', hZ', hT1, hT2⟩ := lemmaP_aux n
      ((Z0⁻¹ * E1 * Y0).submatrix Fin.succ Fin.succ) ((Z0⁻¹ * E2 * Y0).submatrix Fin.succ Fin.succ)
    have key : ∀ (E : Matrix (Fin (n + 1)) (Fin (n + 1)) K) (c : K), E *ᵥ v = c • w →
        ((X' * (Z0⁻¹ * E * Y0).submatrix Fin.succ Fin.succ * Z').IsUpperTriangular) →
        (ext1 X' * Z0⁻¹ * E * (Y0 * ext1 Z')).IsUpperTriangular := by
      intro E c hE hT i j hij
      have hij' : j < i := hij
      have hre : ext1 X' * Z0⁻¹ * E * (Y0 * ext1 Z') = ext1 X' * (Z0⁻¹ * E * Y0) * ext1 Z' := by
        simp only [Matrix.mul_assoc]
      rw [hre]
      induction i using Fin.cases with
      | zero => exact absurd hij' (Fin.not_lt_zero _)
      | succ i =>
        induction j using Fin.cases with
        | zero => exact ext1_mul_mul_succ_zero _ _ _ (hcol E c hE) i
        | succ j =>
          rw [ext1_mul_mul_succ_succ]
          exact hT (Fin.succ_lt_succ_iff.1 hij')
    refine ⟨ext1 X' * Z0⁻¹, Y0 * ext1 Z', ?_, ?_, key E1 c1 h1 hT1, key E2 c2 h2 hT2⟩
    · rw [Matrix.det_mul, det_ext1]
      exact hX'.mul (Matrix.isUnit_nonsing_inv_det _ hZ0)
    · rw [Matrix.det_mul, det_ext1]
      exact hY0.mul hZ'

/-- **Lemma P** (pencils): over an algebraically closed field, any two `N × N` matrices are
simultaneously equivalent to upper triangular ones. -/
theorem lemmaP [IsAlgClosed K] (E1 E2 : Matrix (Fin N) (Fin N) K) :
    ∃ X Y : Matrix (Fin N) (Fin N) K, IsUnit X.det ∧ IsUnit Y.det ∧
      (∀ i j, j < i → (X * E1 * Y⁻¹) i j = 0) ∧ (∀ i j, j < i → (X * E2 * Y⁻¹) i j = 0) := by
  obtain ⟨X, Z, hX, hZ, h1, h2⟩ := lemmaP_aux N E1 E2
  refine ⟨X, Z⁻¹, hX, Matrix.isUnit_nonsing_inv_det _ hZ, ?_, ?_⟩
  · intro i j hij; rw [Matrix.nonsing_inv_nonsing_inv _ hZ]; exact h1 hij
  · intro i j hij; rw [Matrix.nonsing_inv_nonsing_inv _ hZ]; exact h2 hij

/-- **Theorem 1, r = 2** (algebraically closed field):
`rank (Swap - E₁ ⊗ C₁ - E₂ ⊗ C₂) ≥ N(N-1)/2`. -/
theorem swap_sub_two_kron_rank [IsAlgClosed K] (E1 C1 E2 C2 : Matrix (Fin N) (Fin N) K) :
    N * (N - 1) / 2 ≤ (Swap K N - (E1 ⊗ₖ C1 + E2 ⊗ₖ C2)).rank := by
  obtain ⟨X, Y, hX, hY, h1, h2⟩ := lemmaP E1 E2
  have := lemmaT ![E1, E2] ![C1, C2] X Y hX hY (by
    intro l; fin_cases l
    · exact h1
    · exact h2)
  simpa [Fin.sum_univ_two] using this


/-- Rank does not grow under a field extension. -/
theorem rank_map_algebraMap_le {L : Type*} [Field L] [Algebra K L] {m n : Type*} [Fintype m]
    [Fintype n] (A : Matrix m n K) : (A.map (algebraMap K L)).rank ≤ A.rank := by
  classical
  rw [← finrank_span_range_row, ← finrank_span_range_row]
  set S := Set.range (A.map (algebraMap K L)).row
  obtain ⟨b, hbS, hspan, hli⟩ := exists_linearIndependent L S
  have hfin : b.Finite := (Set.finite_range _).subset hbS
  let : Fintype b := hfin.fintype
  have h1 : Module.finrank L (Submodule.span L S) = Fintype.card b := by
    rw [← hspan, ← finrank_span_eq_card hli, Subtype.range_coe]
  have hc : ∀ x : b, ∃ i, (A.map (algebraMap K L)).row i = x := fun x => hbS x.2
  choose c hc using hc
  let g : b → n → K := fun x => A.row (c x)
  have hg' : (fun x => algebraMap K L ∘ g x) = ((↑) : b → n → L) := by
    funext x
    rw [← hc x]
    rfl
  have hg : LinearIndependent K g := by
    rw [← linearIndependent_algebraMap_comp_iff (S := L), hg']
    exact hli
  rw [h1, ← finrank_span_eq_card hg]
  apply Submodule.finrank_mono
  apply Submodule.span_mono
  rintro _ ⟨x, rfl⟩
  exact ⟨c x, rfl⟩

/-- **Theorem 1, r = 2**, any field: `rank (Swap - E₁ ⊗ C₁ - E₂ ⊗ C₂) ≥ N(N-1)/2`. -/
theorem swap_sub_two_kron_rank' (E1 C1 E2 C2 : Matrix (Fin N) (Fin N) K) :
    N * (N - 1) / 2 ≤ (Swap K N - (E1 ⊗ₖ C1 + E2 ⊗ₖ C2)).rank := by
  let L := AlgebraicClosure K
  let f := algebraMap K L
  have hmap : (Swap K N - (E1 ⊗ₖ C1 + E2 ⊗ₖ C2)).map f =
      Swap L N - (E1.map f ⊗ₖ C1.map f + E2.map f ⊗ₖ C2.map f) := by
    ext p q
    simp only [Matrix.map_apply, Matrix.sub_apply, Matrix.add_apply, Matrix.kroneckerMap_apply,
      Swap, Matrix.of_apply, map_sub, map_add, map_mul]
    split_ifs <;> simp
  calc N * (N - 1) / 2 ≤ (Swap L N - (E1.map f ⊗ₖ C1.map f + E2.map f ⊗ₖ C2.map f)).rank :=
        swap_sub_two_kron_rank _ _ _ _
    _ = ((Swap K N - (E1 ⊗ₖ C1 + E2 ⊗ₖ C2)).map f).rank := by rw [hmap]
    _ ≤ _ := rank_map_algebraMap_le _

end Triang

#print axioms Triang.card_Om_le_rank_of_upper
#print axioms Triang.lemmaT
#print axioms Triang.swap_sub_kron_rank
#print axioms Triang.lemmaP
#print axioms Triang.swap_sub_two_kron_rank
#print axioms Triang.swap_sub_two_kron_rank'
