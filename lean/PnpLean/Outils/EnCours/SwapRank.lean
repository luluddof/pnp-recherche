import Mathlib.Tactic
import Mathlib.LinearAlgebra.Matrix.Rank

/-!
# Theorem C (swap with a spectator: rank vs operator-Schmidt rank) — G20g

Source: `code/gen19/g19e/notes.md` §1 (statement and proofs) and §6 (Lean statement, kept verbatim:
`swapW`, `realign`, `swap_rank_osr_tradeoff`, `swap_rank_osr_product`).

* `part_a`  : `k Q n ≤ r + k² Q b` for `k ≤ m` (compression to `k` values of the `E` index on both
  sides; the compressed `B` factors through `ρ × k × k × Q`, `ρ = rank (realign B)`).
* `part_b`  : `k n Q ≤ b Q + k² r` for `k ≤ m` (`Ψ = ker`, of codim `≤ k r`; `1_k ⊗ Ψ`).
* `part_a'`, `part_b'` : the mirror forms (`k ≤ n`), by transposition.
* `part_c`  : `b < m n → Q ≤ r` (contraction of site 0 with a functional killing `B`).
* `swap_rank_osr_product` : `m = n = M`, `1 ≤ b < M²` ⇒ `2 Q M² ≤ 9 r b`.
-/

open Matrix Module

namespace SwapRank

section Rank

variable {F : Type*} [Field F]

theorem rank_add_le {m n : Type*} [Fintype n] (A B : Matrix m n F) :
    (A + B).rank ≤ A.rank + B.rank := by
  rw [← finrank_range_mulVecLin, ← finrank_range_mulVecLin, ← finrank_range_mulVecLin,
    Matrix.mulVecLin_add]
  refine le_trans (Submodule.finrank_mono ?_) (Submodule.finrank_add_le_finrank_add_finrank _ _)
  rintro _ ⟨x, rfl⟩
  exact Submodule.add_mem_sup (LinearMap.mem_range_self _ x) (LinearMap.mem_range_self _ x)

theorem rank_sum_le {m n : Type*} [Fintype n] {ι : Type*} (s : Finset ι) (M : ι → Matrix m n F) :
    (∑ k ∈ s, M k).rank ≤ ∑ k ∈ s, (M k).rank := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (rank_add_le _ _).trans (Nat.add_le_add_left ih _)

/-- A matrix of the form `M y x = Fy y (ψ x) * Hx (φ y) x` factors through `P × Q`. -/
theorem rank_le_of_factor {R C P Q : Type*} [Fintype R] [Fintype C] [Fintype P] [Fintype Q]
    (φ : R → P) (ψ : C → Q) (Fy : R → Q → F) (Hx : P → C → F)
    (M : Matrix R C F) (hM : ∀ y x, M y x = Fy y (ψ x) * Hx (φ y) x) :
    M.rank ≤ Fintype.card P * Fintype.card Q := by
  classical
  let U : Matrix R (P × Q) F := fun y pq => if φ y = pq.1 then Fy y pq.2 else 0
  let V : Matrix (P × Q) C F := fun pq x => if ψ x = pq.2 then Hx pq.1 x else 0
  have h : M = U * V := by
    ext y x
    rw [hM, Matrix.mul_apply, Fintype.sum_prod_type]
    rw [Finset.sum_eq_single (φ y), Finset.sum_eq_single (ψ x)]
    · simp [U, V]
    · intro q _ hq
      simp [U, V, Ne.symm hq]
    · simp
    · intro p _ hp
      simp [U, Ne.symm hp]
    · simp
  rw [h]
  calc (U * V).rank ≤ U.rank := Matrix.rank_mul_le_left _ _
    _ ≤ Fintype.card (P × Q) := Matrix.rank_le_card_width _
    _ = Fintype.card P * Fintype.card Q := Fintype.card_prod _ _

/-- Rank factorisation: a matrix of rank `ρ` is a sum of `ρ` products `U i t * V t c`. -/
theorem exists_factor {R C : Type*} [Fintype R] [Fintype C] (M : Matrix R C F) :
    ∃ (ρ : ℕ) (U : R → Fin ρ → F) (V : Fin ρ → C → F), ρ = M.rank ∧
      ∀ i c, M i c = ∑ t, U i t * V t c := by
  classical
  let S := Submodule.span F (Set.range M.col)
  let bs := Module.finBasis F S
  refine ⟨finrank F S, fun i t => (bs t : R → F) i,
    fun t c => bs.repr ⟨M.col c, Submodule.subset_span ⟨c, rfl⟩⟩ t, rfl, ?_⟩
  intro i c
  have h := bs.sum_repr ⟨M.col c, Submodule.subset_span ⟨c, rfl⟩⟩
  have h2 := congrArg (fun v : S => (v : R → F) i) h
  simp only [Submodule.coe_sum, Submodule.coe_smul, Finset.sum_apply, Pi.smul_apply,
    smul_eq_mul] at h2
  have h3 : M.col c i = M i c := rfl
  rw [← h3, ← h2]
  simp only [mul_comm]

/-- `X ⊗ 1_Q`, indexed as `(a × Fin Q) × (b × Fin Q)`. -/
def kronId {a b : Type*} (X : Matrix a b F) (Q : ℕ) : Matrix (a × Fin Q) (b × Fin Q) F :=
  Matrix.of fun x y => if x.2 = y.2 then X x.1 y.1 else 0

/-- `rank (X ⊗ 1_Q) ≥ Q rank X`. -/
theorem rank_mul_le_rank_kronId {a b : Type*} [Fintype a] [Fintype b] [DecidableEq a]
    (X : Matrix a b F) (Q : ℕ) : X.rank * Q ≤ (kronId X Q).rank := by
  classical
  let S := Submodule.span F (Set.range X.col)
  let bs := Module.finBasis F S
  let L : Fin Q → (a → F) →ₗ[F] (a × Fin Q → F) := fun w0 =>
    { toFun := fun v x => if x.2 = w0 then v x.1 else 0
      map_add' := by intro u v; funext x; by_cases hx : x.2 = w0 <;> simp [hx]
      map_smul' := by intro c v; funext x; by_cases hx : x.2 = w0 <;> simp [hx] }
  let fam : Fin (finrank F S) × Fin Q → (a × Fin Q → F) := fun tw => L tw.2 (bs tw.1)
  have hLI : LinearIndependent F fam := by
    rw [Fintype.linearIndependent_iff]
    intro c hc tw
    obtain ⟨t, w⟩ := tw
    have hw : ∑ t, c (t, w) • bs t = 0 := by
      apply Subtype.ext
      funext i
      have := congrFun hc (i, w)
      rw [Finset.sum_apply] at this
      simp only [Pi.smul_apply, smul_eq_mul, Pi.zero_apply, fam, L, LinearMap.coe_mk,
        AddHom.coe_mk, Fintype.sum_prod_type] at this
      rw [Finset.sum_comm, Finset.sum_eq_single w] at this
      · simpa using this
      · intro w2 _ hw2
        simp [Ne.symm hw2]
      · simp
    exact (Fintype.linearIndependent_iff.1 bs.linearIndependent) _ hw t
  have hspan : Submodule.span F (Set.range fam) ≤ Submodule.span F (Set.range (kronId X Q).col) := by
    rw [Submodule.span_le]
    rintro _ ⟨⟨t, w0⟩, rfl⟩
    have h1 : (bs t : a → F) ∈ S := (bs t).2
    have h2 : L w0 (bs t : a → F) ∈ S.map (L w0) := Submodule.mem_map_of_mem h1
    rw [Submodule.map_span] at h2
    refine Submodule.span_mono ?_ h2
    rintro _ ⟨_, ⟨j, rfl⟩, rfl⟩
    refine ⟨(j, w0), ?_⟩
    funext x
    simp [kronId, L, Matrix.col, eq_comm]
  have := Submodule.finrank_mono hspan
  rw [finrank_span_eq_card hLI, Fintype.card_prod, Fintype.card_fin, Fintype.card_fin] at this
  exact this

end Rank

section Swap

variable (K : Type*) [Field K]

/-- swap with a spectator, `E ⊗ G ⊗ W → G ⊗ E ⊗ W` (`dim E = m`, `dim G = n`, `dim W = Q`);
rows = outputs `(g, e, w)`, columns = inputs `(e', g', w')`. -/
def swapW (m n Q : ℕ) : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K :=
  fun o i => if o.1 = i.2.1 ∧ o.2.1 = i.1 ∧ o.2.2 = i.2.2 then 1 else 0

/-- realignment across site 0 | rest: rows `(out₀, in₀) = (g, e')`,
columns `(out_rest, in_rest) = ((e, w), (g', w'))`. -/
def realign {m n Q : ℕ} (T : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K) :
    Matrix (Fin n × Fin m) ((Fin m × Fin Q) × (Fin n × Fin Q)) K :=
  fun a b => T (a.1, b.1.1, b.1.2) (a.2, b.2.1, b.2.2)

/-- `(g, j, w) ↦ (j, g, w)`. -/
def swap12 (a b c : Type*) : a × b × c ≃ b × a × c where
  toFun x := (x.2.1, x.1, x.2.2)
  invFun x := (x.2.1, x.1, x.2.2)
  left_inv _ := rfl
  right_inv _ := rfl

/-- **Theorem C (a)**: `r ≥ k Q (n - k b)` for `k ≤ m`, in the form `k Q n ≤ r + k² Q b`. -/
theorem part_a (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (k : ℕ) (hk : k ≤ m) : k * Q * n ≤ r + k * k * Q * b := by
  classical
  let f : Fin k → Fin m := Fin.castLE hk
  let rw' : Fin n × Fin k × Fin Q → Fin n × Fin m × Fin Q := fun x => (x.1, f x.2.1, x.2.2)
  let cl : Fin k × Fin n × Fin Q → Fin m × Fin n × Fin Q := fun x => (f x.1, x.2.1, x.2.2)
  -- the compressed swap is a permutation matrix
  have hS : (swapW K m n Q).submatrix rw' cl =
      (1 : Matrix (Fin k × Fin n × Fin Q) (Fin k × Fin n × Fin Q) K).submatrix
        (swap12 _ _ _) (Equiv.refl _) := by
    ext ⟨g, j, w⟩ ⟨j', g', w'⟩
    simp only [swapW, submatrix_apply, Matrix.one_apply, swap12, Equiv.coe_fn_mk, Equiv.refl_apply, rw', cl,
      f, Fin.castLE_inj, Prod.mk.injEq]
    congr 1
    apply propext
    tauto
  have hSr : ((swapW K m n Q).submatrix rw' cl).rank = k * n * Q := by
    rw [hS, rank_submatrix (1 : Matrix (Fin k × Fin n × Fin Q) (Fin k × Fin n × Fin Q) K)
      (swap12 _ _ _) (Equiv.refl _)]
    simp [Matrix.rank_one, Fintype.card_prod, mul_assoc]
  -- the compressed `B` has rank `≤ ρ k Q k`
  obtain ⟨ρ, U, V, hρ, hUV⟩ := exists_factor (realign K B)
  let T : Fin ρ → Matrix (Fin n × Fin k × Fin Q) (Fin k × Fin n × Fin Q) K :=
    fun t y x => U (y.1, f x.1) t * V t ((f y.2.1, y.2.2), (x.2.1, x.2.2))
  have hBT : B.submatrix rw' cl = ∑ t, T t := by
    ext y x
    rw [Matrix.sum_apply]
    exact hUV (y.1, f x.1) ((f y.2.1, y.2.2), (x.2.1, x.2.2))
  have hT : ∀ t, (T t).rank ≤ k * Q * k := by
    intro t
    have := rank_le_of_factor (fun y : Fin n × Fin k × Fin Q => (y.2.1, y.2.2))
      (fun x : Fin k × Fin n × Fin Q => x.1) (fun y j' => U (y.1, f j') t)
      (fun p x => V t ((f p.1, p.2), (x.2.1, x.2.2))) (T t) (fun _ _ => rfl)
    simpa [Fintype.card_prod] using this
  have hBr : (B.submatrix rw' cl).rank ≤ ρ * (k * Q * k) := by
    rw [hBT]
    refine (rank_sum_le _ _).trans ?_
    calc ∑ t, (T t).rank ≤ ∑ _t : Fin ρ, k * Q * k := Finset.sum_le_sum fun t _ => hT t
      _ = ρ * (k * Q * k) := by simp
  have hAr : (A.submatrix rw' cl).rank ≤ r := (rank_submatrix_le A rw' cl).trans hA
  have key : k * n * Q ≤ r + ρ * (k * Q * k) := by
    rw [← hSr, ← h, Matrix.submatrix_add]
    exact (rank_add_le _ _).trans (Nat.add_le_add hAr hBr)
  have hρb : ρ * (k * Q * k) ≤ b * (k * Q * k) := Nat.mul_le_mul_right _ (hρ ▸ hB)
  calc k * Q * n = k * n * Q := by ring
    _ ≤ r + b * (k * Q * k) := key.trans (Nat.add_le_add_left hρb _)
    _ = r + k * k * Q * b := by ring

/-- **Theorem C (b)**: `b Q ≥ k (n Q - k r)` for `k ≤ m`, in the form `k n Q ≤ b Q + k² r`. -/
theorem part_b (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (k : ℕ) (hk : k ≤ m) : k * n * Q ≤ b * Q + k * k * r := by
  classical
  let f : Fin k → Fin m := Fin.castLE hk
  -- `ψ ↦ (A (e_{f j} ⊗ ψ))_j`
  let At : Matrix ((Fin n × Fin m × Fin Q) × Fin k) (Fin n × Fin Q) K :=
    fun x c => A x.1 (f x.2, c.1, c.2)
  let D : Fin k → Matrix ((Fin n × Fin m × Fin Q) × Fin k) (Fin n × Fin Q) K := fun j =>
    (diagonal (fun x : (Fin n × Fin m × Fin Q) × Fin k => if x.2 = j then (1 : K) else 0)) *
      A.submatrix Prod.fst (fun c : Fin n × Fin Q => (f j, c.1, c.2))
  have hAt : At.rank ≤ k * r := by
    have hsum : At = ∑ j, D j := by
      ext x c
      simp [Matrix.sum_apply, diagonal_mul, At, D]
    have hD : ∀ j, (D j).rank ≤ r := fun j =>
      (rank_mul_le_right _ _).trans ((rank_submatrix_le _ _ _).trans hA)
    rw [hsum]
    refine (rank_sum_le _ _).trans ?_
    calc ∑ j, (D j).rank ≤ ∑ _j : Fin k, r := Finset.sum_le_sum fun j _ => hD j
      _ = k * r := by simp
  -- `Ψ = ker At`, with a basis `bs`
  let S := LinearMap.ker At.mulVecLin
  let bs := Module.finBasis K S
  set d := finrank K S with hd_def
  have hd : At.rank + d = n * Q := by
    rw [← finrank_range_mulVecLin]
    have := LinearMap.finrank_range_add_finrank_ker At.mulVecLin
    simpa [Module.finrank_fintype_fun_eq_card, Fintype.card_prod] using this
  let Z : Fin n × Fin Q → Fin d → K := fun c s => (bs s : Fin n × Fin Q → K) c
  have hZ : ∀ x s, ∑ c, At x c * Z c s = 0 := by
    intro x s
    have h2 : (At *ᵥ (bs s : Fin n × Fin Q → K)) x = 0 := by
      have := LinearMap.mem_ker.mp (bs s).2
      rw [Matrix.mulVecLin_apply] at this
      rw [this]; rfl
    simp only [Matrix.mulVec, dotProduct] at h2
    exact h2
  -- `M = 1_k ⊗ Z` has rank `k d`
  let M : Matrix (Fin k × Fin n × Fin Q) (Fin k × Fin d) K :=
    fun y x => if y.1 = x.1 then Z (y.2.1, y.2.2) x.2 else 0
  have hMr : M.rank = k * d := by
    rw [← rank_transpose]
    have hli : LinearIndependent K Mᵀ.row := by
      rw [Fintype.linearIndependent_iff]
      intro c hc js
      obtain ⟨j, s⟩ := js
      have hj : ∑ s, c (j, s) • bs s = 0 := by
        apply Subtype.ext
        funext gw
        obtain ⟨g, w⟩ := gw
        have := congrFun hc (j, g, w)
        rw [Finset.sum_apply] at this
        simp only [Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at this
        have e : ∀ x : Fin k × Fin d, Mᵀ.row x (j, g, w) = if j = x.1 then Z (g, w) x.2 else 0 :=
          fun _ => rfl
        simp only [e, mul_ite, mul_zero, Fintype.sum_prod_type] at this
        rw [Finset.sum_eq_single j] at this
        · simpa [Z] using this
        · intro j2 _ hj2
          simp [Ne.symm hj2]
        · simp
      exact (Fintype.linearIndependent_iff.1 bs.linearIndependent) _ hj s
    rw [hli.rank_matrix]
    simp
  -- factorisation of `B`
  obtain ⟨ρ, U, V, hρ, hUV⟩ := exists_factor (realign K B)
  let V' : Fin ρ → (Fin m × Fin Q) × Fin d → K :=
    fun t x => ∑ c : Fin n × Fin Q, V t (x.1, c) * Z c x.2
  let T : Fin ρ → Matrix (Fin k × Fin n × Fin Q) (Fin k × Fin d) K :=
    fun t y x => U (y.2.1, f y.1) t * V' t ((f x.1, y.2.2), x.2)
  have hMT : M = ∑ t, T t := by
    ext ⟨j, g, w⟩ ⟨j2, s⟩
    have e1 : M (j, g, w) (j2, s) =
        ∑ c : Fin n × Fin Q, swapW K m n Q (g, f j2, w) (f j, c.1, c.2) * Z c s := by
      rw [Finset.sum_eq_single (g, w)]
      · by_cases hj : j = j2
        · subst hj; simp [M, swapW]
        · have : f j2 ≠ f j := fun he => hj (Fin.castLE_injective hk he).symm
          simp [M, swapW, hj, this]
      · intro c _ hc
        have : ¬ (g = c.1 ∧ f j2 = f j ∧ w = c.2) := fun hh => hc (by
          obtain ⟨c1, c2⟩ := c
          simp only at hh
          rw [hh.1, hh.2.2])
        simp [swapW, this]
      · simp
    rw [e1, ← h]
    simp only [Matrix.add_apply, add_mul, Finset.sum_add_distrib]
    have hz := hZ ((g, f j2, w), j) s
    simp only [At] at hz
    rw [hz, zero_add, Matrix.sum_apply]
    simp only [T, V']
    have hB' : ∀ c : Fin n × Fin Q, B (g, f j2, w) (f j, c.1, c.2) =
        ∑ t, U (g, f j) t * V t ((f j2, w), c) := fun c => hUV (g, f j) ((f j2, w), c)
    simp only [hB', Finset.sum_mul, Finset.mul_sum, mul_assoc]
    exact Finset.sum_comm
  have hT : ∀ t, (T t).rank ≤ Q * 1 := by
    intro t
    have := rank_le_of_factor (fun y : Fin k × Fin n × Fin Q => y.2.2)
      (fun _ : Fin k × Fin d => ()) (fun y _ => U (y.2.1, f y.1) t)
      (fun w x => V' t ((f x.1, w), x.2)) (T t) (fun _ _ => rfl)
    simpa using this
  have hMB : k * d ≤ ρ * Q := by
    rw [← hMr, hMT]
    refine (rank_sum_le _ _).trans ?_
    calc ∑ t, (T t).rank ≤ ∑ _t : Fin ρ, Q * 1 := Finset.sum_le_sum fun t _ => hT t
      _ = ρ * Q := by simp
  have hρb : ρ * Q ≤ b * Q := Nat.mul_le_mul_right _ (hρ ▸ hB)
  calc k * n * Q = k * (n * Q) := by ring
    _ = k * At.rank + k * d := by rw [← hd]; ring
    _ ≤ k * (k * r) + b * Q := Nat.add_le_add (Nat.mul_le_mul_left _ hAt) (hMB.trans hρb)
    _ = b * Q + k * k * r := by ring

/-! ### Mirror statements, by transposition -/

theorem swapW_transpose (m n Q : ℕ) : (swapW K m n Q)ᵀ = swapW K n m Q := by
  ext x y
  simp only [transpose_apply, swapW]
  congr 1
  apply propext
  constructor <;> rintro ⟨h1, h2, h3⟩ <;> exact ⟨h2.symm, h1.symm, h3.symm⟩

theorem rank_realign_transpose {m n Q : ℕ}
    (B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K) :
    (realign K Bᵀ).rank = (realign K B).rank := by
  have : realign K Bᵀ = (realign K B).submatrix (Equiv.prodComm _ _) (Equiv.prodComm _ _) := by
    ext a c
    rfl
  rw [this, rank_submatrix]

/-- **Theorem C (a), mirror**: `r ≥ k Q (m - k b)` for `k ≤ n`. -/
theorem part_a' (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (k : ℕ) (hk : k ≤ n) : k * Q * m ≤ r + k * k * Q * b := by
  refine part_a K n m Q r b Aᵀ Bᵀ ?_ (by rwa [rank_transpose]) (by rwa [rank_realign_transpose])
    k hk
  rw [← transpose_add, h, swapW_transpose]

/-- **Theorem C (b), mirror**: `b Q ≥ k (m Q - k r)` for `k ≤ n`. -/
theorem part_b' (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (k : ℕ) (hk : k ≤ n) : k * m * Q ≤ b * Q + k * k * r := by
  refine part_b K n m Q r b Aᵀ Bᵀ ?_ (by rwa [rank_transpose]) (by rwa [rank_realign_transpose])
    k hk
  rw [← transpose_add, h, swapW_transpose]

/-- **Theorem C (c)** (exact case): if `osr(B) < m n` then `r ≥ Q`. -/
theorem part_c (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (hb : b < m * n) : Q ≤ r := by
  classical
  -- a nonzero functional `c` on `Hom(E,G)` killing the realigned `B`
  obtain ⟨c, hc, hc0⟩ : ∃ c : Fin n × Fin m → K,
      c ∈ LinearMap.ker (realign K B)ᵀ.mulVecLin ∧ c ≠ 0 := by
    apply Submodule.exists_mem_ne_zero_of_ne_bot
    intro hbot
    have h1 := LinearMap.finrank_range_add_finrank_ker (realign K B)ᵀ.mulVecLin
    rw [hbot, finrank_bot, finrank_range_mulVecLin, rank_transpose,
      Module.finrank_fintype_fun_eq_card, Fintype.card_prod, Fintype.card_fin,
      Fintype.card_fin, Nat.mul_comm n m] at h1
    omega
  have hcB : ∀ y, ∑ x, c x * realign K B x y = 0 := by
    intro y
    have h2 : ((realign K B)ᵀ *ᵥ c) y = 0 := by
      have := LinearMap.mem_ker.mp hc
      rw [Matrix.mulVecLin_apply] at this
      rw [this]; rfl
    simp only [Matrix.mulVec, dotProduct, transpose_apply] at h2
    rw [← h2]
    exact Finset.sum_congr rfl fun x _ => mul_comm _ _
  -- contraction of site 0 with `c`
  let Φ : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K →
      Matrix (Fin m × Fin Q) (Fin n × Fin Q) K := fun X =>
    Matrix.of fun y z => ∑ x : Fin n × Fin m, c x * X (x.1, y.1, y.2) (x.2, z.1, z.2)
  have hΦadd : Φ (A + B) = Φ A + Φ B := by
    ext y z
    simp [Φ, mul_add, Finset.sum_add_distrib]
  have hΦB : Φ B = 0 := by
    ext y z
    exact hcB ((y.1, y.2), (z.1, z.2))
  let C : Matrix (Fin n) (Fin m) K := Matrix.of fun g e => c (g, e)
  have hΦP : Φ (swapW K m n Q) = kronId Cᵀ Q := by
    ext y z
    simp only [Φ, swapW, kronId, Matrix.of_apply, transpose_apply, C, mul_ite, mul_one, mul_zero]
    rw [Finset.sum_eq_single (z.1, y.1)]
    · by_cases hw : y.2 = z.2 <;> simp [hw]
    · intro x _ hx
      have : ¬ (x.1 = z.1 ∧ y.1 = x.2 ∧ y.2 = z.2) := fun hh => hx (by
        obtain ⟨x1, x2⟩ := x
        simp only at hh
        rw [hh.1, hh.2.1])
      simp [this]
    · simp
  -- `rank (Φ A) ≤ ρ r` with `ρ = rank C`
  obtain ⟨ρ, U, V, hρ, hUV⟩ := exists_factor C
  let Lt : Fin ρ → Matrix (Fin m × Fin Q) (Fin n × Fin m × Fin Q) K := fun t =>
    Matrix.of fun y o => if o.2.1 = y.1 ∧ o.2.2 = y.2 then U o.1 t else 0
  let Rt : Fin ρ → Matrix (Fin m × Fin n × Fin Q) (Fin n × Fin Q) K := fun t =>
    Matrix.of fun i z => if i.2.1 = z.1 ∧ i.2.2 = z.2 then V t i.1 else 0
  have hΦA : Φ A = ∑ t, Lt t * A * Rt t := by
    ext y z
    rw [Matrix.sum_apply]
    simp only [Φ, Matrix.of_apply]
    have hc' : ∀ x : Fin n × Fin m, c x = ∑ t, U x.1 t * V t x.2 := fun x => hUV x.1 x.2
    simp only [hc', Finset.sum_mul]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun t _ => ?_
    simp only [Matrix.mul_apply, Lt, Rt, Matrix.of_apply, Fintype.sum_prod_type, ite_mul,
      zero_mul, mul_ite, mul_zero]
    simp only [Finset.sum_mul]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun g _ => ?_
    -- remaining: collapse the delta sums
    simp [Finset.sum_ite_eq', ite_and, mul_comm, mul_left_comm, mul_assoc]
  have hLt : ∀ t, (Lt t * A * Rt t).rank ≤ r := fun t =>
    (rank_mul_le_left _ _).trans ((rank_mul_le_right _ _).trans hA)
  have hAr : (Φ A).rank ≤ ρ * r := by
    rw [hΦA]
    refine (rank_sum_le _ _).trans ?_
    calc ∑ t, (Lt t * A * Rt t).rank ≤ ∑ _t : Fin ρ, r := Finset.sum_le_sum fun t _ => hLt t
      _ = ρ * r := by simp
  have hρpos : 0 < ρ := by
    rcases Nat.eq_zero_or_pos ρ with h0 | h0
    · exfalso
      apply hc0
      funext x
      have := hUV x.1 x.2
      subst h0
      simpa [C] using this
    · exact h0
  have hkr : ρ * Q ≤ (Φ (swapW K m n Q)).rank := by
    rw [hΦP]
    have := rank_mul_le_rank_kronId Cᵀ Q
    rwa [rank_transpose, ← hρ] at this
  rw [← h, hΦadd, hΦB, add_zero] at hkr
  exact Nat.le_of_mul_le_mul_left (hkr.trans hAr) hρpos

/-- **Theorem C**, all parts, in the form of `g19e/notes.md` §6. -/
theorem swap_rank_osr_tradeoff (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b) :
    (∀ k ≤ m, k * Q * n ≤ r + k * k * Q * b ∧ k * n * Q ≤ b * Q + k * k * r) ∧
    (∀ k ≤ n, k * Q * m ≤ r + k * k * Q * b ∧ k * m * Q ≤ b * Q + k * k * r) ∧
    (b < m * n → Q ≤ r) :=
  ⟨fun k hk => ⟨part_a K m n Q r b A B h hA hB k hk, part_b K m n Q r b A B h hA hB k hk⟩,
   fun k hk => ⟨part_a' K m n Q r b A B h hA hB k hk, part_b' K m n Q r b A B h hA hB k hk⟩,
   part_c K m n Q r b A B h hA hB⟩

end Swap

/-! ### Arithmetic corollary -/

theorem exists_ceil_mul (M d : ℕ) (hd : 0 < d) : ∃ k, M ≤ d * k ∧ d * k < M + d := by
  refine ⟨(M + d - 1) / d, ?_, ?_⟩
  · have h1 := Nat.div_add_mod (M + d - 1) d
    have h2 := Nat.mod_lt (M + d - 1) hd
    omega
  · have h1 := Nat.div_mul_le_self (M + d - 1) d
    rw [Nat.mul_comm]
    omega

theorem arith_cor (M Q r b : ℕ) (hb1 : 1 ≤ b) (hQr : Q ≤ r)
    (ha : ∀ k ≤ M, k * Q * M ≤ r + k * k * Q * b)
    (hb : ∀ k ≤ M, k * M * Q ≤ b * Q + k * k * r) : 2 * Q * M ^ 2 ≤ 9 * (r * b) := by
  rcases Nat.eq_zero_or_pos Q with hQ | hQ
  · simp [hQ]
  rcases Nat.eq_zero_or_pos M with hM | hM
  · simp [hM]
  by_cases h1 : b * Q ≤ r
  · by_cases h2 : 3 * b ≤ M
    · obtain ⟨k, hk1, hk2⟩ := exists_ceil_mul M (3 * b) (by omega)
      have hkM : k ≤ M := by nlinarith
      have hak := ha k hkM
      have hu2 : 3 * b * k ≤ 2 * M := by omega
      zify at *
      nlinarith [mul_nonneg (mul_nonneg (sub_nonneg.2 hk1) (sub_nonneg.2 hu2)) hQ.le]
    · have ha1 := ha 1 (by omega)
      by_cases h3 : 3 * b ≤ 2 * M
      · zify at *
        nlinarith [mul_nonneg (mul_nonneg (sub_nonneg.2 (le_of_lt (not_le.1 h2)))
          (sub_nonneg.2 h3)) hQ.le]
      · zify at *
        nlinarith
  · push Not at h1
    have hr : 0 < r := lt_of_lt_of_le hQ hQr
    suffices hs : Q * (2 * Q * M ^ 2) ≤ Q * (9 * (r * b)) from Nat.le_of_mul_le_mul_left hs hQ
    by_cases h2 : 3 * r ≤ M * Q
    · obtain ⟨k, hk1, hk2⟩ := exists_ceil_mul (M * Q) (3 * r) (by omega)
      have hkM : k ≤ M := by
        by_contra hc
        push Not at hc
        nlinarith
      have hbk := hb k hkM
      have hu2 : 3 * r * k ≤ 2 * (M * Q) := by omega
      zify at *
      nlinarith [mul_nonneg (sub_nonneg.2 hk1) (sub_nonneg.2 hu2)]
    · have hb1' := hb 1 (by omega)
      by_cases h3 : 3 * r ≤ 2 * (M * Q)
      · zify at *
        nlinarith [mul_nonneg (sub_nonneg.2 (le_of_lt (not_le.1 h2))) (sub_nonneg.2 h3)]
      · have h3' : 2 * (M * Q) < 3 * r := not_le.1 h3
        have e1 : r * r < r * (b * Q) := Nat.mul_lt_mul_of_pos_left h1 hr
        have e2 : (2 * (M * Q)) * (2 * (M * Q)) < (3 * r) * (3 * r) :=
          Nat.mul_self_lt_mul_self h3'
        zify at *
        nlinarith

/-- arithmetic corollary (`m = n = M`): `r b ≥ (2/9) Q M²` for `1 ≤ b < M²`. -/
theorem swap_rank_osr_product (K : Type*) [Field K] (M Q r b : ℕ)
    (A B : Matrix (Fin M × Fin M × Fin Q) (Fin M × Fin M × Fin Q) K)
    (h : A + B = swapW K M M Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (hb : 1 ≤ b) (hb' : b < M * M) : 2 * Q * M ^ 2 ≤ 9 * (r * b) :=
  arith_cor M Q r b hb (part_c K M M Q r b A B h hA hB hb')
    (fun k hk => part_a K M M Q r b A B h hA hB k hk)
    (fun k hk => part_b K M M Q r b A B h hA hB k hk)

end SwapRank

#print axioms SwapRank.part_a
#print axioms SwapRank.part_b
#print axioms SwapRank.part_a'
#print axioms SwapRank.part_b'
#print axioms SwapRank.part_c
#print axioms SwapRank.swap_rank_osr_tradeoff
#print axioms SwapRank.swap_rank_osr_product
