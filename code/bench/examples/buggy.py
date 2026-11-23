"""Deliberately buggy candidates for testing the harness."""
def say_unsat_if_dense(n, clauses):
    # wrong heuristic: claims UNSAT when density > 4.5, else greedy SAT (may be wrong)
    if len(clauses) > 4.5 * n:
        return ("UNSAT", None)
    a = {v: True for v in range(1, n + 1)}
    return ("SAT", a)  # often a bad assignment
