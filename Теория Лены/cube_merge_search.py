#!/usr/bin/env python3
"""
Exploration of the "merging rule" from Lena's Theory.

Physical framing (from the issue): elementary lattice points may merge into a
common volume only when the equation

    a_1^3 + a_2^3 + ... + a_k^3 = z^3 ,   k in {3, 4, 5, 6, 7, 8}

holds (a "merge" of k smaller cubes into one cube of side z). The reverse is a
"separation". We search for the smallest such mergers for every k, count how
the number of distinct mergers grows with z (the claim "the larger the number,
the more variations"), and check the special status of k = 8.

All output here is mechanically verified, not asserted.
"""

from collections import defaultdict
from itertools import combinations_with_replacement


def representations(z, k, allow_equal=True):
    """All ways to write z^3 as a sum of exactly k positive cubes a_i (a_i <= z).
    Sorted tuples, a_1 <= a_2 <= ... <= a_k, so each multiset counted once."""
    target = z ** 3
    results = []
    bound = z  # any part >= z would already reach/exceed z^3 with k>=2

    def rec(remaining, parts_left, max_val, acc):
        if parts_left == 0:
            if remaining == 0:
                results.append(tuple(acc))
            return
        if remaining <= 0:
            return
        # prune: even using max_val for all remaining parts must reach target
        if max_val ** 3 * parts_left < remaining:
            return
        # prune: smallest usable part is 1; remaining must be >= parts_left
        for v in range(min(max_val, bound), 0, -1):
            c = v ** 3
            if c * parts_left < remaining:
                break
            if c > remaining:
                continue
            acc.append(v)
            rec(remaining - c, parts_left - 1, v, acc)
            acc.pop()

    rec(target, k, bound - 1, [])
    return [tuple(reversed(r)) for r in results]


def first_merger(k, zmax=120):
    """Smallest z (and one representation) such that z^3 is a sum of k cubes."""
    for z in range(2, zmax + 1):
        reps = representations(z, k)
        if reps:
            return z, reps[0], len(reps)
    return None


def growth_table(k, zmax=80):
    """Count distinct mergers (representations) of z^3 as k cubes, for each z."""
    return [(z, len(representations(z, k))) for z in range(2, zmax + 1)]


if __name__ == "__main__":
    print("=" * 70)
    print("SMALLEST MERGER FOR EACH NUMBER OF TERMS k")
    print("=" * 70)
    for k in range(2, 9):
        res = first_merger(k)
        if res is None:
            print(f"k={k}: no solution found up to bound "
                  f"(expected for k=2 by Fermat's Last Theorem)")
            continue
        z, rep, count = res
        lhs = " + ".join(f"{a}^3" for a in rep)
        print(f"k={k}: {lhs} = {z}^3   "
              f"(={z**3}; {count} distinct merger(s) at z={z})")

    print()
    print("=" * 70)
    print("THE k=8 CLAIM:  z^3 as sum of EIGHT cubes, smallest cases")
    print("=" * 70)
    for z in range(2, 16):
        reps = representations(z, 8)
        if reps:
            for r in reps[:3]:
                lhs = " + ".join(f"{a}^3" for a in r)
                print(f"  {lhs} = {z}^3")
            if len(reps) > 3:
                print(f"  ... and {len(reps) - 3} more at z={z}")
            print(f"  -> first z with an 8-cube merger: {z} "
                  f"({len(reps)} total)\n")
            break

    print("=" * 70)
    print("GROWTH OF VARIATIONS: number of distinct mergers vs z")
    print("(demonstrates 'the larger the number, the more variations')")
    print("=" * 70)
    header = "  z  | " + " | ".join(f"k={k}" for k in range(3, 9))
    print(header)
    print("  " + "-" * (len(header) - 2))
    totals = defaultdict(int)
    for z in range(2, 41):
        counts = [len(representations(z, k)) for k in range(3, 9)]
        for k, c in zip(range(3, 9), counts):
            totals[k] += c
        if any(counts):
            row = f"  {z:>2} | " + " | ".join(f"{c:>3}" for c in counts)
            print(row)
    print("  " + "-" * (len(header) - 2))
    print("  cumulative (z<=40): " +
          " | ".join(f"{totals[k]:>3}" for k in range(3, 9)))
