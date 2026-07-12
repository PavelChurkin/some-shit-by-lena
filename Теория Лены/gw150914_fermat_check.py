#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Проверка «салфетки» к слиянию GW150914 на основе ОУФ (issue #12).

Скрипт механически проверяет три группы утверждений из присланного расчёта:

  A. Физическая бухгалтерия: перевод масс чёрных дыр в «число узлов» N = M / m_P
     и наблюдаемая доля потерь энергии для GW150914.
  B. Ядро гипотезы: спектр «стабильных узлов» N = 3a + 8b (a, b >= 0) — это
     числовая полугруппа <3, 8>. Считаем её лакуны (недостижимые N) и
     асимптотическую плотность «запрещённых» состояний.
  C. Заявленный «фундаментальный предел потерь» 1/НОК(3,8) = 1/24 и сравнение
     его с реальной плотностью запрещённых состояний.

Вывод скрипта — это факты; их интерпретация приведена в отчёте рядом.
Без внешних зависимостей: python3 gw150914_fermat_check.py
"""

from math import gcd


def lcm(a, b):
    return a * b // gcd(a, b)


# --- A. Физическая бухгалтерия -------------------------------------------------

M_SUN = 1.98892e30       # кг, масса Солнца
M_PLANCK = 2.176434e-8   # кг, планковская масса (CODATA)

M1, M2, MF = 36.0, 29.0, 62.0  # массы в M_sun (LIGO, округлённо)


def solar_to_nodes(m_solar):
    """Число «узлов» N = M / m_P для массы в солнечных единицах."""
    return m_solar * M_SUN / M_PLANCK


def physics_bookkeeping():
    N1, N2, NF = (solar_to_nodes(m) for m in (M1, M2, MF))
    N_init = N1 + N2
    dN = N_init - NF
    frac = dN / N_init
    return {
        "N1": N1, "N2": N2, "NF": NF,
        "N_init": N_init, "dN": dN, "loss_fraction": frac,
        "radiated_solar": (M1 + M2) - MF,
    }


# --- B. Числовая полугруппа <3, 8> --------------------------------------------

def semigroup_gaps(gens, bound=1000):
    """Недостижимые целые (лакуны) для полугруппы, порождённой gens."""
    reachable = set()
    # достаточно перебрать коэффициенты до bound // min(gens)
    span = bound // min(gens) + 1
    a_gen, b_gen = gens
    for a in range(span + 1):
        base = a_gen * a
        if base > bound:
            break
        b = 0
        while base + b_gen * b <= bound:
            reachable.add(base + b_gen * b)
            b += 1
    return [n for n in range(bound + 1) if n not in reachable]


def forbidden_density(gaps, N):
    """Доля запрещённых состояний среди 1..N (все лакуны конечны и малы)."""
    return sum(1 for g in gaps if 1 <= g <= N) / N


# --- C. Заявленный предел ------------------------------------------------------

def main():
    print("=" * 70)
    print("A. Физическая бухгалтерия GW150914")
    print("=" * 70)
    p = physics_bookkeeping()
    for m in (M1, M2, MF):
        print(f"  {m:>4.0f} M_sun -> N = {solar_to_nodes(m):.3e}")
    print(f"  N_нач = N1 + N2      = {p['N_init']:.3e}")
    print(f"  N_кон (M_f)         = {p['NF']:.3e}")
    print(f"  ΔN                  = {p['dN']:.3e}")
    print(f"  Излучено массы       = {p['radiated_solar']:.1f} M_sun")
    print(f"  Доля потерь ΔN/N_нач = {p['loss_fraction']*100:.3f}%  "
          f"(в «салфетке»: ~4.7%)")

    print()
    print("=" * 70)
    print("B. Спектр N = 3a + 8b  —  числовая полугруппа <3, 8>")
    print("=" * 70)
    gaps = semigroup_gaps((3, 8), bound=1000)
    genus_formula = (3 - 1) * (8 - 1) // 2  # для двух взаимно простых
    frobenius = 3 * 8 - 3 - 8               # число Фробениуса
    print(f"  Лакуны (недостижимые N): {gaps}")
    print(f"  Число лакун (род)       : {len(gaps)}  "
          f"(формула (3-1)(8-1)/2 = {genus_formula})")
    print(f"  Число Фробениуса        : {max(gaps)}  "
          f"(формула 3·8-3-8 = {frobenius})")
    print("  => все N > 13 представимы; запрещённых больших состояний НЕТ.")
    print()
    print("  Плотность запрещённых состояний среди 1..N:")
    for N in (50, 100, 1000, 10**6, 10**39):
        # для N=10^39 плотность = (число лакун)/N, лакуны все <= 13
        dens = len(gaps) / N if N > max(gaps) else forbidden_density(gaps, N)
        print(f"    N <= {N:<8}: {dens:.3e}")

    print()
    print("=" * 70)
    print("C. Заявленный «фундаментальный предел» 1/НОК(3,8)")
    print("=" * 70)
    L = lcm(3, 8)
    print(f"  НОК(3, 8)      = {L}")
    print(f"  1/НОК(3, 8)    = {1/L:.4f}  ({1/L*100:.2f}%)  — заявлено в салфетке")
    print(f"  Реальная плотность запрещённых при N~10^39 = "
          f"{len(gaps)/10**39:.3e}")
    print()
    print("  ВЫВОД: 1/НОК(3,8) НЕ равно плотности запрещённых состояний")
    print("  полугруппы <3,8>. Совпадение 4.17% ≈ 4.7% — численное, а не")
    print("  выводимое из спектра. Подробности — в отчёте.")


if __name__ == "__main__":
    main()
