#!/usr/bin/env python3
"""Optional publication figure; calculations/tests do not require matplotlib."""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import checks


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    result = checks.run_all()
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    ax = axes[0, 0]
    rows = checks.firas_rows()
    frequency = [r[0] * 100 * checks.C / 1e9 for r in rows]
    ax.errorbar(frequency, [r[2] for r in rows], yerr=[r[3] for r in rows],
                fmt='.', capsize=2, label='FIRAS: остаток ± σ (без ковариаций)')
    for tau in [1e-14, 1e-13]:
        correction = [1000 * (r[1] - r[2] / 1000)
                      * checks.spectral_fraction(nu * 1e9, tau)
                      for nu, r in zip(frequency, rows)]
        ax.plot(frequency, correction, label=f'Пример τ = {tau:.0e} s')
    ax.axhline(0, color='gray', linewidth=0.7)
    ax.set(xlabel='Частота, GHz', ylabel='Остаток, kJy/sr',
           title='АЧТ: опубликованные данные и выбранная поправка')
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    temperatures = list(range(100, 273))
    expansion = [checks.ice(t, 101325)['alpha_volume_K_inverse'] * 1e6
                 for t in temperatures]
    ax.plot(temperatures, expansion)
    ax.set(xlabel='Температура, K', ylabel='Объёмный αV, 10⁻⁶ K⁻¹',
           title='Лёд Ih: уравнение IAPWS, p = 101325 Pa')

    ax = axes[1, 0]
    atoms = result['atoms']['radii']
    ax.plot([r['Z'] for r in atoms], [r['reference_radius_pm'] for r in atoms],
            '.-', markersize=3, label='Ковалентные радиусы Pyykkö–Atsumi')
    ax.plot([r['Z'] for r in atoms], [r['unscreened_radius_pm'] for r in atoms],
            label='a₀ n²/Z; n = период (грубая оценка)')
    ax.plot([r['Z'] for r in atoms[:20]], [r['slater_radius_pm'] for r in atoms[:20]],
            'x--', label='Слэтер, Z ≤ 20 (оценка)')
    ax.set(xlabel='Атомный номер Z', ylabel='Радиус, pm',
           title='118 элементов: разные определения радиуса')
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    delays = list(range(-400, 401, 2))
    for visibility in [1, 0.5, 0]:
        ax.plot(delays, [checks.hom(t, 100, visibility) for t in delays],
                label=f'Задано V = {visibility}')
    ax.set(xlabel='Задержка, fs', ylabel='Вероятность совпадения', ylim=(-0.02, 0.55),
           title='HOM: стандартный модельный пример, σ = 100 fs')
    ax.legend(fontsize=8)

    fig.suptitle('Проверки к issue #14: данные, эталонные формулы и иллюстрации\n'
                 'Кривые не являются предсказаниями не заданного в переписке оператора',
                 fontsize=13)
    path = checks.ROOT / 'results/overview.png'
    fig.savefig(path, dpi=160, metadata={'Software': 'matplotlib; render_figures.py'})
    plt.close(fig)
    print('Wrote', path)


if __name__ == '__main__':
    main()
