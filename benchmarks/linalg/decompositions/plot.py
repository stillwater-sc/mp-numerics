#!/usr/bin/env python3
"""Plot the decompositions study's CSVs as the SVG figures in figures/.

    python3 plot.py [results-dir] [figures-dir]

Needs matplotlib. The SVGs carry no date and a fixed hash salt, so the same
CSVs give the same files.
"""
import csv
import math
import os
import sys
from collections import OrderedDict

import matplotlib

matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'mp-numerics-14'
matplotlib.rcParams['font.size'] = 8
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
QR = ['householder_mtl5', 'householder_beta', 'householder_textbook', 'givens_mtl5']
LS = ['qr_householder', 'qr_householder_beta', 'qr_givens', 'normal_equations']
COLOR = {
    'householder_mtl5': '#2ca02c', 'householder_beta': '#ff7f0e', 'householder_textbook': '#1f77b4',
    'givens_mtl5': '#d62728',
    'qr_householder': '#2ca02c', 'qr_householder_beta': '#ff7f0e', 'qr_givens': '#d62728',
    'normal_equations': '#9467bd',
}
COLS = 4


def load(d, name):
    with open(os.path.join(d, name + '.csv'), newline='') as fh:
        return list(csv.DictReader(fh))


def save(fig, d, name):
    path = os.path.join(d, name + '.svg')
    fig.savefig(path, metadata={'Date': None, 'Creator': None}, bbox_inches='tight')
    plt.close(fig)
    print('wrote', path)


def value(x):
    """A CSV number as a float, or None for a row with no finite samples."""
    v = float(x)
    return v if math.isfinite(v) else None


def panels(types):
    rows = (len(types) + COLS - 1) // COLS
    fig, axes = plt.subplots(rows, COLS, figsize=(2.6 * COLS, 2.1 * rows), squeeze=False, sharex=True)
    for ax in axes.flat[len(types):]:
        ax.set_visible(False)
    return fig, list(axes.flat)


def qr_figure(rows, n, key, title, d, name):
    """key in eps against the generated kappa: randsvd solid, graded dashed."""
    rows = [r for r in rows if int(r['n']) == n and r['matrix'] in ('randsvd', 'graded')]
    types = list(OrderedDict.fromkeys(r['type'] for r in rows))
    fig, axes = panels(types)
    for ax, t in zip(axes, types):
        for m in QR:
            for cls, ls in (('randsvd', '-'), ('graded', '--')):
                rs = sorted((r for r in rows if r['type'] == t and r['method'] == m and r['matrix'] == cls),
                            key=lambda r: float(r['param']))
                pts = [(float(r['param']), value(r[key])) for r in rs]
                pts = [(x, max(y / float(rs[0]['eps']), 0.1)) for x, y in pts if y is not None]
                if pts:
                    ax.plot(*zip(*pts), ls, color=COLOR[m], lw=1,
                            label=m if cls == 'randsvd' else None)
        ax.set_yscale('log')
        ax.set_title(t)
        ax.axvline(math.log10(1 / float(next(r['eps'] for r in rows if r['type'] == t))), color='0.7', lw=0.6)
        ax.grid(True, which='major', lw=0.3)
    for ax in axes[len(types) - COLS:len(types)]:
        ax.set_xlabel('log10 kappa')
    fig.suptitle('%s, n = %d (eps; solid randsvd, dashed graded, grey line 1/eps)' % (title, n))
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=4,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, d, name)


def ls_figure(rows, m, n, d, name):
    """Forward error against the stored kappa, randsvd: zero residual solid,
    large residual dashed; grey lines kappa eps and kappa^2 eps."""
    rows = [r for r in rows if int(r['m']) == m and int(r['n']) == n and r['matrix'] == 'randsvd']
    types = list(OrderedDict.fromkeys(r['type'] for r in rows))
    fig, axes = panels(types)
    for ax, t in zip(axes, types):
        eps = float(next(r['eps'] for r in rows if r['type'] == t))
        for x in LS:
            for res, ls in (('zero', '-'), ('large', '--')):
                rs = sorted((r for r in rows if r['type'] == t and r['method'] == x and r['residual'] == res),
                            key=lambda r: float(r['param']))
                pts = [(value(r['kappa_stored_p50']), value(r['fwd_max'])) for r in rs]
                pts = [(math.log10(k), max(y, 1e-40)) for k, y in pts if k is not None and y is not None]
                if pts:
                    ax.plot(*zip(*pts), ls, color=COLOR[x], lw=1, label=x if res == 'zero' else None)
        ks = [value(r['kappa_stored_p50']) for r in rows if r['type'] == t]
        top = math.log10(max(k for k in ks if k is not None))
        xs = [0, top]
        ax.plot(xs, [eps * 10 ** v for v in xs], color='0.6', lw=0.6)
        ax.plot(xs, [eps * 10 ** (2 * v) for v in xs], ':', color='0.6', lw=0.6)
        ax.set_yscale('log')
        ax.set_ylim(max(eps / 100, 1e-40), 1e3)
        ax.set_title(t)
        ax.grid(True, which='major', lw=0.3)
    for ax in axes[len(types) - COLS:len(types)]:
        ax.set_xlabel('log10 kappa of the stored A')
    fig.suptitle('Least squares %d x %d, randsvd: relative forward error '
                 '(solid zero residual, dashed large; grey kappa eps, dotted kappa^2 eps)' % (m, n))
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=4,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, d, name)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results')
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'figures')
    os.makedirs(out, exist_ok=True)
    qr = load(d, 'qr')
    ls = load(d, 'lstsq')
    n = max(int(r['n']) for r in qr)
    qr_figure(qr, n, 'bwd_col_max', 'QR columnwise backward error', out, 'qr_backward_n%d' % n)
    qr_figure(qr, n, 'orth_max', 'QR orthogonality ||Q^T Q - I||', out, 'qr_orthogonality_n%d' % n)
    m, k = max({(int(r['m']), int(r['n'])) for r in ls})
    ls_figure(ls, m, k, out, 'lstsq_forward_%dx%d' % (m, k))


if __name__ == '__main__':
    main()
