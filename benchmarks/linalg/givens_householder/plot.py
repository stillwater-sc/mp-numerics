#!/usr/bin/env python3
"""Plot the givens_householder study's CSVs as the SVG figures in figures/.

    python3 plot.py [results-dir] [figures-dir]

Needs matplotlib. The SVGs carry no date and a fixed hash salt, so the same
CSVs give the same files.
"""
import csv
import math
import os
import sys
from collections import OrderedDict, defaultdict

import matplotlib

matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'mp-numerics-13'
matplotlib.rcParams['font.size'] = 9
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OPS = ['givens_mtl5', 'givens_max_scaled', 'givens_safe_scaled', 'givens_pow2_scaled',
       'householder_mtl5', 'householder_textbook', 'householder_lapack_pow2']
ROTATIONS = OPS[:4]
REFLECTORS = OPS[4:]
STYLE = {
    'givens_mtl5': ('#d62728', '-'), 'givens_max_scaled': ('#ff7f0e', '--'),
    'givens_safe_scaled': ('#9467bd', '-.'), 'givens_pow2_scaled': ('#2ca02c', '-'),
    'householder_mtl5': ('#d62728', '-'), 'householder_textbook': ('#1f77b4', '--'),
    'householder_lapack_pow2': ('#2ca02c', '-'),
}
FLOOR = 1e-1   # errors below 0.1 eps are drawn at the floor of the log axis


def load(d, name):
    with open(os.path.join(d, name + '.csv'), newline='') as fh:
        return list(csv.DictReader(fh))


def save(fig, d, name):
    path = os.path.join(d, name + '.svg')
    fig.savefig(path, metadata={'Date': None, 'Creator': None}, bbox_inches='tight')
    plt.close(fig)
    print('wrote', path)


def by_type(rows):
    g = OrderedDict()
    for r in rows:
        if int(r['samples']) == 0:   # the bottom of a range that flushes to zero
            continue
        g.setdefault(r['type'], defaultdict(list))[r['operator']].append(r)
    return g


def clean(r):
    return int(r['gross']) + int(r['nonfinite']) == 0


def finite_eps(r):
    """orth_max in eps of a clean row, or None."""
    if not clean(r) or r['orth_max'] == 'nan':
        return None
    return max(float(r['orth_max']) / float(r['eps']), FLOOR)


def gross_level(rows):
    """Where gross and non-finite rows are drawn: two decades above the
    largest finite error among the rows plotted, so none is clipped."""
    finite = [v for v in map(finite_eps, rows) if v is not None]
    top = max(finite, default=1.0)
    return 10.0 ** (math.ceil(math.log10(max(top, 1.0))) + 2)


def clean_share(rs):
    """The share of the exponent range on which every sample is clean. Each
    sampled k stands for the half-intervals to its neighbours, so the closely
    spaced exponents scale_grid adds at each end weigh what they cover."""
    rs = sorted(rs, key=lambda r: float(r['param']))
    ks = [float(r['param']) for r in rs]
    if len(ks) == 1:
        return float(clean(rs[0]))
    w = [((ks[min(i + 1, len(ks) - 1)] - ks[max(i - 1, 0)]) / 2) for i in range(len(ks))]
    return sum(wi for wi, r in zip(w, rs) if clean(r)) / sum(w)


def heatmap(scale_rows, d):
    """Share of each type's exponent range on which every sample is clean."""
    g = by_type(scale_rows)
    types = list(g)
    data = []
    for t in types:
        row = []
        for o in OPS:
            rs = g[t].get(o, [])
            row.append(clean_share(rs) if rs else float('nan'))
        data.append(row)
    fig, ax = plt.subplots(figsize=(7.2, 0.2 * len(types) + 1.6))
    im = ax.imshow(data, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1)
    ax.set_xticks(range(len(OPS)), [o.replace('_', '\n', 1) for o in OPS], fontsize=7)
    ax.set_yticks(range(len(types)), types, fontsize=7)
    for i, row in enumerate(data):
        for j, v in enumerate(row):
            if v == v:
                ax.text(j, i, '%d' % round(100 * v), ha='center', va='center', fontsize=6)
    ax.set_title('Scale sweep: % of the exponent range where every sample is clean')
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    save(fig, d, 'scale_clean_range')


def family_panels(rows, families, ops, xkey, xlabel, title, name, d, normalize_x=False):
    """One panel per type, one line per operator: worst orthogonality in eps."""
    g = by_type(rows)
    types = [t for t in g if g[t][ops[0]] and g[t][ops[0]][0]['family'] in families]
    if not types:
        return
    cols = 4
    nrows = math.ceil(len(types) / cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(10, 1.9 * nrows + 1.0), squeeze=False, sharey=True,
                             layout='constrained')
    gross = gross_level([r for t in types for o in ops for r in g[t].get(o, [])])
    for ax, t in zip(axes.flat, types):
        for o in ops:
            rs = g[t].get(o, [])
            if not rs:
                continue
            xs = [float(r[xkey]) for r in rs]
            if normalize_x:
                lo, hi = min(xs), max(xs)
                xs = [(x - lo) / (hi - lo) if hi > lo else 0.0 for x in xs]
            color, ls = STYLE[o]
            ys = [finite_eps(r) for r in rs]
            ax.plot(xs, [gross if y is None else y for y in ys], color=color, ls=ls, lw=1, label=o)
        ax.set_yscale('log')
        ax.set_ylim(FLOOR / 2, gross * 10)
        ax.axhline(gross / 10, color='#999', lw=0.5)
        ax.set_title(t, fontsize=8)
        ax.tick_params(labelsize=6)
    for ax in axes.flat[len(types):]:
        ax.axis('off')
    fig.supxlabel(xlabel, fontsize=8)
    fig.supylabel('worst ||Q^T Q - I|| / eps  (top band: gross or non-finite)', fontsize=8)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncol=len(ops), fontsize=7, frameon=False,
               title=title, title_fontsize=10)
    save(fig, d, name)


def repeat_panels(rows, d):
    """Worst deviation in eps against k, per operator, one panel per mode;
    one line per type, coloured by family."""
    fam_color = {'ieee': '#1f77b4', 'posit': '#d62728', 'fixpnt': '#8c564b', 'cfloat': '#7f7f7f',
                 'lns': '#9467bd', 'takum': '#e377c2', 'dd': '#17becf', 'qd': '#bcbd22'}
    for mode in ('power', 'product'):
        fig, axes = plt.subplots(2, 4, figsize=(10, 5), squeeze=False, sharey=True, layout='constrained')
        series = defaultdict(list)
        for r in rows:
            if r['mode'] == mode:
                series[(r['operator'], r['type'], r['family'])].append(r)
        for ax, o in zip(axes.flat, OPS):
            for (op, t, fam), rs in series.items():
                if op != o:
                    continue
                # dev_max is over the trials that stayed finite; a checkpoint
                # where some overflowed is still drawn, and marked with an x.
                pts = [(int(r['k']), float(r['dev_max']) / float(r['eps']), int(r['nonfinite']) > 0)
                       for r in rs if r['dev_max'] != 'nan']
                if pts:
                    color = fam_color.get(fam, 'k')
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=0.7, alpha=0.7)
                    over = [p for p in pts if p[2]]
                    if over:
                        ax.plot([p[0] for p in over], [p[1] for p in over], ls='none', marker='x',
                                ms=3, color=color)
            ax.set_xscale('log', base=2)
            ax.set_yscale('log')
            ax.set_title(o, fontsize=8)
            ax.tick_params(labelsize=6)
        axes.flat[-1].axis('off')
        for fam, c in fam_color.items():
            axes.flat[-1].plot([], [], color=c, label=fam)
        axes.flat[-1].legend(loc='center', fontsize=8, frameon=False)
        what = ('one operator applied k times to unit vectors' if mode == 'power'
                else 'product of k fresh operators on an 8 x 8 Q')
        fig.suptitle('Repeated application (%s): worst deviation / eps over the finite trials; '
                     'x: some trials overflowed' % what, fontsize=9)
        fig.supxlabel('k', fontsize=8)
        save(fig, d, 'repeat_' + mode)


def main():
    rdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results')
    fdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'figures')
    os.makedirs(fdir, exist_ok=True)
    scale, tail, repeat = load(rdir, 'scale'), load(rdir, 'tail'), load(rdir, 'repeat')
    heatmap(scale, fdir)
    groups = [('ieee_posit_fixpnt', {'ieee', 'posit', 'fixpnt'}), ('cfloat', {'cfloat'}),
              ('lns_takum_dd_qd', {'lns', 'takum', 'dd', 'qd'})]
    for tag, fams in groups:
        family_panels(scale, fams, ROTATIONS, 'param', 'position of ||x|| in the exponent range (0 = denorm_min, 1 = max)',
                      'Scale sweep, rotations: worst orthogonality', 'scale_rotations_' + tag, fdir, True)
        family_panels(scale, fams, REFLECTORS, 'param', 'position of ||x|| in the exponent range (0 = denorm_min, 1 = max)',
                      'Scale sweep, reflectors: worst orthogonality', 'scale_reflectors_' + tag, fdir, True)
        family_panels(tail, fams, REFLECTORS, 'param', 'j, where x = (1, 2^-j)',
                      'Tail sweep, reflectors: worst orthogonality', 'tail_reflectors_' + tag, fdir)
    repeat_panels(repeat, fdir)


if __name__ == '__main__':
    main()
