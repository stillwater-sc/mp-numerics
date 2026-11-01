#!/usr/bin/env python3
"""Summarize the decompositions study's CSVs as the Markdown tables in
README.md. Standard library only.

    python3 summarize.py [results-dir]      (default: results/ next to this file)
"""
import csv
import math
import os
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
QR = ['householder_mtl5', 'householder_beta', 'householder_textbook', 'givens_mtl5']
LS = ['qr_householder', 'qr_householder_beta', 'qr_givens', 'normal_equations']
USABLE = 0.1   # a least-squares solution is usable while its relative forward error is at most 10%


def load(d, name):
    with open(os.path.join(d, name + '.csv'), newline='') as fh:
        return list(csv.DictReader(fh))


def f(x):
    return float(x)


def eps_units(x, eps):
    if x != x:
        return '-'
    if math.isinf(x):
        return 'inf'
    v = x / eps
    return '%.3g' % v if v < 1e4 else '%.1e' % v


def failures(r):
    """Samples that gave no measurable result, and why."""
    parts = []
    for key, label in (('nonfinite', 'non-finite'), ('thrown', 'thrown'), ('breakdown', 'breakdown')):
        if key in r and int(r[key]):
            parts.append('%s %s' % (r[key], label))
    return parts


def types(rows):
    return list(OrderedDict.fromkeys(r['type'] for r in rows))


def qr_table(rows, n):
    """Per type and method, over randsvd and graded at every kappa the type
    runs: the worst columnwise backward error and orthogonality, in eps, and
    the samples that failed."""
    out = ['| type | ' + ' | '.join('`%s`' % m for m in QR) + ' |', '|---|' + '---|' * len(QR)]
    for t in types(rows):
        cells = []
        for m in QR:
            rs = [r for r in rows if r['type'] == t and r['method'] == m and int(r['n']) == n
                  and r['matrix'] in ('randsvd', 'graded')]
            eps = f(rs[0]['eps'])
            col = max((f(r['bwd_col_max']) for r in rs if r['bwd_col_max'] != 'nan'), default=float('nan'))
            orth = max((f(r['orth_max']) for r in rs if r['orth_max'] != 'nan'), default=float('nan'))
            bad = sum(int(r['nonfinite']) + int(r['thrown']) for r in rs)
            ident = sum(int(r['q_identity']) for r in rs)
            cell = '%s / %s' % (eps_units(col, eps), eps_units(orth, eps))
            if bad:
                cell += ' (%d failed)' % bad
            if ident:
                cell += ' (%d Q = I)' % ident
            cells.append(cell)
        out.append('| %s | %s |' % (t, ' | '.join(cells)))
    return '\n'.join(out)


def reach(rs):
    """The largest kappa_stored at which every sample is usable, scanning up
    from kappa = 1 and stopping at the first row that is not: the condition
    number up to which the method can be trusted in this type."""
    best = None
    for r in sorted(rs, key=lambda r: f(r['param'])):
        ok = (int(r['singular']) == 0 and not failures(r) and r['fwd_max'] != 'nan'
              and f(r['fwd_max']) <= USABLE)
        if not ok:
            break
        best = f(r['kappa_stored_p50'])
    return best


def kappa_cell(k):
    if k is None:
        return 'none'
    return '1e%d' % math.floor(math.log10(k) + 1e-9) if k >= 10 else '%.3g' % k


def ls_table(rows, m, n, residual):
    """Per type and method, randsvd: the largest stored condition number up to
    which every solution has a forward error of at most 10%."""
    out = ['| type | 1/eps | ' + ' | '.join('`%s`' % x for x in LS) + ' |', '|---|---|' + '---|' * len(LS)]
    for t in types(rows):
        cells = []
        eps = None
        for x in LS:
            rs = [r for r in rows if r['type'] == t and r['method'] == x and r['matrix'] == 'randsvd'
                  and r['residual'] == residual and int(r['m']) == m and int(r['n']) == n]
            eps = f(rs[0]['eps'])
            cells.append(kappa_cell(reach(rs)))
        out.append('| %s | %.1e | %s |' % (t, 1 / eps, ' | '.join(cells)))
    return '\n'.join(out)


def fixed_table(rows, cls, m, n):
    """Forward error of each method on one deterministic-in-shape class
    (hilbert, or kahan at its smallest theta), zero residual."""
    rs = [r for r in rows if r['matrix'] == cls and r['residual'] == 'zero' and int(r['m']) == m
          and int(r['n']) == n]
    if cls == 'kahan':
        theta = min(f(r['param']) for r in rs)
        rs = [r for r in rs if f(r['param']) == theta]
    out = ['| type | kappa stored | ' + ' | '.join('`%s`' % x for x in LS) + ' |',
           '|---|---|' + '---|' * len(LS)]
    for t in types(rs):
        cells = []
        ks = None
        for x in LS:
            r = next(r for r in rs if r['type'] == t and r['method'] == x)
            ks = r['kappa_stored_p50']
            why = failures(r)
            if int(r['singular']):
                cells.append('singular')
            elif why:
                cells.append(', '.join(why))
            else:
                cells.append('%.2g' % f(r['fwd_max']))
        out.append('| %s | %.2g | %s |' % (t, f(ks), ' | '.join(cells)))
    return '\n'.join(out)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results')
    qr = load(d, 'qr')
    ls = load(d, 'lstsq')
    for n in sorted({int(r['n']) for r in qr}):
        print('## QR, n = %d: worst columnwise backward error / orthogonality (eps)\n' % n)
        print(qr_table(qr, n))
        print()
    for m, n in sorted({(int(r['m']), int(r['n'])) for r in ls}):
        for residual in ('zero', 'large'):
            print('## Least squares, %d x %d, %s residual: largest kappa with forward error <= 10%%\n'
                  % (m, n, residual))
            print(ls_table(ls, m, n, residual))
            print()
        for cls in ('hilbert', 'kahan'):
            print('## Least squares, %d x %d, %s: forward error\n' % (m, n, cls))
            print(fixed_table(ls, cls, m, n))
            print()


if __name__ == '__main__':
    main()
