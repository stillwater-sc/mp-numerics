#!/usr/bin/env python3
"""Summarize the givens_householder study's CSVs as the Markdown tables in
README.md. Standard library only.

    python3 summarize.py [results-dir]      (default: results/ next to this file)
"""
import csv
import os
import sys
from collections import OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OPS = ['givens_mtl5', 'givens_max_scaled', 'givens_safe_scaled', 'givens_pow2_scaled',
       'householder_mtl5', 'householder_textbook', 'householder_lapack_pow2']


def load(d, name):
    with open(os.path.join(d, name + '.csv'), newline='') as fh:
        return list(csv.DictReader(fh))


def num(x):
    return float(x)


def clean(r):
    return int(r['gross']) + int(r['nonfinite']) == 0


def eps_units(x, eps):
    if x != x:
        return '-'
    v = x / eps
    return '%.3g' % v if v < 1e4 else '%.1e' % v


def grouped(rows):
    g = OrderedDict()
    for r in rows:
        g.setdefault(r['type'], defaultdict(list))[r['operator']].append(r)
    return g


def scale_table(rows):
    """Per type and operator: the exponents k (||x|| = 2^k) where every sample
    is clean, out of the type's range, and the worst orthogonality there."""
    out = ['| type | range of k | ' + ' | '.join('`%s`' % o for o in OPS) + ' |',
           '|---|---|' + '---|' * len(OPS)]
    for t, ops in grouped(rows).items():
        ks = [int(num(r['param'])) for r in next(iter(ops.values()))]
        cells = []
        for o in OPS:
            rs = ops.get(o, [])
            good = [int(num(r['param'])) for r in rs if clean(r)]
            if len(good) == len(rs):
                span = 'all'
            elif not good:
                span = 'none'
            else:
                span = '%d..%d' % (min(good), max(good))
            eps = num(rs[0]['eps'])
            worst = max((num(r['orth_max']) for r in rs if clean(r)), default=float('nan'))
            cells.append('%s, %s' % (span, eps_units(worst, eps)))
        out.append('| %s | %d..%d | %s |' % (t, min(ks), max(ks), ' | '.join(cells)))
    return '\n'.join(out)


def tail_table(rows):
    """Per type: worst orthogonality over the tail sweep, and the first j at
    which the operator is the identity."""
    out = ['| type | ' + ' | '.join('`%s`' % o for o in OPS) + ' |', '|---|' + '---|' * len(OPS)]
    for t, ops in grouped(rows).items():
        cells = []
        for o in OPS:
            rs = ops.get(o, [])
            eps = num(rs[0]['eps'])
            worst = max(num(r['orth_max']) for r in rs if r['orth_max'] != 'nan')
            ident = [int(num(r['param'])) for r in rs if int(r['identity']) > 0]
            cells.append(eps_units(worst, eps) + ('' if not ident else ' (I from j=%d)' % min(ident)))
        out.append('| %s | %s |' % (t, ' | '.join(cells)))
    return '\n'.join(out)


def repeat_table(rows):
    """dev_max at the largest k, in eps: power (one operator k times) and
    product (k fresh operators on an 8 x 8 Q)."""
    last = OrderedDict()
    for r in rows:
        last[(r['type'], r['operator'], r['mode'])] = r
    types = list(OrderedDict.fromkeys(r['type'] for r in rows))
    out = []
    for mode in ('power', 'product'):
        ks = {r['k'] for (t, o, m), r in last.items() if m == mode}
        out.append('%s, k = %s:\n' % (mode, '/'.join(sorted(ks))))
        out.append('| type | ' + ' | '.join('`%s`' % o for o in OPS) + ' |')
        out.append('|---|' + '---|' * len(OPS))
        for t in types:
            cells = []
            for o in OPS:
                r = last.get((t, o, mode))
                cells.append('-' if r is None else eps_units(num(r['dev_max']), num(r['eps'])))
            out.append('| %s | %s |' % (t, ' | '.join(cells)))
        out.append('')
    return '\n'.join(out)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results')
    print('## Scale sweep: clean exponents k, and worst orthogonality there (eps)\n')
    print(scale_table(load(d, 'scale')))
    print('\n## Tail sweep: worst orthogonality (eps)\n')
    print(tail_table(load(d, 'tail')))
    print('\n## Repeated application: worst deviation (eps)\n')
    print(repeat_table(load(d, 'repeat')))


if __name__ == '__main__':
    main()
