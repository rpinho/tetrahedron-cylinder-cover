#!/usr/bin/env python3
"""Exact obligations for the arbitrary-positive-radius extension.

Alexander Richard Towell, with substantial OpenAI Codex assistance.
MIT license; predecessor geometry/certificate: Ricardo Pinho.
This reuses the accepted exact integral checker and adds the endpoint
polynomial check. The continuum argument is written in INDEPENDENT_RADII.md
(PROOF.md in the research archive); it is not formalized by this program.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import structural_independent as previous


def endpoint_margin(sector):
    """The exact E=4D-3d(1+w) at q=l, tau=1."""
    a, f = sector.a, sector.a-sector.b*sector.l
    return 1-3*a*a+f*(4*a-3-3*a*a)


def check_endpoints(sectors):
    """Check the new hypothesis; the caller supplies valid geometric sectors."""
    if not sectors:
        raise previous.VerificationError('No sectors for the endpoint check')
    margins = []
    for sector in sectors:
        margin = endpoint_margin(sector)
        if margin < 0:
            raise previous.VerificationError(
                f'Endpoint estimate fails on [{sector.l},{sector.r}]: {margin}')
        margins.append(margin)
    index = min(range(len(sectors)), key=margins.__getitem__)
    return {'sectors_checked': len(sectors),
            'endpoint_ratio_lower': '3/2',
            'minimum_endpoint_margin': str(margins[index]),
            'minimum_sector': [str(sectors[index].l), str(sectors[index].r)],
            'continuum_justification': 'E>=0 plus the written q/tau monotonicity argument'}


def verify(data):
    # Validate coverage before any numerical replay. Source lower values are ignored.
    previous.read_partition(data)
    family = previous.FINAL_FAMILY
    sectors = previous.make_sectors(family.n, family.shape)
    endpoint_check = check_endpoints(sectors)
    result = previous.verify(data)
    result['scope'] = ('fixed final secant axes; arbitrary positive node-radius '
                       'profiles in the two families, without evenness; 0<=tau<=1')
    result['endpoint_check'] = endpoint_check
    result['proof_kind'] = 'written continuum argument plus exact finite arithmetic; not a formal proof'
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', type=Path,
                        default=Path(__file__).with_name('structural_final_family.json'))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--negative-control', choices=['endpoint', 'partition'])
    args = parser.parse_args(argv)
    try:
        if args.negative_control == 'endpoint':
            # This sector satisfies the old geometric inequalities, but not E>=0.
            bad = previous.Sector(F(3, 4), F(1), F(9, 10), F(9, 10))
            check_endpoints((bad,))
            print('RESULT: negative control was accepted')
            return 0
        raw = args.certificate.read_bytes()
        data = json.loads(raw)
        if args.negative_control == 'partition':
            data['accepted_intervals'].pop(0)
        result = verify(data)
        result['source_certificate_sha256'] = hashlib.sha256(raw).hexdigest()
        result['integral_checker_sha256'] = hashlib.sha256(
            Path(previous.__file__).read_bytes()).hexdigest()
        result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        if args.output:
            args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k: v for k, v in result.items() if k != 'interval_bounds'}, indent=2))
        print('RESULT: certified C >= 4993/10000 for arbitrary positive radius profiles')
        return 0
    except (ValueError, TypeError, KeyError, OSError) as error:
        print(f'FAIL: {error}')
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
