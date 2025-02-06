import os
import _find_a_factor
from enum import IntEnum


class FactoringMethod(IntEnum):
    PRIME_PROVER = 0
    FACTOR_FINDER = 1


def find_a_factor(n,
                  method=FactoringMethod(int(os.environ.get('FINDAFACTOR_METHOD'))) if os.environ.get('FINDAFACTOR_METHOD') else FactoringMethod.PRIME_PROVER,
                  node_count=int(os.environ.get('FINDAFACTOR_NODE_COUNT')) if os.environ.get('FINDAFACTOR_NODE_COUNT') else 1,
                  node_id=int(os.environ.get('FINDAFACTOR_NODE_ID')) if os.environ.get('FINDAFACTOR_NODE_ID') else 0,
                  gear_factorization_level=int(os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL') else 23,
                  wheel_factorization_level=int(os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL') else 13,
                  sieving_bound_multiplier=float(os.environ.get('FINDAFACTOR_SIEVING_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SIEVING_BOUND_MULTIPLIER') else 1.0,
                  smoothness_bound_multiplier=float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0,
                  gaussian_elimination_row_offset=int(os.environ.get('FINDAFACTOR_GAUSSIAN_ELIMINATION_ROW_OFFSET')) if os.environ.get('FINDAFACTOR_GAUSSIAN_ELIMINATION_ROW_OFFSET') else 1,
                  check_small_factors=True if os.environ.get('FINDAFACTOR_CHECK_SMALL_FACTORS') else False,
                  wheel_primes_excluded=[int(i) for i in os.environ.get('FINDAFACTOR_WHEEL_PRIMES_EXCLUDED').split()] if os.environ.get('FINDAFACTOR_GAUSSIAN_ELIMINATION_ROW_OFFSET') else []):
    return int(_find_a_factor._find_a_factor(str(n),
                                             int(method),
                                             node_count, node_id,
                                             gear_factorization_level,
                                             wheel_factorization_level,
                                             sieving_bound_multiplier,
                                             smoothness_bound_multiplier,
                                             gaussian_elimination_row_offset,
                                             check_small_factors,
                                             wheel_primes_excluded))