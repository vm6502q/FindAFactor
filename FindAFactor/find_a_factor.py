import os
import _find_a_factor
from enum import IntEnum


class FactoringMethod(IntEnum):
    PRIME_SOLVER = 0
    MIXED = 1
    FACTOR_FINDER = 2


def find_a_factor(n,
                  method=FactoringMethod(int(os.environ.get('FINDAFACTOR_METHOD'))) if os.environ.get('FINDAFACTOR_METHOD') else FactoringMethod.PRIME_SOLVER,
                  node_count=int(os.environ.get('FINDAFACTOR_NODE_COUNT')) if os.environ.get('FINDAFACTOR_NODE_COUNT') else 1,
                  node_id=int(os.environ.get('FINDAFACTOR_NODE_ID')) if os.environ.get('FINDAFACTOR_NODE_ID') else 0,
                  gear_factorization_level=int(os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL') else 11,
                  wheel_factorization_level=int(os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL') else 11,
                  sieving_bound_multiplier=float(os.environ.get('FINDAFACTOR_SIEVING_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SIEVING_BOUND_MULTIPLIER') else 1.0,
                  smoothness_bound_multiplier=float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0,
                  gaussian_elimination_row_multiplier=float(os.environ.get('FINDAFACTOR_GAUSSIAN_ELIMINATION_ROW_MULTIPLIER')) if os.environ.get('FINDAFACTOR_GAUSSIAN_ELIMINATION_ROW_MULTIPLIER') else 2.0,
                  skip_trial_division=True if os.environ.get('FINDAFACTOR_SKIP_TRIAL_DIVISION') else False):
    return int(_find_a_factor._find_a_factor(str(n),
                                             int(method),
                                             node_count, node_id,
                                             gear_factorization_level,
                                             wheel_factorization_level,
                                             sieving_bound_multiplier,
                                             smoothness_bound_multiplier,
                                             gaussian_elimination_row_multiplier,
                                             skip_trial_division))