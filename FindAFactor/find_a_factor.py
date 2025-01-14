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
                  trial_division_level=int(os.environ.get('FINDAFACTOR_TRIAL_DIVISION_LEVEL')) if os.environ.get('FINDAFACTOR_TRIAL_DIVISION_LEVEL') else (1<<20),
                  gear_factorization_level=int(os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL') else 13,
                  wheel_factorization_level=int(os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL') else 11,
                  smoothness_bound_multiplier=float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0,
                  batch_size_multiplier=float(os.environ.get('FINDAFACTOR_BATCH_SIZE_MULTIPLIER')) if os.environ.get('FINDAFACTOR_BATCH_SIZE_MULTIPLIER') else 128.0):
    return int(_find_a_factor._find_a_factor(str(n),
                                             int(method),
                                             node_count, node_id,
                                             trial_division_level,
                                             gear_factorization_level,
                                             wheel_factorization_level,
                                             smoothness_bound_multiplier,
                                             batch_size_multiplier))