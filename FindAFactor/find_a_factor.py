import os
import _find_a_factor

def find_a_factor(n,
                  use_congruence_of_squares=True if os.environ.get('FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES') else False,
                  node_count=int(os.environ.get('FINDAFACTOR_NODE_COUNT')) if os.environ.get('FINDAFACTOR_NODE_COUNT') else 1,
                  node_id=int(os.environ.get('FINDAFACTOR_NODE_ID')) if os.environ.get('FINDAFACTOR_NODE_ID') else 0,
                  gear_factorization_level=int(os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL') else 11,
                  wheel_factorization_level=int(os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL') else 7,
                  smoothness_bound_multiplier=float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0,
                  batch_size_multiplier=float(os.environ.get('FINDAFACTOR_BATCH_SIZE_MULTIPLIER')) if os.environ.get('FINDAFACTOR_BATCH_SIZE_MULTIPLIER') else 512.0):
    return int(_find_a_factor._find_a_factor(str(n),
                                             use_congruence_of_squares,
                                             node_count, node_id,
                                             gear_factorization_level,
                                             wheel_factorization_level,
                                             smoothness_bound_multiplier,
                                             batch_size_multiplier))