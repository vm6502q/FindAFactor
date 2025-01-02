import os
import _find_a_factor

def find_a_factor(n,
                  use_congruence_of_squares=bool(os.environ.get('FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES')) if os.environ.get('FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES') else True,
                  node_count=int(os.environ.get('FINDAFACTOR_NODE_COUNT')) if os.environ.get('FINDAFACTOR_NODE_COUNT') else 1,
                  node_id=int(os.environ.get('FINDAFACTOR_NODE_ID')) if os.environ.get('FINDAFACTOR_NODE_ID') else 0,
                  smoothness_bound_multiplier=float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0):
    return int(_find_a_factor._find_a_factor(str(n), use_congruence_of_squares, node_count, node_id, smoothness_bound_multiplier))