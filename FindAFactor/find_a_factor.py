import _find_a_factor

def find_a_factor(n, use_congruence_of_squares=True, node_count=1, node_id=0, wheel_factorization_level=17):
    return int(_find_a_factor._find_a_factor(str(n), use_congruence_of_squares, node_count, node_id, wheel_factorization_level))