import _find_a_factor

def find_a_factor(n, use_congruence_of_squares=True, node_count=1, node_id=0, gear_factorization_level=13, wheel_factorization_level=11, smoothness_bound_multiplier=1.0):
    return int(_find_a_factor._find_a_factor(str(n), use_congruence_of_squares, node_count, node_id, gear_factorization_level, wheel_factorization_level, smoothness_bound_multiplier))