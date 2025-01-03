import os
import sys
import time

from FindAFactor import find_a_factor

def main():
    argv_len = len(sys.argv)

    if argv_len < 2:
        print("Enter a number to factor on the command line when calling " + sys.argv[0])
        return 1

    to_factor = int(sys.argv[1])
    use_congruence_of_squares = bool(os.environ.get('FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES')) if os.environ.get('FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES') else True
    node_count = int(os.environ.get('FINDAFACTOR_NODE_COUNT')) if os.environ.get('FINDAFACTOR_NODE_COUNT') else 1
    node_id = int(os.environ.get('FINDAFACTOR_NODE_ID')) if os.environ.get('FINDAFACTOR_NODE_ID') else 0
    gear_factorization_level = int(os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL') else 13
    wheel_factorization_level = int(os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL') else 7
    smoothness_bound_multiplier = float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0

    if argv_len > 2:
        use_congruence_of_squares = bool(sys.argv[2])
    if argv_len > 4:
        node_count = int(sys.argv[3])
        node_id = int(sys.argv[4])
    if argv_len > 5:
        gear_factorization_level = int(sys.argv[5])
    if argv_len > 6:
        wheel_factorization_level = int(sys.argv[6])
    if argv_len > 7:
        smoothness_bound_multiplier = float(sys.argv[7])

    start = time.perf_counter()
    result = find_a_factor(
        to_factor,
        use_congruence_of_squares = use_congruence_of_squares,
        node_count = node_count,
        node_id = node_id,
        gear_factorization_level = gear_factorization_level,
        wheel_factorization_level = wheel_factorization_level,
        smoothness_bound_multiplier = smoothness_bound_multiplier
    )
    print(time.perf_counter() - start)
    print(str(result) + " * " + str(to_factor // result) + " == " + str(to_factor))
    print((result * (to_factor // result)) == to_factor)

    return 0

if __name__ == '__main__':
    sys.exit(main())