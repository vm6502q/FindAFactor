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
    use_congruence_of_squares = True if os.environ.get('FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES') else False
    use_gaussian_elimination = True if os.environ.get('FINDAFACTOR_USE_GAUSSIAN_ELIMINATION') else False
    node_count = int(os.environ.get('FINDAFACTOR_NODE_COUNT')) if os.environ.get('FINDAFACTOR_NODE_COUNT') else 1
    node_id = int(os.environ.get('FINDAFACTOR_NODE_ID')) if os.environ.get('FINDAFACTOR_NODE_ID') else 0
    gear_factorization_level = int(os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_GEAR_FACTORIZATION_LEVEL') else 11
    wheel_factorization_level = int(os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL')) if os.environ.get('FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL') else 7
    smoothness_bound_multiplier = float(os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER')) if os.environ.get('FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER') else 1.0
    batch_size_multiplier=float(os.environ.get('FINDAFACTOR_BATCH_SIZE_MULTIPLIER')) if os.environ.get('FINDAFACTOR_BATCH_SIZE_MULTIPLIER') else 1.0

    if argv_len > 2:
        use_congruence_of_squares = ((sys.argv[2] != "False") and (sys.argv[2] != "0"))
    if argv_len > 3:
        use_gaussian_elimination = ((sys.argv[3] != "False") and (sys.argv[3] != "0"))
    if argv_len > 5:
        node_count = int(sys.argv[4])
        node_id = int(sys.argv[5])
    if argv_len > 6:
        gear_factorization_level = int(sys.argv[6])
    if argv_len > 7:
        wheel_factorization_level = int(sys.argv[7])
    if argv_len > 8:
        smoothness_bound_multiplier = float(sys.argv[8])
    if argv_len > 9:
        batch_size_multiplier = float(sys.argv[9])

    start = time.perf_counter()
    result = find_a_factor(
        to_factor,
        use_congruence_of_squares = use_congruence_of_squares,
        use_gaussian_elimination = use_gaussian_elimination,
        node_count = node_count,
        node_id = node_id,
        gear_factorization_level = gear_factorization_level,
        wheel_factorization_level = wheel_factorization_level,
        smoothness_bound_multiplier = smoothness_bound_multiplier,
        batch_size_multiplier = batch_size_multiplier
    )
    print(time.perf_counter() - start)
    print(str(result) + " * " + str(to_factor // result) + " == " + str(to_factor))
    print((result * (to_factor // result)) == to_factor)

    return 0

if __name__ == '__main__':
    sys.exit(main())