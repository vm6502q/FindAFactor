import sys
import time

from FindAFactor import find_a_factor

def main():
    argv_len = len(sys.argv)

    if argv_len < 2:
        print("Enter a number to factor on the command line when calling " + sys.argv[0])
        return 1

    to_factor = int(sys.argv[1])
    use_congruence_of_squares = True
    node_count = 1
    node_id = 0
    wheel_factorization_level = 17

    if argv_len > 2:
        use_congruence_of_squares = bool(sys.argv[5])
    if argv_len > 4:
        node_count = int(sys.argv[2])
        node_id = int(sys.argv[3])
    if argv_len > 5:
        wheel_factorization_level = int(sys.argv[4])

    start = time.perf_counter()
    result = find_a_factor(
        to_factor,
        use_congruence_of_squares = use_congruence_of_squares,
        node_count = node_count,
        node_id = node_id,
        wheel_factorization_level = wheel_factorization_level
    )
    print(time.perf_counter() - start)
    print(str(result) + " * " + str(to_factor // result) + " == " + str(to_factor))
    print((result * (to_factor // result)) == to_factor)

    return 0

if __name__ == '__main__':
    sys.exit(main())