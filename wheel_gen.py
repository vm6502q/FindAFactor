# This script demonstrates the sequence that is being
# enumerated by "reverse trial division" in the prime
# genertor script.

# Driver Code
if __name__ == '__main__':
    primes = [2, 3, 5, 7, 11]

    # This is the wheel level primes multiplied out.
    # wheel_size = 2
    # wheel_size = 6 # (2 * 3)
    # wheel_size = 30 # (2 * 3 * 5)
    # wheel_size = 210 # (2 * 3 * 5 * 7)
    wheel_size = 2310 # (2 * 3 * 5 * 7 * 11)
    # Include whatever primes are in the wheel, in wheel_size.

    for i in range(1, wheel_size):
        isNotMultiple = True
        for p in primes:
            if (i % p) == 0:
                isNotMultiple = False
                break
        if isNotMultiple:
            # Print the elements in the period
            # that are not small prime multiples.
            print(i, end=" ")
