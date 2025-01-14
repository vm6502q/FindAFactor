# Tuner code written by Elara (OpenAI GPT) and improved by Dan Strano
# (Install scikit-optimize to use)

import sys
import time
from skopt import gp_minimize
from skopt.space import Real, Integer
from FindAFactor import find_a_factor

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

# Define the optimization objective
def optimization_objective(to_factor, params):
    # Unpack parameters
    method, trial_division_level, gear_factorization_level, wheel_factorization_level, smoothness_bound_multiplier, batch_size_multiplier = params

    # Start timing
    start_time = time.perf_counter()

    # Run the factorization function
    result = find_a_factor(
        to_factor,
        use_congruence_of_squares=(method > 0),
        node_count=1,
        node_id=0,
        trial_division_level=(1<<int(trial_division_level)),
        gear_factorization_level=primes[int(gear_factorization_level)],
        wheel_factorization_level=primes[int(wheel_factorization_level)],
        smoothness_bound_multiplier=2**smoothness_bound_multiplier,
        batch_size_multiplier=2**(batch_size_multiplier-gear_factorization_level)
    )

    # Measure elapsed time
    elapsed_time = time.perf_counter() - start_time

    # Penalize failed factorizations (return large time if the factorization fails)
    if result == 1 or result == to_factor or (to_factor % result) != 0:
        return 1e6  # Arbitrary large penalty for failure

    return elapsed_time

def main():
    # Define the semiprime to factor
    to_factor = 799158200572966361  # Arbitrary default semiprime

    argv_len = len(sys.argv)
    if argv_len > 1:
        to_factor = int(sys.argv[1])

    # Define the parameter space for optimization
    param_space = [
        Integer(0, 1, name="use_congruence_of_squares"),          # Enumeration of bool
        Integer(12, 24, name="trial_division_level"),             # Range for trial division level
        Integer(4, 6, name="gear_factorization_level"),           # Gear factorization level
        Integer(4, 5, name="wheel_factorization_level"),          # Wheel factorization level
        Real(-1.0, 1.0, name="smoothness_bound_multiplier"),      # Smoothness bound multiplier
        Real(10.0, 15.0, name="batch_size_multiplier")            # Batch size multiplier
    ]

    # Run Bayesian optimization
    result = gp_minimize(
        func=(lambda params : optimization_objective(to_factor, params)),
        dimensions=param_space,
        n_calls=50,  # Number of evaluations
        random_state=42
    )

    # Print the results
    print("Optimal Parameters:")
    print(f"Use Congruence of Squaes: {(result.x[0] > 0)}")
    print(f"Trial Division Level: {(1<<result.x[1])}")
    print(f"Gear Factorization Level: {primes[result.x[2]]}")
    print(f"Wheel Factorization Level: {primes[result.x[3]]}")
    print(f"Smoothness Bound Multiplier: {2**result.x[4]}")
    print(f"Batch Size Multiplier: {2**(result.x[5]-result.x[2])}")
    print(f"Minimum Factorization Time: {result.fun:.4f} seconds")

    return 0

if __name__ == '__main__':
    sys.exit(main())