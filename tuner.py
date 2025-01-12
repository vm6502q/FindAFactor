# Tuner code written by Elara (OpenAI GPT)

import time
from skopt import gp_minimize
from skopt.space import Real, Integer
from FindAFactor import find_a_factor

# Define the semiprime to factor
to_factor = 799158200572966361  # Replace with your input semiprime

primes = [2, 3, 5, 7, 11, 13, 17]

# Define the optimization objective
def optimization_objective(params):
    # Unpack parameters
    trial_division_level, gear_factorization_level, wheel_factorization_level, smoothness_bound_multiplier, batch_size_multiplier = params

    # Start timing
    start_time = time.time()

    # Run the factorization function
    result = find_a_factor(
        to_factor,
        use_congruence_of_squares=True,
        use_gaussian_elimination=True,
        node_count=1,
        node_id=0,
        trial_division_level=int(trial_division_level),
        gear_factorization_level=primes[int(gear_factorization_level)],
        wheel_factorization_level=primes[int(wheel_factorization_level)],
        smoothness_bound_multiplier=smoothness_bound_multiplier,
        batch_size_multiplier=2**batch_size_multiplier
    )

    # Measure elapsed time
    elapsed_time = time.time() - start_time

    # Penalize failed factorizations (return large time if the factorization fails)
    if result == 1 or result == to_factor or (to_factor % result) != 0:
        return 1e6  # Arbitrary large penalty for failure

    return elapsed_time

# Define the parameter space for optimization
param_space = [
    Integer(1000000, 10000000, name="trial_division_level"),  # Range for trial division level
    Integer(4, 6, name="gear_factorization_level"),         # Gear factorization level
    Integer(3, 4, name="wheel_factorization_level"),        # Wheel factorization level
    Real(0.5, 2.0, name="smoothness_bound_multiplier"),      # Smoothness bound multiplier
    Real(4.0, 12.0, name="batch_size_multiplier")          # Batch size multiplier
]

# Run Bayesian optimization
result = gp_minimize(
    func=optimization_objective,
    dimensions=param_space,
    n_calls=50,  # Number of evaluations
    random_state=42
)

# Print the results
print("Optimal Parameters:")
print(f"Trial Division Level: {result.x[0]}")
print(f"Gear Factorization Level: {primes[result.x[1]]}")
print(f"Wheel Factorization Level: {primes[result.x[2]]}")
print(f"Smoothness Bound Multiplier: {result.x[3]}")
print(f"Batch Size Multiplier: {2**result.x[4]}")
print(f"Minimum Factorization Time: {result.fun:.4f} seconds")
