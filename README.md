# FindAFactor
Find any nontrivial factor of a number

[![PyPI Downloads](https://static.pepy.tech/badge/findafactor)](https://pepy.tech/projects/findafactor)

## Copyright and license
(c) Daniel Strano and the Qrack contributors 2017-2025. All rights reserved.

## Installation
From PyPi:
```
pip3 install FindAFactor
```

From Source: install `pybind11`, then
```
pip3 install .
```
in the root source directory (with `setup.py`).

Windows users might find Windows Subsystem Linux (WSL) to be the easier and preferred choice for installation.

## Usage

```py
from FindAFactor import find_a_factor, FactoringMethod

to_factor = 1000

factor = find_a_factor(
    to_factor,
    method=FactoringMethod.PRIME_SOLVER,
    node_count=1, node_id=0,
    gear_factorization_level=11,
    wheel_factorization_level=11,
    sieving_bound_multiplier=1.0,
    smoothness_bound_multiplier=1.0,
    gaussian_elimination_row_offset=3,
    check_small_factors=False
)
```

The `find_a_factor()` function should return any nontrivial factor of `to_factor` (that is, any factor besides `1` or `to_factor`) if it exists. If a nontrivial factor does _not_ exist (i.e., the number to factor is prime), the function will return `1` or the original `to_factor`.

- `method` (default value: `PRIME_SOLVER`/`0`): `PRIME_SOLVER`/`0` will prove that a number is prime (by failing to find any factors with wheel and gear factorization). `FACTOR_FINDER`/`1` is optimized for the assumption that the number has at least two nontrivial factors.
- `node_count` (default value: `1`): `FindAFactor` can perform factorization in a _distributed_ manner, across nodes, without network communication! When `node_count` is set higher than `1`, the search space for factors is segmented equally per node. If the number to factor is semiprime, and brute-force search is used instead of congruence of squares, for example, all nodes except the one that happens to contain the (unknown) prime factor less than the square root of `to_factor` will ultimately return `1`, while one node will find and return this factor. For best performance, every node involved in factorization should have roughly the same CPU throughput capacity. For `FACTOR_FINDER` mode, this splits the sieving range between nodes, but it does not actually coordinate Gaussian elimination rows between nodes.
- `node_id` (default value: `0`): This is the identifier of this node, when performing distributed factorization with `node_count` higher than `1`. `node_id` values start at `0` and go as high as `(node_count - 1)`.
- `gear_factorization_level` (default value: `1`): This is the value up to which "wheel (and gear) factorization" are applied to "brute force." A value of `11` includes all prime factors of `11` and below and works well for `PRIME_PROVER`, though significantly higher might be preferred in certain cases.
- `wheel_factorization_level` (default value: `1`): "Wheel" vs. "gear" factorization balances two types of factorization wheel ("wheel" vs. "gear" design) that often work best when the "wheel" is only a few prime number levels lower than gear factorization. Optimized implementation for wheels is only available up to `13`. The primes above "wheel" level, up to "gear" level, are the primes used specifically for "gear" factorization. Wheel factorization is also applied to map the sieving intervale of `FACTOR_FINDER` mode onto non-multiples on the wheel, if the level is set above `1` (which might not actually pay dividends in practical complexity, but we leave it for your experimentation). In `FACTOR_FINDER` mode, wheel factorization multiples are systematically avoided by construction, while gear factorization multiples are just rejected after-the-fact, without special construction to avoid them. (It is possible in theory to implement handling for gear factorization by construction, though, and that might be added in a future release.)
- `sieving_bound_multiplier` (default value: `1.0`): This controls the sieving bound and is calibrated such that it linearly multiplies the square root of the number to factor (for each `1.0` increment). While this might be a huge bound, remember that sieving termination is primarily controlled by when `gaussian_elimination_row_multiplier` is exactly satisfied.
- `smoothness_bound_multiplier` (default value: `1.0`): This controls smoothness bound and is calibrated such that it linearliy multiplies `pow(exp(0.5 * sqrt(log(N) * log(log(N)))), sqrt(2.0)/4)` for `N` being the number to factor (for each `1.0` increment). This was a heuristic suggested by Elara (an OpenAI custom GPT).
- `gaussian_elimination_row_offset` (default value: `1`): This controls the number of rows greater than the count of smooth primes that are sieved before Gaussian elimination. Basically, for each increment starting with `1`, the chance of finding at least one solution in Gaussian elimination goes like `(1 - 2^(-m))` for a setting value of `m`: `1` value is a 50% chance of success, and the chance of failure is halved for each unit of `1` added. So long as this setting is appropriately low enough, `sieving_bound_multiplier` can be set basically arbitrarily high.
- `check_small_factors` (default value: `False`): `True` performs initial-phase trial division up to the smoothness bound, and `False` skips it.

All variables defaults can also be controlled by environment variables:
- `FINDAFACTOR_METHOD` (integer value)
- `FINDAFACTOR_NODE_COUNT`
- `FINDAFACTOR_NODE_ID`
- `FINDAFACTOR_GEAR_FACTORIZATION_LEVEL`
- `FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL`
- `FINDAFACTOR_SIEVING_BOUND_MULTIPLIER`
- `FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER`
- `FINDAFACTOR_GAUSSIAN_ELIMINATION_ROW_OFFSET`
- `FINDAFACTOR_CHECK_SMALL_FACTORS` (`True` if set at all, otherwise `False`)

## About 
This library was originally called ["Qimcifa"](https://github.com/vm6502q/qimcifa) and demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has since been developed into a general factoring algorithm and tool.

**Special thanks to OpenAI GPT "Elara," for help with indicated region of contributed code!**
