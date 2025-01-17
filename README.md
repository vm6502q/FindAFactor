# FindAFactor
Find any nontrivial factor of a number

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
    trial_division_level=2**20,
    gear_factorization_level=11,
    wheel_factorization_level=11,
    smoothness_bound_multiplier=1.0,
    batch_size_multiplier=512.0,
    ladder_multiple=4
)
```

The `find_a_factor()` function should return any nontrivial factor of `to_factor` (that is, any factor besides `1` or `to_factor`) if it exists. If a nontrivial factor does _not_ exist (i.e., the number to factor is prime), the function will return `1` or the original `to_factor`.

- `method` (default value: `PRIME_SOLVER`/`0`): `PRIME_SOLVER`/`0` will prove that a number is prime (by failing to find any factors with wheel and gear factorization). `FACTOR_FINDER`/`1` is optimized for the assumption that the number has at least two nontrivial factors.
- `node_count` (default value: `1`): `FindAFactor` can perform factorization in a _distributed_ manner, across nodes, without network communication! When `node_count` is set higher than `1`, the search space for factors is segmented equally per node. If the number to factor is semiprime, and brute-force search is used instead of congruence of squares, for example, all nodes except the one that happens to contain the (unknown) prime factor less than the square root of `to_factor` will ultimately return `1`, while one node will find and return this factor. For best performance, every node involved in factorization should have roughly the same CPU throughput capacity.
- `node_id` (default value: `0`): This is the identifier of this node, when performing distributed factorization with `node_count` higher than `1`. `node_id` values start at `0` and go as high as `(node_count - 1)`.
- `trial_division_level` (default value: `2**20`): Trial division is carried out as a preliminary round for all primes up this number. If you need more primes for your smoothness bound, increase this level.
- `gear_factorization_level` (default value: `11`): This is the value up to which "wheel (and gear) factorization" and trial division are used to check factors and optimize "brute force," in general. The default value of `11` includes all prime factors of `11` and below and works well in general, though significantly higher might be preferred in certain cases.
- `wheel_factorization_level` (default value: `11`): "Wheel" vs. "gear" factorization balances two types of factorization wheel ("wheel" vs. "gear" design) that often work best when the "wheel" is only a few prime number levels lower than gear factorization. Optimized implementation for wheels is only available up to `13`. The primes above "wheel" level, up to "gear" level, are the primes used specifically for "gear" factorization.
- `smoothness_bound_multiplier` (default value: `1.0`): starting with the first prime number after wheel factorization, the congruence of squares approach (with Quadratic Sieve) has a "smoothness bound" unit with as many distinct prime numbers as bits in the number to factor (for argument of `1.0` multiplier). To increase or decrease this number, consider it multiplied by the value of `smoothness_bound_multiplier`.
- `batch_size_multiplier` (default value: `512.0`): For `FACTOR_FINDER`/`1` method, each `1.0` increment of the multiplier adds `k ln k` Markov mixing replacement steps for `k` count of smooth primes, before reseeding Monte Carlo.
- `batch_size_variance` (default value: `4`): For `FACTOR_FINDER`/`1` method, `k ln k` is the right proportionality for a Markov mixing process, but a linear factor in front is hard to predict. As such, it can be useful to dynamically vary the batch size, as if to cover and amortize the cost of several different batch sizes at once. In sequence, each batch size will be multiplied by `2 ** i` for `i` in `range(batch_size_variance)`, repeating from `0`.
- `ladder_multiple` (default value: `4`): Controls how many times randomly-selected square prime multiplication is repeated with the same square prime per random selection, in ascending a "ladder" of smooth perfect squares. (Any smooth perfect square can be multiplied by any square prime in the factor base, or any other smooth perfect square, and produce a different smooth perfect square.)

All variables defaults can also be controlled by environment variables:
- `FINDAFACTOR_METHOD` (integer value)
- `FINDAFACTOR_NODE_COUNT`
- `FINDAFACTOR_NODE_ID`
- `FINDAFACTOR_TRIAL_DIVISION_LEVEL`
- `FINDAFACTOR_GEAR_FACTORIZATION_LEVEL`
- `FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL`
- `FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER`
- `FINDAFACTOR_BATCH_SIZE_MULTIPLIER`
- `FINDAFACTOR_BATCH_SIZE_VARIANCE`

## About 
This library was originally called ["Qimcifa"](https://github.com/vm6502q/qimcifa) and demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has since been developed into a general factoring algorithm and tool.

`FindAFactor` uses heavily wheel-factorized brute-force "exhaust" numbers as "smooth" inputs to Quadratic Sieve, widely regarded as the asymptotically second fastest algorithm class known for cryptographically relevant semiprime factoring. Actually, the primary difference between Quadratic Sieve (QS, regarded second-fastest) and General Number Field Sieve (GNFS, fastest) is based in how "smooth" numbers are generated as intermediate inputs to Gaussian elimination, and the "brute-force exhaust" of Qimcifa provides smooth numbers rather than canonical polynomial generators for QS or GNFS, so whether `FindAFactor` is theoretically fastest depends on how good its smooth number generation is (which is an open question). `FindAFactor` is C++ based, with `pybind11`, which tends to make it faster than pure Python approaches. For the quick-and-dirty application of finding _any single_ nontrivial factor, something like at least 80% of positive integers will factorize in a fraction of a second, but the most interesting cases to consider are semiprime numbers, for which `FindAFactor` should be about as asymptotically competitive as similar Quadratic Sieve implementations.

Our original contribution to Quadratic Sieve seems to be wheel factorization to 13 or 17 and maybe the idea of using the "exhaust" of a brute-force search for smooth number inputs for Quadratic Sieve. For wheel factorization (or "gear factorization"), we collect a short list of the first primes and remove all of their multiples from a "brute-force" guessing range by mapping a dense contiguous integer set, to a set without these multiples, relying on both a traditional "wheel," up to a middle prime number (of `11`), and a "gear-box" that stores increment values per prime according to the principles of wheel factorization, but operating semi-independently, to reduce space of storing the full wheel.

Beyond this, we gain a functional advantage of a square-root over a more naive approach, by setting the brute force guessing range only between the highest prime in wheel factorization and the (modular) square root of the number to factor: if the number is semiprime, there is exactly one correct answer in this range, but including both factors in the range to search would cost us the square root advantage.

Factoring this way is surprisingly easy to distribute: basically 0 network communication is needed to coordinate an arbitrarily high amount of parallelism to factor a single number. Each brute-force trial division instance is effectively 100% independent of all others (i.e. entirely "embarrassingly parallel"), and these guesses can seed independent Gaussian elimination matrices, so `FindAFactor` offers an extremely simply interface that allows work to be split between an arbitrarily high number of nodes with absolutely no network communication at all. In terms of incentives of those running different, cooperating nodes in the context of this specific number of integer factoring, all one ultimately cares about is knowing the correct factorization answer _by any means._ For pratical applications, there is no point at all in factoring a number whose factors are already known. When a hypothetical answer is forwarded to the (0-communication) "network" of collaborating nodes, _it is trivial to check whether the answer is correct_ (such as by simply entering the multiplication and equality check with the original number into a Python shell console)! Hence, collaborating node operators only need to trust that all participants in the "network" are actually performing their alloted segment of guesses and would actually communicate the correct answer to the entire group of collaborating nodes if any specific invidual happened to find the answer, but any purported answer is still trivial to verify.

**Special thanks to OpenAI GPT "Elara," for indicated region of contributed code!**
