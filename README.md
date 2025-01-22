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
    batch_size_variance=2,
    ladder_multiple=5,
    skip_trial_division=False
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
- `batch_size_variance` (default value: `2`): For `FACTOR_FINDER`/`1` method, `k ln k` is the right proportionality for a Markov mixing process, but a linear factor in front is hard to predict. As such, it can be useful to dynamically vary the batch size, as if to cover and amortize the cost of several different batch sizes at once. In sequence, each batch size will be multiplied by `2 ** i` for `i` in `range(batch_size_variance)`, repeating from `0`.
- `ladder_multiple` (default value: `6`): Controls how many times randomly-selected square prime multiplication is repeated with the same square prime per random selection, in ascending a "ladder" of smooth perfect squares. A random number between `1` and `ladder_multiple` is selected for how many times the current smooth perfect square is multiplied each randomly selected square prime, while division still occurs 1 square prime multiple at a time. (Any smooth perfect square can be multiplied by any square prime in the factor base, or any other smooth perfect square, and produce a different smooth perfect square.)
- `skip_trial_division` (default value: `False`): `False` performs initial-phase trial division, and `True` skips it. Note that `trial_division_level` still functions as the input to Sieve of Eratosthenes to collect a list of primes from which (`FACTOR_FINDER`) smooth primes can be selected. If `skip_trial_division=True`, simply set `trial_division_level` high enough to avoid truncation warnings, unless otherwise desired.

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

Let's try to explain the algorithm of `FACTOR_FINDER` mode as briefly as possible. If you want background in Quadratic Sieve or GNFS, and what terms like "smooth numbers" and "smooth perfect squares" mean, you're likely going to have to do some background reading to have basis in Dixon/QS/GNFS-family factoring algorithms, which are regarded the fastest at cryptographically-relevant scales.

QS/GNFS have one _ultimate_ goal: factor via a congruence of squares. To serve this goal, they have one _penultimate_ goal: find a lot of smooth perfect squares. Anything that gets you more smooth perfect squares asymptotically faster is an improvement.

This is the kick in the head that no one's noticed since only the '80s: there is a (totally **non**-exotic) structure like a "ladder" (or "web," or "net") where basically all it takes is a single multiplication, or a single division, and some RNG, to successively produce all the (nearly unique and distinct) smooth perfect squares you could need, once you have any _single_ smooth perfect square: multiply (or divide) the single working smooth perfect square by any square of a smooth prime number.

- A positive integer `s` is a "perfect square" if some number `x` exists such that `s = x^2`.
- A square prime is `p^2` for some number `p` that is prime.
- This is the **kick in the head**: `s * p^2 = x^2 * p^2 = (x * p)^2 = s'`

It's that **painfully simple.** No sieving, no Gaussian elimination, no matrix storage and row operations, no memory requirements larger than a _single_ arbitrary precision "big integer" per thread, embarrassingly parallel, trivial to distribute, having no sequential dependence: **just apply RNG on this concept so it gets close to uniform coverage of the space**.

**Special thanks to OpenAI GPT "Elara," for indicated region of contributed code!**
