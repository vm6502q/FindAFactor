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
from FindAFactor import find_a_factor

to_factor = 1000

factor = find_a_factor(
    to_factor,
    use_congruence_of_squares=False,
    node_count=1, node_id=0,
    trial_division_level=2**20,
    gear_factorization_level=13,
    wheel_factorization_level=11,
    smoothness_bound_multiplier=6.35,
    batch_size_multiplier=1024.0
)
```

The `find_a_factor()` function should return any nontrivial factor of `to_factor` (that is, any factor besides `1` or `to_factor`) if it exists. If a nontrivial factor does _not_ exist (i.e., the number to factor is prime), the function will return `1` or the original `to_factor`.

- `use_congruence_of_squares` (default value: `False`): This attempts to check congruence of squares with "smooth" numbers, if `true`.
- `node_count` (default value: `1`): `FindAFactor` can perform factorization in a _distributed_ manner, across nodes, without network communication! When `node_count` is set higher than `1`, the search space for factors is segmented equally per node. If the number to factor is semiprime, and brute-force search is used instead of congruence of squares, for example, all nodes except the one that happens to contain the (unknown) prime factor less than the square root of `to_factor` will ultimately return `1`, while one node will find and return this factor. For best performance, every node involved in factorization should have roughly the same CPU throughput capacity.
- `node_id` (default value: `0`): This is the identifier of this node, when performing distributed factorization with `node_count` higher than `1`. `node_id` values start at `0` and go as high as `(node_count - 1)`.
- `trial_division_level` (default value: `2**20`): Trial division is carried out as a preliminary round for all primes up this number. If you need more primes for your smoothness bound, increase this level.
- `gear_factorization_level` (default value: `13`): This is the value up to which "wheel (and gear) factorization" and trial division are used to check factors and optimize "brute force," in general. The default value of `13` includes all prime factors of `13` and below and works well for `use_congruence_of_squares=True`, though significantly higher might be preferred in certain cases.
- `wheel_factorization_level` (default value: `11`): "Wheel" vs. "gear" factorization balances two types of factorization wheel ("wheel" vs. "gear" design) that often work best when the "wheel" is only a few prime number levels lower than gear factorization. Optimized implementation for wheels is only available up to `13`. The primes above "wheel" level, up to "gear" level, are the primes used specifically for "gear" factorization.
- `smoothness_bound_multiplier` (default value: `6.35`): starting with the first prime number after wheel factorization, the congruence of squares approach (with Quadratic Sieve) has a "smoothness bound" unit with as many distinct prime numbers as bits in the number to facto0r (for argument of `1.0` multiplier). To increase or decrease this number, consider it multiplied by the value of `smoothness_bound_multiplier`.
- `batch_size_multiplier` (default value: `1024.0`): Each `1.0` increment of the multiplier is 2 cycles of gear and wheel factorization, alternating every other cycle between bottom of guessing range and top of guessing range, for every thread in use.

All variables defaults can also be controlled by environment variables:
- `FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES` (any value makes `True`, while default is `False`)
- `FINDAFACTOR_NODE_COUNT`
- `FINDAFACTOR_NODE_ID`
- `FINDAFACTOR_TRIAL_DIVISION_LEVEL`
- `FINDAFACTOR_GEAR_FACTORIZATION_LEVEL`
- `FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL`
- `FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER`
- `FINDAFACTOR_BATCH_SIZE_MULTIPLIER`

## Factoring parameter strategy

The developer anticipates this single-function set of parameters, as API, is the absolutely most complicated FindAFactor likely ever needs to get. The _only_ required argument is `to_factor`, and it _just works,_ for any number that should reasonably take about less than a second to a few minutes. However, if you're running larger numbers for longer than that, of course, it's worth investing 15 to 30 minutes to read the explanation above in the `README` of every argument and play with the settings, a little. But, with these options, you have _total control to tune_ the algorithm in any all ways necessary to adapt to system resource footprint.

Advantage for `use_congruence_of_squares` is beyond the hardware scale of the developer's experiments, in practicality, but it can be shown to work correctly (at disadvantage, at small factoring bit-width scales). The anticipated use case is to turn this option on when approaching the size of modern-day RSA semiprimes in use.

If this is your use case, you want to specifically consider `smoothness_bound_multiplier`, `batch_size_multiplier`, and potentially `thread_count` for managing memory. By default, as many primes are kept for "smooth" number sieving as bits in the number to factor. This is multiplied by `smooth_bound_multiplier` (and cast to a discrete number of primes in total). This multiplier tends to predominate memory, but `batch_size_multiplier` can also cause problems if set too high or low, as a high value might exhaust memory, while a low value increases potentially nonproductive Gaussian elimination checks, which might be more efficient if batched higher. Our expectation is that most systems will benefit from significant experimentation with and fine tuning of `batch_size_multiplier` to something other than default, potentially to the point of using about half of available memory.

`wheel_factorization_level` and `gear_factorization_level` are common to both `use_congruence_of_squares` (i.e., Gaussian elimination for perfect squares) and "brute force." `11` for gear and `5` for wheel limit works well for small numbers. You'll definitely want to consider (gear/wheel) `13`/`7` or `17`/`11` (or even other values, maybe system-dependent) as your numbers to factor approach cryptographic relevance.

As for `node_count` and `node_id`, believe it or not, factoring parallelism can truly be that simple: just run different node IDs in the set on different (roughly homogenous) isolated CPUs, without networking. The only caveat is that you must manually detect when any single node has found an answer (which is trivial to verify) and manually interrupt other nodes still working (or leave them to complete on their own).

## About 
This library was originally called ["Qimcifa"](https://github.com/vm6502q/qimcifa) and demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has since been developed into a general factoring algorithm and tool.

`FindAFactor` uses heavily wheel-factorized brute-force "exhaust" numbers as "smooth" inputs to Quadratic Sieve, widely regarded as the asymptotically second fastest algorithm class known for cryptographically relevant semiprime factoring. Actually, the primary difference between Quadratic Sieve (QS, regarded second-fastest) and General Number Field Sieve (GNFS, fastest) is based in how "smooth" numbers are generated as intermediate inputs to Gaussian elimination, and the "brute-force exhaust" of Qimcifa provides smooth numbers rather than canonical polynomial generators for QS or GNFS, so whether `FindAFactor` is theoretically fastest depends on how good its smooth number generation is (which is an open question). `FindAFactor` is C++ based, with `pybind11`, which tends to make it faster than pure Python approaches. For the quick-and-dirty application of finding _any single_ nontrivial factor, something like at least 80% of positive integers will factorize in a fraction of a second, but the most interesting cases to consider are semiprime numbers, for which `FindAFactor` should be about as asymptotically competitive as similar Quadratic Sieve implementations.

Our original contribution to Quadratic Sieve seems to be wheel factorization to 13 or 17 and maybe the idea of using the "exhaust" of a brute-force search for smooth number inputs for Quadratic Sieve. For wheel factorization (or "gear factorization"), we collect a short list of the first primes and remove all of their multiples from a "brute-force" guessing range by mapping a dense contiguous integer set, to a set without these multiples, relying on both a traditional "wheel," up to a middle prime number (of `11`), and a "gear-box" that stores increment values per prime according to the principles of wheel factorization, but operating semi-independently, to reduce space of storing the full wheel.

Beyond this, we gain a functional advantage of a square-root over a more naive approach, by setting the brute force guessing range only between the highest prime in wheel factorization and the (modular) square root of the number to factor: if the number is semiprime, there is exactly one correct answer in this range, but including both factors in the range to search would cost us the square root advantage.

Factoring this way is surprisingly easy to distribute: basically 0 network communication is needed to coordinate an arbitrarily high amount of parallelism to factor a single number. Each brute-force trial division instance is effectively 100% independent of all others (i.e. entirely "embarrassingly parallel"), and these guesses can seed independent Gaussian elimination matrices, so `FindAFactor` offers an extremely simply interface that allows work to be split between an arbitrarily high number of nodes with absolutely no network communication at all. In terms of incentives of those running different, cooperating nodes in the context of this specific number of integer factoring, all one ultimately cares about is knowing the correct factorization answer _by any means._ For pratical applications, there is no point at all in factoring a number whose factors are already known. When a hypothetical answer is forwarded to the (0-communication) "network" of collaborating nodes, _it is trivial to check whether the answer is correct_ (such as by simply entering the multiplication and equality check with the original number into a Python shell console)! Hence, collaborating node operators only need to trust that all participants in the "network" are actually performing their alloted segment of guesses and would actually communicate the correct answer to the entire group of collaborating nodes if any specific invidual happened to find the answer, but any purported answer is still trivial to verify.

**Special thanks to OpenAI GPT "Elara," for indicated region of contributed code!**
