# FindAFactor
Find any nontrivial factor of a number

## Copyright and license
(c) Daniel Strano and the Qrack contributors 2017-2024. All rights reserved.

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

factor = find_a_factor(to_factor, use_congruence_of_squares=True, node_count=1, node_id=0, gear_factorization_level=13, wheel_factorization_level=11, smoothness_bound_multiplier=1.0)
```

The `find_a_factor()` function should return any nontrivial factor of `to_factor` (that is, any factor besides `1` or `to_factor`) if it exists. If a nontrivial factor does _not_ exist (i.e., the number to factor is prime), the function will return `1` or the original `to_factor`.

- `use_congruence_of_squares` (default value: `True`): This attempts to check congruence of squares with Gaussian elimination for Quadratic Sieve.
- `node_count` (default value: `1`): `FindAFactor` can perform factorization in a _distributed_ manner, across nodes, without network communication! When `node_count` is set higher than `1`, the search space for factors is segmented equally per node. If the number to factor is semiprime, and brute-force search is used instead of congruence of squares, for example, all nodes except the one that happens to contain the (unknown) prime factor less than the square root of `to_factor` will ultimately return `1`, while one node will find and return this factor. For best performance, every node involved in factorization should have roughly the same CPU throughput capacity.
- `node_id` (default value: `0`): This is the identifier of this node, when performing distributed factorization with `node_count` higher than `1`. `node_id` values start at `0` and go as high as `(node_count - 1)`.
- `gear_factorization_level` (default value: `13`): This is the value up to which "wheel (and gear) factorization" and trial division are used to check factors and optimize "brute force," in general. The default value of `13` includes all prime factors of `13` and below and works well in general, though `17` or higher might be preferred in certain cases.
- `wheel_factorization_level` (default value: `11`): "Wheel" vs. "gear" factorization balances two types of factorization wheel ("wheel" vs. "gear" design) that usually work best when the "wheel" is one prime smaller than the largest "gear." (In other words, as `11` and `13` are consecutive primes, "wheels" and "gears" tend to work best with consecutive limits.) Optimized implementation for wheels is only available up to 11.
- `smoothness_bound_multiplier` (default value: `1.0`): starting with the first prime number after wheel factorization, the congruence of squares approach (with Quadratic Sieve) takes a default "smoothness bound" with as many distinct prime numbers as bits in the number to factor (for default argument of `1.0` multiplier). To increase or decrease this number, consider it multiplied by the value of `smoothness_bound_multiplier`.

All variables defaults can also be controlled by environment variables:
- `FINDAFACTOR_USE_CONGRUENCE_OF_SQUARES`
- `FINDAFACTOR_NODE_COUNT`
- `FINDAFACTOR_NODE_ID`
- `FINDAFACTOR_GEAR_FACTORIZATION_LEVEL`
- `FINDAFACTOR_WHEEL_FACTORIZATION_LEVEL`
- `FINDAFACTOR_SMOOTHNESS_BOUND_MULTIPLIER`

## About 
This library was originally called ["Qimcifa"](https://github.com/vm6502q/qimcifa) and demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has since been developed into a general factoring algorithm and tool.

`FindAFactor` uses heavily wheel-factorized brute-force "exhaust" numbers as "smooth" inputs to Quadratic Sieve, widely regarded as the asymptotically second fastest algorithm class known for cryptographically relevant semiprime factoring. Actually, the primary difference between Quadratic Sieve (QS, regarded second-fastest) and General Number Field Sieve (GNFS, fastest) is based in how "smooth" numbers are generated as intermediate inputs to Gaussian elimination, and the "brute-force exhaust" of Qimcifa provides smooth numbers rather than canonical polynomial generators for QS or GNFS, so whether `FindAFactor` is theoretically fastest depends on how good its smooth number generation is (which is an open question). `FindAFactor` is C++ based, with `pybind11`, which tends to make it faster than pure Python approaches. For the quick-and-dirty application of finding _any single_ nontrivial factor, something like at least 80% of positive integers will factorize in a fraction of a second, but the most interesting cases to consider are semiprime numbers, for which `FindAFactor` should be about as asymptotically competitive as similar Quadratic Sieve implementations.

Our original contribution to Quadratic Sieve seems to be wheel factorization to 13 or 17 and maybe the idea of using the "exhaust" of a brute-force search for smooth number inputs for Quadratic Sieve. For wheel factorization (or "gear factorization"), we collect a short list of the first primes and remove all of their multiples from a "brute-force" guessing range by mapping a dense contiguous integer set, to a set without these multiples, relying on both a traditional "wheel," up to a middle prime number (of `11`), and a "gear-box" that stores increment values per prime according to the principles of wheel factorization, but operating semi-independently, to reduce space of storing the full wheel.

Beyond this, we gain a functional advantage of a square-root over a more naive approach, by setting the brute force guessing range only between the highest prime in wheel factorization and the (modular) square root of the number to factor: if the number is semiprime, there is exactly one correct answer in this range, but including both factors in the range to search would cost us the square root advantage.

Factoring this way is surprisingly easy to distribute: basically 0 network communication is needed to coordinate an arbitrarily high amount of parallelism to factor a single number. Each brute-force trial division instance is effectively 100% independent of all others (i.e. entirely "embarrassingly parallel"), and these guesses can seed independent Gaussian elimination matrices, so `FindAFactor` offers an extremely simply interface that allows work to be split between an arbitrarily high number of nodes with absolutely no network communication at all. In terms of incentives of those running different, cooperating nodes in the context of this specific number of integer factoring, all one ultimately cares about is knowing the correct factorization answer _by any means._ For pratical applications, there is no point at all in factoring a number whose factors are already known. When a hypothetical answer is forwarded to the (0-communication) "network" of collaborating nodes, _it is trivial to check whether the answer is correct_ (such as by simply entering the multiplication and equality check with the original number into a Python shell console)! Hence, collaborating node operators only need to trust that all participants in the "network" are actually performing their alloted segment of guesses and would actually communicate the correct answer to the entire group of collaborating nodes if any specific invidual happened to find the answer, but any purported answer is still trivial to verify.

**Special thanks to OpenAI GPT "Elara," for indicated region of contributed code!**
