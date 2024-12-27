# FindAFactor
Find any nontrivial factor of a number

## Copyright and license
(c) Daniel Strano and the Qrack contributors 2017-2024. All rights reserved.

## Usage

```py
from FindAFactor import find_a_factor

to_factor = 1000

factor = find_a_factor(to_factor, node_count=1, node_id=0, wheel_factorization_level=17, use_congruence_of_squares=False)
```

The `find_a_factor()` function should return any nontrivial factor of `to_factor` (that is, any factor besides `1` or `to_factor`) if it exists. If a nontrivial factor does _not_ exist (i.e., the number to factor is prime), the function will return `1` or the original `to_factor`. We do not (yet) guarantee that `find_a_factor()` can be used 100% reliably for primality proving (if it returns `1`), but, by design intention, it ultimately should. (It's possible that edge cases might occassionally be missed, but this ultimately shouldn't happen, as we improve the library.)

- `node_count` (default value: `1`): `FindAFactor` can perform factorization in a _distributed_ manner, across nodes, without network communication! When `node_count` is set higher than `1`, the search space for factors is segmented equally per node. If the number to factor is semiprime, for example, all nodes except the one that happens to contain the (unknown) prime factor less than the square root of `to_factor` will ultimately return `1`, while one node will find and return this factor. For best performance, every node involved in factorization should have roughly the same CPU throughput capacity.
- `node_id` (default value: `0`): This is the identifier of this node, when performing distributed factorization with `node_count` higher than `1`. `node_id` values start at `0` and go as high as `(node_count - 1)`.
- `wheel_factorization_level` (default value: `17`): This is the value up to which "wheel factorization" and trial division are used to check factors and optimize "brute force," in general. The default value of `17` includes all prime factors of `17` and below and works well in general, though `19` or higher might be slightly preferred in certain cases.
- `use_congruence_of_squares` (default value: `False`): This attempts to check congruence of squares. (This mode will ultimately become Gaussian elimination, with future development.)

## About 
This library was originally called ["Qimcifa"](https://github.com/vm6502q/qimcifa) and demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has since been developed into a general factoring algorithm and tool.

Admittedly, `FindAFactor` is _not yet_ anything particularly "groundbreaking" for factorization algorithms. Its biggest advantage over certain other similar software is basically only that it is C++ based, with `pybind11`, which tends to make it faster than pure Python approaches, sometimes. For the quick-and-dirty application of finding _any single_ nontrivial factor, something like at least 80% of positive integers will factorize in a fraction of a second, but the most interesting cases to consider are semiprime numbers, for which `FindAFactor` is not yet particularly competitive, maybe at all. Before the v1.0.0 release, we hope to take the potentially novel "tricks" `FindAFactor` uses and incoporate them with an industry-standard factorization algorithm like Quadratic Sieve or General Number Field Sieve (GNFS). When this is done, `FindAFactor` might ultimately constitute a "competitive" approach for semiprime number factoring, so stay tuned for the v1.0.0 release.

The only potentially "original" part of this factoring algorithm is the "reverse wheel factorization," as far as I can tell. The idea is, instead of performing typical wheel factorization or trial division, we collect a short list of the first primes and remove all of their multiples from a "brute-force" guessing range by mapping a dense contiguous integer set, to a set without these multiples, by successively applying `guess = guess + guess / (p[i] - 1U) + 1U` for prime "`p`" in ascending (or any) order. Each prime applied this way effectively multiplies the brute-force guessing cardinality by a fraction (p-1)/p. Whatever "level" of primes we use, the cost per "guess" becomes higher.

Then, we have a tuner that empirically estimates the cost per guess, and we multiply this by the (known) total cardinality of potential guesses. Whichever reverse wheel factorization level has the lowest product of average cost per guess times guessing set cardinality should have the best performance, and the best level increases with the scale of the problem.

Beyond this, we gain a functional advantage of a square-root over a more naive approach, by setting the brute force guessing range only between the highest prime in reverse wheel factorization and the (modular) square root of the number to factor: if the number is semiprime, there is exactly one correct answer in this range, but including both factors in the range to search would cost us the square root advantage.

Beyond that, we observed that many simple and well-known factoring techniques just don't pay dividends, for semiprime factoring. There's basically no point in checking either congruence of squares or even for a greatest common divisor, as these techniques require some dynamically-variable overhead, and it tends to be faster (for semiprimes) just to check if a guess is an exact factor, on the smallest range we can identify that contains at least one correct answer.

So, this is actually quite rudimentary and just "brute force," except for "reverse wheel factorization" and the upper bound on the guessing range. It just better work entirely in CPU cache, then, but it only requires de minimis maximum memory footprint. (There are congruence of squares and greatest common divisor checks available for numbers besides semiprimes.)

Theoretically, this algorithm might return to its original "quantum-inspired" design with the availability of a high-quality, high-throughput generator of uniform random bit strings. If we were to use the algorithm as-is, except guessing according to a uniform random distribution instead of systematically ascending through every possible "guess," then the average time to solution can be realized in any case, unlike the deterministic version of the algorithm. Then, no work towards the solution can ever be lost in event of interruption of the program, because every single guess (even the first) has the same probability (in the ideal) of leading to successful factoring.

A nearly "brute-force" technique like this has a surprising advantage: basically 0 network communication is needed to coordinate an arbitrarily high amount of parallelism to factor a single number. Each trial division instance is effectively 100% independent of all others (i.e. entirely "embarrassingly parallel"), so Qimcifa offers an interface that allows work to be split between an arbitrarily high number of nodes with absolutely no network communication at all. In terms of incentives of those running different, cooperating nodes in the context of this specific number of integer factoring, all one ultimately cares about is knowing the correct factorization answer _by any means._ For pratical applications, there is no point at all in factoring a number whose factors are already known. When a hypothetical answer is forwarded to the (0-communication) "network" of collaborating nodes, _it is trivial to check whether the answer is correct_ (such as by simply entering the multiplication and equality check with the original number into a Python shell console)! Hence, collaborating node operators only need to trust that all participants in the "network" are actually performing their alloted segment of guesses and would actually communicate the correct answer to the entire group of collaborating nodes if any specific invidual happened to find the answer, but any purported answer is still trivial to verify.
