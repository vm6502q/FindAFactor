///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2025. All rights reserved.
//
// "A quantum-inspired Monte Carlo integer factoring algorithm"
//
// Special thanks to https://github.com/NachiketUN/Quadratic-Sieve-Algorithm, for providing an example implementation of Quadratic sieve.
//
//**Special thanks to OpenAI GPT "Elara," for help with indicated region of contributed code!**
//
// Licensed under the MIT License.
// See LICENSE.md in the project root or
// https://opensource.org/license/mit for details.

#include "dispatchqueue.hpp"
#include "wheel_factorization.hpp"

#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdlib.h>
#include <string>

#include <boost/dynamic_bitset.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Qimcifa {

typedef boost::multiprecision::cpp_int BigInteger;

const unsigned CpuCount = std::thread::hardware_concurrency();
DispatchQueue dispatch(CpuCount);

size_t biggestWheel = 1U;
std::vector<size_t> wheel;

BigInteger smoothForwardFn(const BigInteger &p) {
  return wheel[(size_t)(p % wheel.size())] + (p / wheel.size()) * biggestWheel;
}
BigInteger smoothBackwardFn(const BigInteger &p) {
  return std::distance(wheel.begin(), std::lower_bound(wheel.begin(), wheel.end(), (size_t)(p % biggestWheel))) + wheel.size() * (p / biggestWheel) + 1U;
}


// See https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
BigInteger ipow(BigInteger base, size_t exp) {
  BigInteger result = 1U;
  for (;;) {
    if (exp & 1U) {
      result *= base;
    }
    exp >>= 1U;
    if (!exp) {
      break;
    }
    base *= base;
  }

  return result;
}

inline size_t log2(BigInteger n) {
  size_t pow = 0U;
  while (n >>= 1U) {
    ++pow;
  }
  return pow;
}

inline BigInteger gcd(const BigInteger& n1, const BigInteger& n2) {
  if (!n2) {
    return n1;
  }
  return gcd(n2, n1 % n2);
}

BigInteger sqrt(const BigInteger &toTest) {
  BigInteger start = 1U, end = toTest >> 1U, ans = 0U;
  do {
    const BigInteger mid = (start + end) >> 1U;

    // If toTest is a perfect square
    const BigInteger sqr = mid * mid;
    if (sqr == toTest) {
      return mid;
    }

    if (sqr < toTest) {
      // Since we need floor, we update answer when mid*mid is smaller than p, and move closer to sqrt(p).
      start = mid + 1U;
      ans = mid;
    } else {
      // If mid*mid is greater than p
      end = mid - 1U;
    }
  } while (start <= end);

  return ans;
}

size_t _sqrt(const size_t &toTest) {
  size_t start = 1U, end = toTest >> 1U, ans = 0U;
  do {
    const size_t mid = (start + end) >> 1U;

    // If toTest is a perfect square
    const size_t sqr = mid * mid;
    if (sqr == toTest) {
      return mid;
    }

    if (sqr < toTest) {
      // Since we need floor, we update answer when mid*mid is smaller than p, and move closer to sqrt(p).
      start = mid + 1U;
      ans = mid;
    } else {
      // If mid*mid is greater than p
      end = mid - 1U;
    }
  } while (start <= end);

  return ans;
}

inline size_t GetWheel5and7Increment(unsigned short &wheel5, unsigned long long &wheel7) {
  constexpr unsigned short wheel5Back = 1U << 9U;
  constexpr unsigned long long wheel7Back = 1ULL << 55U;
  size_t wheelIncrement = 0U;
  bool is_wheel_multiple = false;
  do {
    is_wheel_multiple = (bool)(wheel5 & 1U);
    wheel5 >>= 1U;
    if (is_wheel_multiple) {
      wheel5 |= wheel5Back;
      ++wheelIncrement;
      continue;
    }

    is_wheel_multiple = (bool)(wheel7 & 1U);
    wheel7 >>= 1U;
    if (is_wheel_multiple) {
      wheel7 |= wheel7Back;
    }
    ++wheelIncrement;
  } while (is_wheel_multiple);

  return wheelIncrement;
}

std::vector<size_t> SieveOfEratosthenes(const size_t &n) {
  std::vector<size_t> knownPrimes = {2U, 3U, 5U, 7U};
  if (n < 2U) {
    return std::vector<size_t>();
  }

  if (n < (knownPrimes.back() + 2U)) {
    const auto highestPrimeIt = std::upper_bound(knownPrimes.begin(), knownPrimes.end(), n);
    return std::vector<size_t>(knownPrimes.begin(), highestPrimeIt);
  }

  knownPrimes.reserve((size_t)(((double)n) / log((double)n)));

  // We are excluding multiples of the first few
  // small primes from outset. For multiples of
  // 2, 3, and 5 this reduces complexity to 4/15.
  const size_t cardinality = backward5(n);

  // Create a boolean array "prime[0..cardinality]"
  // and initialize all entries it as true. Rather,
  // reverse the true/false meaning, so we can use
  // default initialization. A value in notPrime[i]
  // will finally be false only if i is a prime.
  std::unique_ptr<bool[]> uNotPrime(new bool[cardinality + 1U]());
  bool *notPrime = uNotPrime.get();

  // Get the remaining prime numbers.
  // These wheel initializations are simply correct and optimal.
  // The integral form is rather a distinguishable bit set form.
  unsigned short wheel5 = 129U;
  unsigned long long wheel7 = 9009416540524545ULL;
  size_t o = 1U;
  for (;;) {
    o += GetWheel5and7Increment(wheel5, wheel7);

    const size_t p = forward3(o);
    if ((p * p) > n) {
      break;
    }

    if (notPrime[backward5(p)]) {
      continue;
    }

    knownPrimes.push_back(p);

    // We are skipping multiples of 2, 3, and 5
    // for space complexity, for 4/15 the bits.
    // More are skipped by the wheel for time.
    const size_t p2 = p << 1U;
    const size_t p4 = p << 2U;
    size_t i = p * p;

    // "p" already definitely not a multiple of 3.
    // Its remainder when divided by 3 can be 1 or 2.
    // If it is 2, we can do a "half iteration" of the
    // loop that would handle remainder of 1, and then
    // we can proceed with the 1 remainder loop.
    // This saves 2/3 of updates (or modulo).
    if ((p % 3U) == 2U) {
      notPrime[backward5(i)] = true;
      i += p2;
      if (i > n) {
        continue;
      }
    }

    for (;;) {
      if (i % 5U) {
        notPrime[backward5(i)] = true;
      }
      i += p4;
      if (i > n) {
        break;
      }

      if (i % 5U) {
        notPrime[backward5(i)] = true;
      }
      i += p2;
      if (i > n) {
        break;
      }
    }
  }

  for (;;) {
    const size_t p = forward3(o);
    if (p > n) {
      break;
    }

    o += GetWheel5and7Increment(wheel5, wheel7);

    if (notPrime[backward5(p)]) {
      continue;
    }

    knownPrimes.push_back(p);
  }

  return knownPrimes;
}

bool isMultiple(const BigInteger &p, const std::vector<size_t> &knownPrimes) {
  for (const size_t &prime : knownPrimes) {
    if (!(p % prime)) {
      return true;
    }
  }

  return false;
}

boost::dynamic_bitset<size_t> nestGearGeneration(std::vector<size_t> primes) {
  BigInteger radius = 1U;
  for (const size_t &i : primes) {
    radius *= i;
  }
  const size_t prime = primes.back();
  primes.pop_back();
  boost::dynamic_bitset<size_t> o;
  for (BigInteger i = 1U; i <= radius; ++i) {
    if (!isMultiple(i, primes)) {
      o.push_back(!(i % prime));
    }
  }
  o >>= 1U;

  return o;
}

std::vector<boost::dynamic_bitset<size_t>> generateGears(const std::vector<size_t> &primes) {
  std::vector<boost::dynamic_bitset<size_t>> output;
  std::vector<size_t> wheelPrimes;
  for (const size_t &p : primes) {
    wheelPrimes.push_back(p);
    output.push_back(nestGearGeneration(wheelPrimes));
  }

  return output;
}

size_t GetGearIncrement(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
  size_t wheelIncrement = 0U;
  bool is_wheel_multiple = false;
  do {
    for (size_t i = 0U; i < inc_seqs->size(); ++i) {
      boost::dynamic_bitset<size_t> &wheel = (*inc_seqs)[i];
      is_wheel_multiple = wheel.test(0U);
      wheel >>= 1U;
      if (is_wheel_multiple) {
        wheel.set(wheel.size() - 1U);
        break;
      }
    }
    ++wheelIncrement;
  } while (is_wheel_multiple);

  return wheelIncrement;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WRITTEN WITH HELP FROM ELARA (GPT) BELOW                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This was taken basically whole-cloth from Elara, with thanks.
BigInteger mod_exp(BigInteger base, BigInteger exp, BigInteger mod) {
  BigInteger result = 1U;
  base = base % mod;
  while (exp) {
    // If exp is odd, multiply base with result
    if (exp & 1U) {
      result = (result * base) % mod;
    }
    // Divide by 2
    exp = exp >> 1U;
    base = (base * base) % mod;
  }
  return result;
}

// Function to compute the Legendre symbol (N / p)
int legendreSymbol(BigInteger N, size_t p) {
  BigInteger result = mod_exp(N, (p - 1U) >> 1U, p);

  if (result == 0U) {
    return 0;  // N is divisible by p
  }

  if (result == 1U) {
    return 1;  // N is a quadratic residue mod p
  }

  return -1; // N is a non-quadratic residue mod p
}

// Function to generate factor base
std::vector<size_t> selectFactorBase(const BigInteger N, const std::vector<size_t>& primes) {
  std::vector<size_t> factorBase;
  for (size_t p : primes) {
    // Select only primes where (N/p) = 1
    if (legendreSymbol(N, p) == 1) {
      factorBase.push_back(p);
    }
  }
  return factorBase;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WRITTEN WITH HELP FROM ELARA (GPT) ABOVE                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Factorizer {
  std::mutex batchMutex;
  BigInteger toFactor;
  BigInteger toFactorSqrt;
  BigInteger qsBackwardLowBound;
  BigInteger batchRange;
  BigInteger batchNumber;
  BigInteger batchOffset;
  BigInteger batchTotal;
  BigInteger smoothWheelRadius;
  size_t wheelEntryCount;
  size_t rowLimit;
  bool isIncomplete;
  std::vector<size_t> smoothPrimes;
  std::vector<BigInteger> smoothNumberKeys;
  std::vector<boost::dynamic_bitset<size_t>> smoothNumberValues;
  ForwardFn forwardFn;
  ForwardFn backwardFn;

  Factorizer(const BigInteger &tf, const BigInteger &tfsqrt, const BigInteger &lb, const BigInteger &range, size_t nodeCount, size_t nodeId, size_t w, size_t rl, const BigInteger& bn,
             const std::vector<size_t> &sp, ForwardFn ffn, ForwardFn bfn)
    : toFactor(tf), toFactorSqrt(tfsqrt), qsBackwardLowBound(lb), batchRange(range), batchNumber(bn), batchOffset(nodeId * range), batchTotal(nodeCount * range),
    smoothWheelRadius(1U), wheelEntryCount(w), rowLimit(rl), isIncomplete(true), smoothPrimes(sp), forwardFn(ffn), backwardFn(bfn)
  {
    smoothNumberKeys.reserve(rowLimit);
    smoothNumberValues.reserve(rowLimit);
    for (const size_t p : smoothPrimes) {
      smoothWheelRadius *= p;
    }
  }

  BigInteger getNextBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);

    if (batchNumber >= batchRange) {
      isIncomplete = false;
    }

    return batchOffset + batchNumber++;
  }

  BigInteger getNextAltBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);

    if (batchNumber >= batchRange) {
      isIncomplete = false;
    }

    const BigInteger halfIndex = batchOffset + (batchNumber++ >> 1U);

    return ((batchNumber & 1U) ? batchTotal - (halfIndex + 1U) : halfIndex);
  }

  BigInteger bruteForce(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
    // Up to wheel factorization, try all batches up to the square root of toFactor.
    for (BigInteger batchNum = getNextAltBatch(); isIncomplete; batchNum = getNextAltBatch()) {
      const BigInteger batchStart = batchNum * wheelEntryCount;
      for (size_t batchItem = 1U; batchItem <= wheelEntryCount;) {
        const BigInteger n = forwardFn(batchStart + batchItem);
        if (!(toFactor % n) && (n != 1U) && (n != toFactor)) {
          isIncomplete = false;
          return n;
        }
        batchItem += GetGearIncrement(inc_seqs);
      }
    }

    return 1U;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                              WRITTEN WITH HELP FROM ELARA (GPT) BELOW                                  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Sieving function
  BigInteger sievePolynomials(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
    for (BigInteger batchNum = getNextBatch(); isIncomplete; batchNum = getNextBatch()) {
      // NOTE: If you want to add gear factorization back in, realize that these bounds
      // do not yet properly align to exact wheel boundaries, for full repetitions.
      // (They cycle through every validate candidate, but potentially with an offset.)
      const BigInteger batchStart = batchNum * wheelEntryCount + qsBackwardLowBound;
      for (size_t batchItem = 0U; batchItem < wheelEntryCount; ++batchItem) {
        // Make the candidate NOT a multiple on the wheels.
        const BigInteger x = forwardFn(batchStart + batchItem);
        // Make the candidate a perfect square.
        // The residue (mod N) needs to be smooth (but not a perfect square).
        // The candidate is guaranteed to be between toFactor and its square,
        // so subtracting toFactor is equivalent to % toFactor.
        const BigInteger ySqr = (x * x) - toFactor;
        const boost::dynamic_bitset<size_t> rfv = factorizationParityVector(ySqr);
        if (rfv.empty()) {
          // The number is useless to us.
          // batchItem += GetGearIncrement(inc_seqs);
          continue;
        }
        // We have a successful candidate.

        // If the candidate is already a perfect square,
        // we got lucky, and we might be done already.
        if (rfv.none()) {
          // x^2 % toFactor = y^2
          const BigInteger y = sqrt(ySqr);

          // Check x + y
          BigInteger factor = gcd(toFactor, x + y);
          if ((factor > 1U) && (factor < toFactor)) {
            isIncomplete = false;

            return factor;
          }

          // Avoid division by 0
          if (x != y) {
            // Check x - y
            factor = gcd(toFactor, x - y);
            if ((factor > 1U) && (factor < toFactor)) {
              isIncomplete = false;

              return factor;
            }
          }
        }

        std::lock_guard<std::mutex> lock(batchMutex);

        const auto& snvIt = std::find(smoothNumberValues.begin(), smoothNumberValues.end(), rfv);

        if (snvIt == smoothNumberValues.end()) {
          // This is a unique factorization parity row.
          smoothNumberValues.push_back(rfv);
          smoothNumberKeys.push_back(x);
          // If we have enough rows for Gaussian elimination already,
          // there's no reason to sieve any further.
          if (smoothNumberKeys.size() > rowLimit) {
            isIncomplete = false;

            return 1U;
          }
        } else {
          // Don't add this duplicate row, but check the square residue.
          // x^2 % toFactor = y^2
          const BigInteger _x = x * smoothNumberKeys[std::distance(smoothNumberValues.begin(), snvIt)];
          const BigInteger y = sqrt((_x * _x) % toFactor);

          // Check x + y
          BigInteger factor = gcd(toFactor, _x + y);
          if ((factor > 1U) && (factor < toFactor)) {
            isIncomplete = false;

            return factor;
          }

          // Avoid division by 0
          if (_x != y) {
            // Check x - y
            factor = gcd(toFactor, _x - y);
            if ((factor > 1U) && (factor < toFactor)) {
              isIncomplete = false;

              return factor;
            }
          }
        }

        // We must manually increment on exiting the loop body.
        // batchItem += GetGearIncrement(inc_seqs);
      }
    }

    return 1U;
  }

  std::vector<std::vector<size_t>> extractSolutionRows(const boost::dynamic_bitset<size_t>& marks) {
    std::vector<std::vector<size_t>> solutions;

    for (size_t col = 0U; col < marks.size(); ++col) {
      if (marks.test(col)) {
        // Skip pivot columns
        continue;
      }

      std::vector<size_t> selectedRows;
      boost::dynamic_bitset<size_t> solutionRow(marks.size(), 0U);

      // Collect rows that have a 1 in this free column
      for (size_t row = 0U; row < smoothNumberValues.size(); ++row) {
        if (smoothNumberValues[row].test(col)) {
          selectedRows.push_back(row);
           // XOR to construct dependency
          solutionRow ^= smoothNumberValues[row];
        }
      }

      // Ensure the dependency is valid (all exponents must sum to even parity)
      if (solutionRow.none()) {
        solutions.push_back(selectedRows);
      }
    }

    return solutions;
  }

  // Perform Gaussian elimination on a binary matrix
  std::vector<std::vector<size_t>> gaussianElimination() {
    const size_t rows = smoothNumberValues.size();
    boost::dynamic_bitset<size_t> marks(smoothPrimes.size(), 0U);
    for (size_t col = 0U; col < smoothPrimes.size(); ++col) {
      // Look for a pivot row in this column
      size_t row = col;
      for (; row < rows; ++row) {
        if (smoothNumberValues[row][col]) {
          // Make sure the rows are in reduced row echelon order.
          if (row != col) {
            std::swap(smoothNumberKeys[row], smoothNumberKeys[col]);
            std::swap(smoothNumberValues[row], smoothNumberValues[col]);
          }

          // Mark this column as having a pivot.
          marks.set(col);
          break;
        }
      }

      if ((col < smoothPrimes.size()) && marks[col]) {
        // Row might have been swapped.
        const boost::dynamic_bitset<size_t> &cm = smoothNumberValues[col];
        // Pivot found, now eliminate entries in this column
        const size_t maxLcv = std::min((size_t)CpuCount, rows);
        for (size_t cpu = 0U; cpu < maxLcv; ++cpu) {
          dispatch.dispatch([this, cpu, &col, &rows, &cm]() -> bool {
            // Notice that each thread updates rows with space increments of cpuCount,
            // based on the same unchanged outer-loop row, and this covers the inner-loop set.
            // We're covering every row except for the one corresponding to "col."
            // We break this into two loops to avoid an inner conditional to check whether row == col.
            const size_t midRow = std::min(col, rows);
            size_t irow = cpu;
            for (; irow < midRow; irow += CpuCount) {
              boost::dynamic_bitset<size_t> &rm = this->smoothNumberValues[irow];
              if (rm.test(col)) {
                // XOR-ing factorization rows
                // is like multiplying the numbers.
                rm ^= cm;
              }
            }
            if (irow == col) {
              irow += CpuCount;
            }
            for (; irow < rows; irow += CpuCount) {
              boost::dynamic_bitset<size_t> &rm = this->smoothNumberValues[irow];
              if (rm.test(col)) {
                // XOR-ing factorization rows
                // is like multiplying the numbers.
                rm ^= cm;
              }
            }
            return false;
          });
        }
        // All dispatched work must complete.
        dispatch.finish();
      }
    }

    const std::vector<std::vector<size_t>> solutions = extractSolutionRows(marks);

    if (solutions.empty()) {
      throw std::runtime_error("Gaussian elimination found no solution (with rank " + std::to_string(smoothPrimes.size()) + "). If your rank is very low, consider increasing the smoothness bound. Otherwise, produce and retain more smooth numbers.");
    }

    return solutions;
  }

  BigInteger solveCongruence(const std::vector<size_t>& solutionVec)
  {
    // x^2 % toFactor = y^2
    BigInteger x = 1U;
    for (const size_t& idx : solutionVec) {
      x *= smoothNumberKeys[idx];
    }
    const BigInteger y = sqrt((x * x) % toFactor);
    // Uncomment to validate our math, overall
    if ((y * y) != ((x * x) % toFactor)) {
      throw std::runtime_error("Math is not self-consistent!");
    }

    // Check x + y
    BigInteger factor = gcd(toFactor, x + y);
    if ((factor > 1U) && (factor < toFactor)) {
      return factor;
    }

    // Avoid division by 0
    if (x != y) {
      // Check x - y
      return gcd(toFactor, x - y);
    }

    return 1U;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                              WRITTEN WITH HELP FROM ELARA (GPT) ABOVE                                  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  BigInteger solveForFactor() {
    // Gaussian elimination is used to create a perfect square of the residues.
    if (smoothNumberKeys.empty()) {
        throw std::runtime_error("No smooth numbers found. Sieve more, or increase smoothness bound to reduce selectiveness. (The sieving bound multiplier is equivalent to that many times the square root of the number to factor, for calculated numerical range above an offset of the square root of the number to factor.)");
    }

    const std::vector<std::vector<size_t>> solutions = gaussianElimination();
    for (const std::vector<size_t>& solution : solutions) {
      const BigInteger factor = solveCongruence(solution);
      if ((factor > 1U) && (factor < toFactor)) {
        return factor;
      }
    }

    // Depending on row count, a successful result should be nearly guaranteed,
    // but we default to no solution.
    throw std::runtime_error("No solution produced a congruence of squares. (We found " + std::to_string(solutions.size()) + " solutions, and even 1 should often be enough.)");
  }

  // Compute the prime factorization modulo 2
  boost::dynamic_bitset<size_t> factorizationParityVector(BigInteger num) {
    boost::dynamic_bitset<size_t> vec(smoothPrimes.size(), 0U);
    std::vector<size_t> spids(smoothPrimes.size());
    std::iota(spids.begin(), spids.end(), 0);
    while (true) {
      // Proceed in steps of the GCD with the smooth prime wheel radius.
      BigInteger factor = gcd(num, smoothWheelRadius);
      if (factor == 1U) {
        break;
      }
      num /= factor;
      // Remove smooth primes from factor.
      // (The GCD is necessarily smooth.)
      std::vector<size_t> nspids(spids);
      for (size_t rpi = 0U; rpi < spids.size(); ++rpi) {
        const size_t pi = spids.size() - (rpi + 1U);
        const size_t& p = smoothPrimes[spids[pi]];
        if (factor % p) {
          // Once a preamble factor is found not to be present,
          // there's no longer use trying for it on the next iteration.
          nspids.erase(nspids.begin() + pi);
          continue;
        }
        factor /= p;
        vec.flip(pi);
        if (factor == 1U) {
          // The step is fully factored.
          // (This case is always reached.)
          nspids.erase(nspids.begin(), nspids.begin() + pi);
          break;
        }
      }
      if (num == 1U) {
        // The number is fully factored and smooth.
        return vec;
      }
      // Update smoothPrimes IDs
      spids = nspids;
    }
    if (num != 1U) {
      // The number was not fully factored, because it is not smooth.
      // We reject it as a sieving candidate.
      return boost::dynamic_bitset<size_t>();
    }

    // This number is smooth, and we return its factorization parity.
    return vec;
  }
};

std::string find_a_factor(std::string toFactorStr, size_t method, size_t nodeCount, size_t nodeId, size_t gearFactorizationLevel, size_t wheelFactorizationLevel,
                          double sievingBoundMultiplier, double smoothnessBoundMultiplier, size_t gaussianEliminationRowOffset, bool checkSmallFactors) {
  // Validation section
  if (method > 1U) {
    std::cout << "Mode number " << method << " not implemented. Defaulting to FACTOR_FINDER." << std::endl;
    method = 1U;
  }
  const bool isFactorFinder = (method > 0U);
  if (!wheelFactorizationLevel) {
    wheelFactorizationLevel = 1U;
  } else if (!method && (wheelFactorizationLevel > 17U)) {
    wheelFactorizationLevel = 13U;
    std::cout << "Warning: Wheel factorization limit for PRIME_PROVER method is 17. (Parameter will be ignored and default to 17.)" << std::endl;
  }
  if (!gearFactorizationLevel) {
    gearFactorizationLevel = 1U;
  } else if (gearFactorizationLevel < wheelFactorizationLevel) {
    gearFactorizationLevel = wheelFactorizationLevel;
    std::cout << "Warning: Gear factorization level must be at least as high as wheel level. (Parameter will be ignored and default to wheel level.)" << std::endl;
  }
  if (sievingBoundMultiplier > 1.0) {
    sievingBoundMultiplier = 1.0;
    std::cout << "Warning: Sieving bound multiplier was set higher than 1.0. A setting of 1.0 indicates to use the full sieving range. (Parameter will be ignored and default to 1.0.)";
  }

  // Convert number to factor from string.
  const BigInteger toFactor(toFactorStr);

  // The largest possible discrete factor of "toFactor" is its square root (as with any integer).
  const BigInteger sqrtN = sqrt(toFactor);
  if (sqrtN * sqrtN == toFactor) {
    return boost::lexical_cast<std::string>(sqrtN);
  }

  // This level default (scaling) was suggested by Elara (OpenAI GPT).
  const double N = toFactor.convert_to<double>();
  const double logN = log(N);
  const BigInteger primeCeilingBigInt = (BigInteger)(smoothnessBoundMultiplier * pow(exp(0.5 * std::sqrt(logN * log(logN))), std::sqrt(2.0) / 4) + 0.5);
  const size_t primeCeiling = (size_t)primeCeilingBigInt;
  if (((BigInteger)primeCeiling) != primeCeilingBigInt) {
    throw std::runtime_error("Your primes are out of size_t range! (Your formula smoothness bound calculates to be " + boost::lexical_cast<std::string>(primeCeilingBigInt) + ".) Consider lowering your smoothness bound, since it's unlikely you want to sieve for primes above 2 to the 64th power, but, if so, you can modify the SieveOfEratosthenes() code slightly to allow for this.");
  }
  // This uses very little memory and time, to find primes.
  std::vector<size_t> primes = SieveOfEratosthenes(primeCeiling);
  // "it" is the end-of-list iterator for a list up-to-and-including wheelFactorizationLevel.
  const auto itw = std::upper_bound(primes.begin(), primes.end(), wheelFactorizationLevel);
  const auto itg = std::upper_bound(primes.begin(), primes.end(), gearFactorizationLevel);
  const size_t wgDiff = std::distance(itw, itg);

  if (checkSmallFactors && !nodeId) {
    // This is simply trial division up to the ceiling.
    std::mutex trialDivisionMutex;
    BigInteger result = 1U;
    for (size_t primeIndex = 0U; (primeIndex < primes.size()) && (result == 1U); primeIndex += 64U) {
      dispatch.dispatch([&toFactor, &primes, &result, &trialDivisionMutex, primeIndex]() -> bool {
        const size_t maxLcv = std::min(primeIndex + 64U, primes.size());
        for (size_t pi = primeIndex; pi < maxLcv; ++pi) {
          const size_t& currentPrime = primes[pi];
          if (!(toFactor % currentPrime)) {
            std::lock_guard<std::mutex> lock(trialDivisionMutex);
            result = currentPrime;
            return true;
          }
        }
        return false;
      });
    }
    dispatch.finish();
    // If we've checked all primes below the square root of toFactor, then it's prime.
    if ((result != 1U) || (toFactor <= (primeCeiling * primeCeiling))) {
      return boost::lexical_cast<std::string>(result);
    }
  }

  // Set up wheel factorization (or "gear" factorization)
  std::vector<size_t> gearFactorizationPrimes(primes.begin(), itg);
  std::vector<size_t> wheelFactorizationPrimes(primes.begin(), itw);
  // Primes are only present in range above wheel factorization level
  std::vector<size_t> smoothPrimes;
  if (isFactorFinder) {
    smoothPrimes = selectFactorBase(toFactor, primes);
    if (smoothPrimes.empty()) {
      throw std::runtime_error("No smooth primes found under bound. (The formula smoothness bound calculates to " + std::to_string(primeCeiling) + ".) Increase the smoothness bound multiplier, unless this is in range of check_small_factors=True.");
    }
    for (const size_t& sp : smoothPrimes) {
      const auto& gfpit = std::find(gearFactorizationPrimes.begin(), gearFactorizationPrimes.end(), sp);
      if (gfpit != gearFactorizationPrimes.end()) {
         gearFactorizationPrimes.erase(gfpit);
      }
      const auto& wfpit = std::find(wheelFactorizationPrimes.begin(), wheelFactorizationPrimes.end(), sp);
      if (wfpit != wheelFactorizationPrimes.end()) {
        wheelFactorizationPrimes.erase(wfpit);
      }
    }
  }
  // From 1, this is a period for wheel factorization
  BigInteger biggestWheelBigInt = 1U;
  for (const size_t &wp : gearFactorizationPrimes) {
    biggestWheelBigInt *= (size_t)wp;
  }
  // This is defined globally:
  biggestWheel = (size_t)biggestWheelBigInt;
  if (((BigInteger)biggestWheel) != biggestWheelBigInt) {
    throw std::runtime_error("Wheel is too big! Turn down wheel and/or gear factorization level. (Max is less than 2^64, while calculated wheel has radius" + boost::lexical_cast<std::string>(biggestWheelBigInt) + ".)");
  }

  // Wheel entry count per largest "gear" scales our brute-force range.
  // This is defined globally:
  wheel.clear();
  for (size_t i = 1U; i <= biggestWheel; ++i) {
    if (!isMultiple(i, gearFactorizationPrimes)) {
      wheel.push_back(i);
    }
  }
  if (wheel.empty()) {
    wheel.push_back(1U);
  }
  size_t batchItemCount = wheel.size();
  const size_t minBatch = 256U;
  if (minBatch > batchItemCount) {
    batchItemCount = ((minBatch + batchItemCount - 1U) / batchItemCount) * batchItemCount;
  }
  wheelFactorizationPrimes.clear();
  // These are "gears," for wheel factorization (on top of a "wheel" already in place up to the selected level).
  std::vector<boost::dynamic_bitset<size_t>> inc_seqs = generateGears(gearFactorizationPrimes);
  // We're done with the lowest primes.
  const size_t MIN_RTD_LEVEL = gearFactorizationPrimes.size() - wgDiff;
  const Wheel SMALLEST_WHEEL = wheelByPrimeCardinal(MIN_RTD_LEVEL);
  // Skip multiples removed by wheel factorization.
  inc_seqs.erase(inc_seqs.begin(), inc_seqs.end() - wgDiff);
  gearFactorizationPrimes.clear();

  // For PRIME_PROVER method
  const auto ppBackwardFn = backward(SMALLEST_WHEEL);
  const auto ppForwardFn = forward(SMALLEST_WHEEL);
  const BigInteger ppNodeRange = (((ppBackwardFn(sqrtN) + batchItemCount - 1U) / batchItemCount) + nodeCount - 1U) / nodeCount;
  const size_t ppStartingBatch = ((size_t)ppBackwardFn(primeCeiling)) / batchItemCount;

  // For FACTOR_FINDER method (Quadratic Sieve)
  const size_t rowLimit = smoothPrimes.size() + gaussianEliminationRowOffset;
  BigInteger qsBackwardLowBound = smoothBackwardFn(sqrtN + 1U);
  if (smoothForwardFn(qsBackwardLowBound) < (sqrtN + 1U)) {
    ++qsBackwardLowBound;
  }
  const BigInteger qsNodeRange =((((smoothBackwardFn(sqrtN + (BigInteger)((toFactor - sqrtN).convert_to<double>() * sievingBoundMultiplier + 0.5)) - qsBackwardLowBound)
                                      + batchItemCount - 1U) / batchItemCount) + nodeCount - 1U) / nodeCount;

  // This manages the work of all threads.
  Factorizer worker(toFactor, sqrtN, qsBackwardLowBound,
                    isFactorFinder ? qsNodeRange : ppNodeRange,
                    nodeCount, nodeId,
                    batchItemCount,
                    rowLimit,
                    isFactorFinder ? 0U : ppStartingBatch,
                    smoothPrimes,
                    isFactorFinder ? ((wheel.size() > 1U) ? smoothForwardFn : forward(WHEEL1)) : ppForwardFn,
                    isFactorFinder ? ((wheel.size() > 1U) ? smoothBackwardFn : backward(WHEEL1)) : ppBackwardFn);
  // Square of count of smooth primes, for FACTOR_FINDER batch multiplier base unit, was suggested by Lyra (OpenAI GPT)

  std::vector<std::future<BigInteger>> futures;
  futures.reserve(CpuCount);

  const auto workerFn = [&inc_seqs, &worker, &isFactorFinder] {
    // inc_seq needs to be independent per thread.
    std::vector<boost::dynamic_bitset<size_t>> inc_seqs_clone;
    inc_seqs_clone.reserve(inc_seqs.size());
    for (const boost::dynamic_bitset<size_t> &b : inc_seqs) {
      inc_seqs_clone.emplace_back(b);
    }

    // "Brute force" includes extensive wheel multiplication and can be faster.
    return isFactorFinder ? worker.sievePolynomials(&inc_seqs_clone) : worker.bruteForce(&inc_seqs_clone);
  };

  for (unsigned cpu = 0U; cpu < CpuCount; ++cpu) {
    futures.push_back(std::async(std::launch::async, workerFn));
  }

  for (unsigned cpu = 0U; cpu < futures.size(); ++cpu) {
    const BigInteger r = futures[cpu].get();
    if ((r > 1U) && (r < toFactor)) {
      return boost::lexical_cast<std::string>(r);
    }
  }

  // It's only convenient that a large part of the `FACTOR_FINDER` work
  // happens in a second phase, after a first phase with identical signature.
  if (isFactorFinder) {
    return boost::lexical_cast<std::string>(worker.solveForFactor());
  }

  // We would have already returned if we found a factor.
  return std::to_string(1);
}
} // namespace Qimcifa

using namespace Qimcifa;

PYBIND11_MODULE(_find_a_factor, m) {
  m.doc() = "pybind11 plugin to find any factor of input";
  m.def("_find_a_factor", &find_a_factor, "Finds any nontrivial factor of input");
}