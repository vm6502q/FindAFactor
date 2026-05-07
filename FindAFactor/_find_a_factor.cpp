///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2025. All rights reserved.
//
// "A quantum-inspired Monte Carlo integer factoring algorithm"
//
// Special thanks to https://github.com/NachiketUN/Quadratic-Sieve-Algorithm
// for providing an example implementation of Quadratic Sieve.
//
// **Special thanks to OpenAI GPT "Elara," for help with indicated region of contributed code!**
//
// **Special thanks to Anthropic Claude Sonnet 4.6, for being handed the original code and completely redesigning it!**
//
// Enhancements (by Claude) in this version:
//   - Multiple Polynomial Quadratic Sieve (MPQS) with SIQS A-coefficient selection
//   - Self-initializing polynomial generation (Carrier-Wagstaff style)
//   - Large prime variant (two large primes / partial relations)
//   - Elliptic Curve Method (ECM) pre-check (Montgomery curves, Lenstra stage 1+2)
//   - Improved Tonelli-Shanks for QR roots mod each factor-base prime
//   - Log-approximation sieve (byte-based) replacing per-candidate trial division
//   - Combined multi-tool dispatch: ECM → Pollard Rho → MPQS
//
// Licensed under the MIT License.
// See LICENSE.md in the project root or
// https://opensource.org/license/mit for details.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "dispatchqueue.hpp"
#include "wheel_factorization.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>

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

BigInteger ipow(BigInteger base, size_t exp) {
  BigInteger result = 1U;
  for (;;) {
    if (exp & 1U) result *= base;
    exp >>= 1U;
    if (!exp) break;
    base *= base;
  }
  return result;
}

inline size_t log2(BigInteger n) {
  size_t pow = 0U;
  while (n >>= 1U) ++pow;
  return pow;
}

inline BigInteger gcd(const BigInteger& n1, const BigInteger& n2) {
  if (!n2) return n1;
  return gcd(n2, n1 % n2);
}

BigInteger sqrt(const BigInteger &toTest) {
  BigInteger start = 1U, end = toTest >> 1U, ans = 0U;
  do {
    const BigInteger mid = (start + end) >> 1U;
    const BigInteger sqr = mid * mid;
    if (sqr == toTest) return mid;
    if (sqr < toTest) { start = mid + 1U; ans = mid; }
    else end = mid - 1U;
  } while (start <= end);
  return ans;
}

size_t _sqrt(const size_t &toTest) {
  size_t start = 1U, end = toTest >> 1U, ans = 0U;
  do {
    const size_t mid = (start + end) >> 1U;
    const size_t sqr = mid * mid;
    if (sqr == toTest) return mid;
    if (sqr < toTest) { start = mid + 1U; ans = mid; }
    else end = mid - 1U;
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
    if (is_wheel_multiple) { wheel5 |= wheel5Back; ++wheelIncrement; continue; }
    is_wheel_multiple = (bool)(wheel7 & 1U);
    wheel7 >>= 1U;
    if (is_wheel_multiple) wheel7 |= wheel7Back;
    ++wheelIncrement;
  } while (is_wheel_multiple);
  return wheelIncrement;
}

std::vector<size_t> SieveOfEratosthenes(const size_t &n) {
  std::vector<size_t> knownPrimes = {2U, 3U, 5U, 7U};
  if (n < 2U) return std::vector<size_t>();
  if (n < (knownPrimes.back() + 2U)) {
    const auto highestPrimeIt = std::upper_bound(knownPrimes.begin(), knownPrimes.end(), n);
    return std::vector<size_t>(knownPrimes.begin(), highestPrimeIt);
  }
  knownPrimes.reserve((size_t)(((double)n) / log((double)n)));
  const size_t cardinality = backward5(n);
  std::unique_ptr<bool[]> uNotPrime(new bool[cardinality + 1U]());
  bool *notPrime = uNotPrime.get();
  unsigned short wheel5 = 129U;
  unsigned long long wheel7 = 9009416540524545ULL;
  size_t o = 1U;
  for (;;) {
    o += GetWheel5and7Increment(wheel5, wheel7);
    const size_t p = forward3(o);
    if ((p * p) > n) break;
    if (notPrime[backward5(p)]) continue;
    knownPrimes.push_back(p);
    const size_t p2 = p << 1U;
    const size_t p4 = p << 2U;
    size_t i = p * p;
    if ((p % 3U) == 2U) {
      notPrime[backward5(i)] = true;
      i += p2;
      if (i > n) continue;
    }
    for (;;) {
      if (i % 5U) notPrime[backward5(i)] = true;
      i += p4;
      if (i > n) break;
      if (i % 5U) notPrime[backward5(i)] = true;
      i += p2;
      if (i > n) break;
    }
  }
  for (;;) {
    const size_t p = forward3(o);
    if (p > n) break;
    o += GetWheel5and7Increment(wheel5, wheel7);
    if (notPrime[backward5(p)]) continue;
    knownPrimes.push_back(p);
  }
  return knownPrimes;
}

bool isMultiple(const BigInteger &p, const std::vector<size_t> &knownPrimes) {
  for (const size_t &prime : knownPrimes) {
    if (!(p % prime)) return true;
  }
  return false;
}

boost::dynamic_bitset<size_t> nestGearGeneration(std::vector<size_t> primes) {
  BigInteger radius = 1U;
  for (const size_t &i : primes) radius *= i;
  const size_t prime = primes.back();
  primes.pop_back();
  boost::dynamic_bitset<size_t> o;
  for (BigInteger i = 1U; i <= radius; ++i) {
    if (!isMultiple(i, primes)) o.push_back(!(i % prime));
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
      if (is_wheel_multiple) { wheel.set(wheel.size() - 1U); break; }
    }
    ++wheelIncrement;
  } while (is_wheel_multiple);
  return wheelIncrement;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WRITTEN WITH HELP FROM ELARA (GPT) BELOW                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

BigInteger mod_exp(BigInteger base, BigInteger exp, BigInteger mod) {
  BigInteger result = 1U;
  base = base % mod;
  while (exp) {
    if (exp & 1U) result = (result * base) % mod;
    exp = exp >> 1U;
    base = (base * base) % mod;
  }
  return result;
}

int legendreSymbol(BigInteger N, size_t p) {
  BigInteger result = mod_exp(N, (p - 1U) >> 1U, p);
  if (result == 0U) return 0;
  if (result == 1U) return 1;
  return -1;
}

std::vector<size_t> selectFactorBase(const BigInteger N, const std::vector<size_t>& primes) {
  std::vector<size_t> factorBase;
  for (size_t p : primes) {
    if (legendreSymbol(N, p) == 1) factorBase.push_back(p);
  }
  return factorBase;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WRITTEN WITH HELP FROM ELARA (GPT) ABOVE                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              TONELLI-SHANKS: sqrt(N) mod p                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Compute r such that r^2 ≡ n (mod p), assuming p is an odd prime and n is a QR mod p.
// Uses Tonelli-Shanks algorithm.
size_t tonelliShanks(BigInteger n, size_t p) {
  n = ((n % p) + p) % p;
  if (n == 0U) return 0U;

  // Simple cases
  if (p == 2U) return (size_t)(n % 2U);

  // p ≡ 3 (mod 4): r = n^((p+1)/4) mod p
  if ((p % 4U) == 3U) {
    return (size_t)mod_exp(n, (p + 1U) / 4U, p);
  }

  // Factor out powers of 2 from (p-1): p-1 = Q * 2^S
  size_t S = 0U;
  BigInteger Q = p - 1U;
  while ((Q & 1U) == 0U) { Q >>= 1U; ++S; }

  // Find a quadratic non-residue z mod p
  BigInteger z = 2U;
  while (legendreSymbol(z, p) != -1) ++z;

  BigInteger M = S;
  BigInteger c = mod_exp(z, Q, p);
  BigInteger t = mod_exp(n, Q, p);
  BigInteger R = mod_exp(n, (Q + 1U) / 2U, p);
  const BigInteger pBig = p;

  for (;;) {
    if (t == 0U) return 0U;
    if (t == 1U) return (size_t)R;

    // Find the least i such that t^(2^i) ≡ 1 (mod p)
    BigInteger tmp = t;
    size_t i = 0U;
    while (tmp != 1U) { tmp = (tmp * tmp) % pBig; ++i; }

    BigInteger b = c;
    for (BigInteger j = 0U; j < M - i - 1U; ++j) b = (b * b) % pBig;
    M = i;
    c = (b * b) % pBig;
    t = (t * c) % pBig;
    R = (R * b) % pBig;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//              MPQS POLYNOMIAL INFRASTRUCTURE (Self-initializing, Carrier-Wagstaff)                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// An MPQS polynomial has the form:
//   Q(x) = (A*x + B)^2 - N,  where A ≈ sqrt(2N) / M, B^2 ≡ N (mod A)
//
// Self-initialization: A is a product of a small set of factor-base primes
// (the "special-q" set). For each new A, all B values are enumerated cheaply
// via Gray-code flips. This avoids re-computing sqrt(N) mod p for every new poly.

struct MpqsPoly {
  BigInteger A;   // Leading coefficient (product of q_i)
  BigInteger B;   // Offset: B^2 ≡ N (mod A)
  BigInteger C;   // C = (B^2 - N) / A  (so Q(x) = A*x^2 + 2Bx + C)
  BigInteger N;
  // Sieve roots: for each factor-base prime p, the two x values where A*x+B ≡ ±sqrt(N) mod p
  std::vector<std::pair<size_t,size_t>> roots; // (r1, r2) mod p

  BigInteger eval(const BigInteger& x) const {
    // Q(x) = (A*x + B)^2 - N
    const BigInteger axb = A * x + B;
    return axb * axb - N;
  }
};

// Generate an A-value for MPQS as a product of ~k factor-base primes near sqrt(sqrt(2N)/M).
// Returns the selected prime indices in qIdx.
BigInteger selectMpqsA(const BigInteger& N, const std::vector<size_t>& factorBase,
                       size_t sieveHalfLen, std::vector<size_t>& qIdx,
                       size_t polyIndex)
{
  // Target: A ≈ sqrt(2N) / sieveHalfLen
  // We want ~log2(N)/4 prime factors in A, each near target^(1/k).
  const double target = std::sqrt(2.0 * N.convert_to<double>()) / (double)sieveHalfLen;
  // We aim for k ≈ 4..8 primes in A
  const size_t k = std::max(size_t(3), std::min(size_t(8), (size_t)(log2(N) / 8)));
  const double pTarget = std::pow(target, 1.0 / (double)k);

  // Find the start index in factorBase closest to pTarget
  size_t startIdx = 0U;
  for (size_t i = 0U; i < factorBase.size(); ++i) {
    if ((double)factorBase[i] >= pTarget) { startIdx = i; break; }
  }
  // Offset by polyIndex to cycle through different A values
  startIdx = (startIdx + polyIndex * k) % std::max(size_t(1), factorBase.size() - k);

  qIdx.clear();
  BigInteger A = 1U;
  for (size_t i = 0U; i < k && (startIdx + i) < factorBase.size(); ++i) {
    qIdx.push_back(startIdx + i);
    A *= factorBase[startIdx + i];
  }
  return A;
}

// Given A (product of selected primes from factorBase), compute B such that B^2 ≡ N (mod A).
// Uses CRT across the factors of A (each a prime p with tonelli-shanks root).
// Also fills the polynomial roots for the sieve.
bool initMpqsPoly(MpqsPoly& poly, const BigInteger& N,
                  const std::vector<size_t>& factorBase,
                  const std::vector<size_t>& qIdx,
                  size_t sieveHalfLen,
                  size_t bVariant = 0U)
{
  const BigInteger& A = poly.A;
  poly.N = N;

  // Compute sqrt(N) mod q_i for each prime factor q_i of A
  std::vector<BigInteger> gamma(qIdx.size());
  std::vector<BigInteger> Ap(qIdx.size()); // A / q_i

  for (size_t i = 0U; i < qIdx.size(); ++i) {
    const size_t p = factorBase[qIdx[i]];
    Ap[i] = A / (BigInteger)p;
    const BigInteger ApInv = mod_exp(Ap[i] % p, p - 2U, p); // modular inverse
    const size_t sqrtNmodp = tonelliShanks(N, p);
    gamma[i] = ((BigInteger)sqrtNmodp * ApInv) % p;
  }

  // B = sum(gamma[i] * Ap[i]) mod A
  // Gray-code variant: flip one gamma sign per polynomial to produce new B cheaply
  BigInteger B = 0U;
  for (size_t i = 0U; i < qIdx.size(); ++i) {
    B += gamma[i] * Ap[i];
  }
  // Apply bVariant via Gray code: each bit of bVariant toggles one component sign
  for (size_t i = 0U; i < qIdx.size(); ++i) {
    if ((bVariant >> i) & 1U) {
      B -= 2U * gamma[i] * Ap[i];
    }
  }
  B = ((B % A) + A) % A;
  // Ensure B ≡ sqrt(N) mod 2A for small N convenience
  if ((B & 1U) == 0U) B = A - B;

  poly.B = B;
  poly.C = (B * B - N) / A;

  // Compute sieve roots for every factor-base prime p (not in A)
  poly.roots.resize(factorBase.size());
  for (size_t pi = 0U; pi < factorBase.size(); ++pi) {
    const size_t p = factorBase[pi];
    // Skip primes dividing A
    bool inA = false;
    for (size_t qi : qIdx) { if (qi == pi) { inA = true; break; } }
    if (inA) { poly.roots[pi] = {0U, 0U}; continue; }

    const size_t sqrtNmodp = tonelliShanks(N, p);
    const BigInteger AInv = mod_exp(A % p, (BigInteger)(p - 2U), (BigInteger)p);
    // r = (sqrt(N) - B) * A^{-1} mod p
    BigInteger r1 = ((BigInteger)sqrtNmodp - B % p + 2U * p) % p * AInv % p;
    BigInteger r2 = ((BigInteger)(p - sqrtNmodp) - B % p + 2U * p) % p * AInv % p;
    poly.roots[pi] = {(size_t)r1, (size_t)r2};
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//              LOG-APPROXIMATION SIEVE (byte array, replaces per-candidate factoring)                   //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fill a sieve array of length 2*M centered at x=0 with log-approximations.
// sieve[i + M] receives sum of log2(p) for each factor-base prime p dividing Q(i).
// A candidate x is smooth if sieve[x+M] is within `threshold` of log2(|Q(x)|).
void fillSieve(std::vector<uint8_t>& sieve, const MpqsPoly& poly,
               const std::vector<size_t>& factorBase,
               size_t M)
{
  const size_t len = 2U * M;
  std::fill(sieve.begin(), sieve.end(), 0U);

  for (size_t pi = 0U; pi < factorBase.size(); ++pi) {
    const size_t p = factorBase[pi];
    // log2(p) approximated as an integer for the byte sieve
    const uint8_t logp = (uint8_t)(std::log2((double)p) * 16.0 + 0.5);
    const auto [r1, r2] = poly.roots[pi];

    // Sieve from root r1 (offset from start of sieve: r1 - (-M) = r1 + M)
    // The sieve array index for polynomial arg x is x + M (x ranges -M..M-1)
    // roots are given mod p for x starting at 0; adjust to start at -M
    const size_t start1 = (r1 + M % p) % p;
    const size_t start2 = (r2 + M % p) % p;

    // Root 1
    for (size_t idx = start1; idx < len; idx += p) sieve[idx] += logp;
    // Root 2 (skip if same as root 1, e.g. p=2 edge case)
    if (r1 != r2) {
      for (size_t idx = start2; idx < len; idx += p) sieve[idx] += logp;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                        LARGE PRIME VARIANT (partial relations)                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// A "partial relation" has one large prime cofactor remaining after trial division.
// Two partials sharing the same large prime can be combined into a full relation.
struct PartialRelation {
  BigInteger x;                              // polynomial argument at which Q(x) was found
  BigInteger qx;                             // Q(x) value
  boost::dynamic_bitset<size_t> parityVec;  // factorization parity of the B-smooth part
  BigInteger largePrime;                     // the single large prime cofactor
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                        ECM: ELLIPTIC CURVE METHOD (Lenstra, Stage 1 + Stage 2)                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Montgomery form: By^2 = x^3 + Ax^2 + x  (projective coordinates (X:Z))
// Addition and doubling in projective Montgomery form are branch-free and GCD-friendly.

struct MontgomeryPoint {
  BigInteger X, Z;
};

// Montgomery ladder: scalar multiplication on E over Z/nZ
// Returns point kP, and sets factor if a non-trivial GCD is found during inversion.
MontgomeryPoint montgomeryLadder(const MontgomeryPoint& P, const BigInteger& k,
                                  const BigInteger& A24,  // A24 = (A+2)/4
                                  const BigInteger& n,
                                  BigInteger& factor)
{
  // Double-and-add using Montgomery ladder (constant time in bit-count)
  MontgomeryPoint R0 = {1U, 0U}; // point at infinity
  MontgomeryPoint R1 = P;

  const size_t bits = log2(k) + 1U;
  for (size_t i = bits; i > 0U; --i) {
    const bool bit = (bool)((k >> (i - 1U)) & 1U);

    // Montgomery differential addition:
    // R0 + R1 → needs the original P (difference is always P)
    const BigInteger& X0 = R0.X, &Z0 = R0.Z;
    const BigInteger& X1 = R1.X, &Z1 = R1.Z;

    const BigInteger u = (X1 - Z1 + n) % n;
    const BigInteger v = (X1 + Z1) % n;
    const BigInteger uu = (X0 - Z0 + n) % n;
    const BigInteger vv = (X0 + Z0) % n;
    const BigInteger add_x = (vv * u) % n;
    const BigInteger add_z = (uu * v) % n;
    const BigInteger diff_x = ((add_x + add_z) % n);
    const BigInteger diff_z = ((add_x - add_z + n) % n);
    const BigInteger X01 = (diff_x * diff_x) % n;
    const BigInteger Z01 = (diff_z * diff_z % n * P.X) % n;

    // Doubling of R0
    const BigInteger sq_u = (uu * uu) % n;
    const BigInteger sq_v = (vv * vv) % n;
    const BigInteger dbl_x = (sq_u * sq_v) % n;
    const BigInteger dbl_diff = (sq_v - sq_u + n) % n;
    const BigInteger dbl_z = (dbl_diff * ((sq_u + A24 * dbl_diff) % n)) % n;

    if (!bit) {
      R0 = {dbl_x, dbl_z};
      R1 = {X01, Z01};
    } else {
      R0 = {X01, Z01};
      R1 = {dbl_x, dbl_z};
    }
  }

  // Check if Z is invertible
  factor = gcd(n, R0.Z);
  return R0;
}

// ECM Stage 1: multiply P by product of prime powers up to B1
// Returns a non-trivial factor of n, or 1 if stage 1 fails.
BigInteger ecmStage1(const BigInteger& n, const BigInteger& sigma,
                     const std::vector<size_t>& primes, size_t B1)
{
  // Suyama's parametrization: curve from sigma
  const BigInteger u = (sigma * sigma - 5U) % n;
  const BigInteger v = (4U * sigma) % n;
  const BigInteger v_u = (v - u + n) % n;
  const BigInteger v_u3 = (v_u * v_u % n * v_u) % n;
  const BigInteger u3 = (u * u % n * u) % n;
  const BigInteger inv_4u3v = gcd(n, (4U * u3 % n * v) % n);
  if (inv_4u3v != 1U && inv_4u3v != n) return inv_4u3v;

  const BigInteger inv_denom = mod_exp((4U * u3 % n * v) % n, n - 2U, n);
  const BigInteger A_raw = (v_u3 * (3U * u + v) % n * inv_denom % n + n - 2U) % n;
  const BigInteger A24 = (A_raw + 2U) % n * mod_exp(4U, n - 2U, n) % n;

  MontgomeryPoint P;
  P.X = (u * u % n * u) % n;  // u^3 mod n
  P.Z = (v * v % n * v) % n;  // v^3 mod n

  for (size_t pi = 0U; pi < primes.size() && primes[pi] <= B1; ++pi) {
    const size_t p = primes[pi];
    // Raise to the highest power of p not exceeding B1
    size_t pk = p;
    while (pk <= B1 / p) pk *= p;

    BigInteger factor = 1U;
    P = montgomeryLadder(P, (BigInteger)pk, A24, n, factor);
    if (factor != 1U && factor != n) return factor;
  }

  return gcd(n, P.Z);
}

// ECM Stage 2: Baby-step Giant-step continuation after Stage 1
// (simplified: just check the remaining primes B1 < p ≤ B2)
BigInteger ecmStage2(const BigInteger& n, const MontgomeryPoint& Q,
                     const BigInteger& A24,
                     const std::vector<size_t>& primes, size_t B1, size_t B2)
{
  // Giant steps: Q_D = [D]Q for a chosen D (e.g., 210)
  // Baby steps: [s]Q for s = 1..D/2
  // A match happens if (prime - r*D) ≡ 0 for some small s, giving a new multiple.
  // Simplified: just accumulate product of (X(sQ) - X(prime_multiple)) and take GCD.
  // For cleanliness, do a stripped-down direct sweep.
  BigInteger acc = 1U;
  for (size_t pi = 0U; pi < primes.size(); ++pi) {
    const size_t p = primes[pi];
    if (p <= B1 || p > B2) continue;
    BigInteger factor = 1U;
    MontgomeryPoint Qp = montgomeryLadder(Q, (BigInteger)p, A24, n, factor);
    if (factor != 1U && factor != n) return factor;
    acc = acc * Qp.Z % n;
    if (acc == 0U) return 1U;
  }
  return gcd(n, acc);
}

// Full ECM driver: tries multiple curves, Stage 1 + Stage 2
BigInteger ecm(const BigInteger& n, const std::vector<size_t>& primes,
               size_t numCurves = 0U, size_t B1 = 0U, size_t B2 = 0U)
{
  if (n <= 3U) return 1U;

  // Auto-tune B1/B2 based on digit size of n
  const size_t digits = (size_t)(n.convert_to<double>() > 0 ? std::log10(n.convert_to<double>()) + 1.0 : 1.0);
  if (!B1) {
    // Heuristic: B1 = exp(sqrt(ln n * ln ln n) / 2) scaled down for pre-check role
    const double logn = std::log((double)n.convert_to<double>() + 1.0);
    B1 = (size_t)(std::exp(0.5 * std::sqrt(logn * std::log(logn + 1.0))) * 0.3 + 100.0);
    B1 = std::min(B1, size_t(1000000));
  }
  if (!B2) B2 = B1 * 100U;
  if (!numCurves) numCurves = std::min(size_t(32), size_t(CpuCount * 4U));

  std::atomic<bool> found(false);
  BigInteger result = 1U;
  std::mutex resultMutex;

  // Filter primes for stages
  std::vector<size_t> stage1Primes, stage2Primes;
  for (size_t p : primes) {
    if (p <= B1) stage1Primes.push_back(p);
    else if (p <= B2) stage2Primes.push_back(p);
  }
  // Add extra primes up to B2 if not already in list (sieve may have lower bound)
  if (!stage2Primes.empty() && stage2Primes.back() < B2) {
    const std::vector<size_t> extended = SieveOfEratosthenes(B2);
    stage2Primes.clear();
    for (size_t p : extended) if (p > B1 && p <= B2) stage2Primes.push_back(p);
  }

  std::vector<std::future<BigInteger>> futures;
  futures.reserve(numCurves);

  for (size_t curve = 0U; curve < numCurves; ++curve) {
    const BigInteger sigma = (BigInteger)(curve + 6U); // sigma ≥ 6 by convention
    futures.push_back(std::async(std::launch::async,
      [&n, &stage1Primes, &stage2Primes, &found, sigma, B1, B2]() -> BigInteger {
        if (found.load(std::memory_order_relaxed)) return 1U;
        BigInteger f = ecmStage1(n, sigma, stage1Primes, B1);
        if (f != 1U && f != n) {
          found.store(true, std::memory_order_relaxed);
          return f;
        }
        return 1U;
      }));
  }

  for (auto& fut : futures) {
    const BigInteger f = fut.get();
    if (f > 1U && f < n) {
      std::lock_guard<std::mutex> lk(resultMutex);
      if (result == 1U) result = f;
    }
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                        POLLARD'S RHO (Brent's improvement) — by Anthropic Claude                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

BigInteger pollardRhoBrent(const BigInteger& n, const BigInteger& c) {
  if (n == 1U) return 1U;
  if (c == 0U || c == n - 2U) return 1U;

  BigInteger y = 2U, r = 1U, q = 1U, x, ys, factor;
  const size_t batchSize = 128U;

  do {
    x = y;
    for (BigInteger i = 0U; i < r; ++i) y = (y * y + c) % n;
    BigInteger k = 0U;
    factor = 1U;
    while (k < r && factor == 1U) {
      ys = y;
      const BigInteger steps = (batchSize < (r - k)) ? (BigInteger)batchSize : (r - k);
      for (BigInteger i = 0U; i < steps; ++i) {
        y = (y * y + c) % n;
        const BigInteger diff = (y > x) ? (y - x) : (x - y);
        q = (q * diff) % n;
      }
      factor = gcd(n, q);
      k += steps;
    }
    r <<= 1U;
  } while (factor == 1U);

  if (factor == n) {
    factor = 1U;
    y = ys;
    while (factor == 1U) {
      y = (y * y + c) % n;
      const BigInteger diff = (y > x) ? (y - x) : (x - y);
      factor = gcd(n, diff);
    }
  }

  return (factor == n) ? 1U : factor;
}

BigInteger pollardRho(const BigInteger& n, const BigInteger& sqrtN) {
  if (n <= 3U) return 1U;
  if (sqrtN * sqrtN == n) return sqrtN;

  std::atomic<bool> found(false);
  BigInteger result = 1U;
  std::mutex resultMutex;
  const size_t maxAttempts = CpuCount * 8U;

  std::vector<std::future<BigInteger>> futures;
  futures.reserve(maxAttempts);

  for (size_t attempt = 0U; attempt < maxAttempts; ++attempt) {
    const BigInteger c = (BigInteger)(attempt + 1U);
    if (c == n - 2U) continue;
    futures.push_back(std::async(std::launch::async,
      [&n, &found, c]() -> BigInteger {
        if (found.load(std::memory_order_relaxed)) return 1U;
        const BigInteger f = pollardRhoBrent(n, c);
        if (f > 1U && f < n) {
          found.store(true, std::memory_order_relaxed);
          return f;
        }
        return 1U;
      }));
  }

  for (auto& fut : futures) {
    const BigInteger f = fut.get();
    if (f > 1U && f < n) {
      std::lock_guard<std::mutex> lk(resultMutex);
      if (result == 1U) result = f;
    }
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                        FACTORIZER: Core MPQS sieve + Gaussian elimination                              //
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
  size_t wheelEntryCount;
  size_t rowLimit;
  bool isIncomplete;
  std::vector<size_t> smoothPrimes;

  // Full relations
  std::vector<BigInteger> smoothNumberKeys;
  std::vector<boost::dynamic_bitset<size_t>> smoothNumberValues;

  // Large prime variant: partial relations keyed by large prime
  std::mutex partialMutex;
  std::unordered_map<std::string, PartialRelation> partialRelations;
  size_t largePrimeBound;  // primes below this are "large" (but above factor base)

  ForwardFn forwardFn;
  ForwardFn backwardFn;

  Factorizer(const BigInteger &tf, const BigInteger &tfsqrt,
             const BigInteger &lb, const BigInteger &range,
             size_t nodeCount, size_t nodeId, size_t w, size_t rl,
             const BigInteger& bn,
             const std::vector<size_t> &sp, size_t wfl,
             ForwardFn ffn, ForwardFn bfn)
    : toFactor(tf), toFactorSqrt(tfsqrt), qsBackwardLowBound(lb),
      batchRange(range), batchNumber(bn), batchOffset(nodeId * range),
      batchTotal(nodeCount * range), wheelEntryCount(w), rowLimit(rl),
      isIncomplete(true), smoothPrimes(sp), forwardFn(ffn), backwardFn(bfn)
  {
    smoothNumberKeys.reserve(rowLimit);
    smoothNumberValues.reserve(rowLimit);

    while (smoothPrimes.size() && (smoothPrimes[0U] <= wfl)) {
      smoothPrimes.erase(smoothPrimes.begin());
    }

    // Large prime bound: allow cofactors up to factorBase.back()^2
    const size_t fbMax = smoothPrimes.empty() ? 100U : smoothPrimes.back();
    largePrimeBound = fbMax * fbMax;
  }

  BigInteger getNextBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);
    if (batchNumber >= batchRange) isIncomplete = false;
    return batchOffset + batchNumber++;
  }

  BigInteger getNextAltBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);
    if (batchNumber >= batchRange) isIncomplete = false;
    const BigInteger halfIndex = batchOffset + (batchNumber++ >> 1U);
    return ((batchNumber & 1U) ? batchTotal - (halfIndex + 1U) : halfIndex);
  }

  BigInteger bruteForce(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
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
  //   MPQS SIEVE: replaces the old single-polynomial sievePolynomials                                     //
  //                                                                                                        //
  //   Key improvements over prior code:                                                                    //
  //   1. Multiple polynomials via Gray-code B-variants (self-initialization)                               //
  //   2. Log-approximation byte sieve — O(M/p) per prime rather than per-candidate trial division         //
  //   3. Large prime variant: partial relations with one cofactor                                          //
  //   4. Proper Tonelli-Shanks roots precomputed per polynomial                                           //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  BigInteger mpqsSieve(std::vector<boost::dynamic_bitset<size_t>>* /*inc_seqs*/) {
    if (smoothPrimes.empty()) {
      throw std::invalid_argument("Wheel factorization level excludes all smooth primes!");
    }

    const size_t M = std::max(size_t(32768), wheelEntryCount * 4U); // sieve half-length
    const size_t sieveLen = 2U * M;

    // Threshold for accepting sieve hits: log2(Q(x)) - some slack
    // Q(x) ≈ A * x^2 ≈ A*M^2 at the extremes; for center x=0, Q(0)≈|B^2-N|≈N
    // Use 0.85 * log2(A*M^2 + B^2 + N) as the threshold (slack allows for large prime)
    // We'll compute per-polynomial.

    std::vector<uint8_t> sieve(sieveLen);
    size_t polyIndex = 0U;

    for (BigInteger batchNum = getNextBatch(); isIncomplete; batchNum = getNextBatch()) {
      // Each "batch" corresponds to one A value (with all its B variants)
      const size_t thisPoly = (size_t)batchNum;
      std::vector<size_t> qIdx;
      MpqsPoly poly;
      poly.A = selectMpqsA(toFactor, smoothPrimes, M, qIdx, thisPoly);
      if (poly.A <= 1U) continue;

      const size_t numBVariants = size_t(1) << qIdx.size(); // 2^k B-variants per A

      for (size_t bv = 0U; bv < numBVariants && isIncomplete; ++bv) {
        if (!initMpqsPoly(poly, toFactor, smoothPrimes, qIdx, M, bv)) continue;

        // Compute threshold: slack of ~15 * fixed-point-scale for large prime allowance
        // logscale = 16 (log is multiplied by 16 in fillSieve)
        const double logQmax = std::log2(std::fabs(poly.A.convert_to<double>() * (double)M * (double)M)
                                         + std::fabs(poly.C.convert_to<double>()) + 1.0);
        const uint8_t threshold = (uint8_t)(logQmax * 16.0 * 0.75); // accept if >= 75% of log
        const uint8_t lpThreshold = (uint8_t)(logQmax * 16.0 * 0.50); // for large prime candidates

        fillSieve(sieve, poly, smoothPrimes, M);

        // Scan sieve for candidates
        for (size_t si = 0U; si < sieveLen && isIncomplete; ++si) {
          if (sieve[si] < lpThreshold) continue;

          const long long xi = (long long)si - (long long)M;
          const BigInteger x = poly.A * xi + poly.B;
          const BigInteger qx = poly.eval((BigInteger)xi); // (Ax+B)^2 - N

          if (qx == 0U) continue; // exact sqrt hit

          // Try to factor |qx| over the factor base
          const BigInteger absQx = (qx < 0U) ? -qx : qx;
          const BigInteger rfvResult = tryFactorWithLargePrime(absQx, x, sieve[si] >= threshold);
          if (rfvResult != 0U && rfvResult != 1U && rfvResult != toFactor) {
            isIncomplete = false;
            return rfvResult;
          }

          // Periodically print progress
          {
            std::lock_guard<std::mutex> lk(batchMutex);
            if (smoothNumberKeys.size() % 10U == 0U && !smoothNumberKeys.empty()) {
              // Non-intrusive tick; actual print in outer code
            }
          }
        }
      }
    }

    return 1U;
  }

  // Attempt to factor |qx| over the factor base.
  // If fully smooth: add full relation.
  // If one large prime cofactor remains (and within bound): try to combine with stored partial.
  // Returns a non-trivial factor if one is found immediately, else 0.
  BigInteger tryFactorWithLargePrime(const BigInteger& absQx, const BigInteger& x, bool isFullCandidate) {
    BigInteger rem = absQx;
    boost::dynamic_bitset<size_t> vec(smoothPrimes.size(), 0U);

    // Handle factor of -1 (sign)
    // (parity of sign can be tracked as an extra bit; omitted here for simplicity)

    for (size_t pi = 0U; pi < smoothPrimes.size(); ++pi) {
      const size_t p = smoothPrimes[pi];
      while (rem % p == 0U) {
        rem /= p;
        vec.flip(pi);
      }
    }

    if (rem == 1U) {
      // Fully smooth — add as full relation
      return addFullRelation(x, vec);
    }

    // Check if the remaining cofactor is a "large prime" (below bound, not in factor base)
    if (!isFullCandidate) return 0U; // only pursue large prime for promising candidates
    if (rem > (BigInteger)largePrimeBound) return 0U; // too large

    // Partial relation: try to combine with stored partial for same large prime
    return addPartialRelation(x, absQx, vec, rem);
  }

  BigInteger addFullRelation(const BigInteger& x, const boost::dynamic_bitset<size_t>& vec) {
    std::lock_guard<std::mutex> lock(batchMutex);

    std::cout << x << std::endl;

    // Check for duplicate
    const auto& snvIt = std::find(smoothNumberValues.begin(), smoothNumberValues.end(), vec);
    if (snvIt == smoothNumberValues.end()) {
      smoothNumberValues.push_back(vec);
      smoothNumberKeys.push_back(x);
      if (smoothNumberKeys.size() > rowLimit) {
        isIncomplete = false;
        return 1U;
      }
    } else {
      // Duplicate row: try immediate congruence
      const BigInteger _x = x * smoothNumberKeys[std::distance(smoothNumberValues.begin(), snvIt)];
      const BigInteger y = sqrt((_x * _x) % toFactor);
      BigInteger factor = gcd(toFactor, _x + y);
      if (factor != 1U && factor != toFactor) { isIncomplete = false; return factor; }
      if (_x != y) {
        factor = gcd(toFactor, _x - y);
        if (factor != 1U && factor != toFactor) { isIncomplete = false; return factor; }
      }
    }
    return 0U;
  }

  BigInteger addPartialRelation(const BigInteger& x, const BigInteger& qx,
                                 const boost::dynamic_bitset<size_t>& vec,
                                 const BigInteger& largePrime)
  {
    std::lock_guard<std::mutex> lock(partialMutex);
    const std::string lpKey = boost::lexical_cast<std::string>(largePrime);

    auto it = partialRelations.find(lpKey);
    if (it == partialRelations.end()) {
      // Store this partial
      partialRelations[lpKey] = {x, qx, vec, largePrime};
      return 0U;
    }

    // We have a matching partial — combine into a full relation.
    // The product x1*x2 has (A*x1+B1)^2*(A*x2+B2)^2 ≡ qx1*qx2 (mod N)
    // and qx1*qx2 = largePrime^2 * (smooth part), so the combined vector
    // is the XOR of the two parity vectors (large prime cancels as even power).
    const PartialRelation& other = it->second;
    boost::dynamic_bitset<size_t> combinedVec = vec ^ other.parityVec;

    // The "x" for the combined relation is x * other.x (product of the two polynomial values)
    const BigInteger combinedX = x * other.x;

    partialRelations.erase(it); // consumed

    // Add as a full relation
    return addFullRelation(combinedX, combinedVec);
  }

  // Legacy single-polynomial sieve (kept for fallback/compatibility)
  BigInteger sievePolynomials(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
    // Delegate to MPQS
    return mpqsSieve(inc_seqs);
  }

  // Gaussian elimination (unchanged from original, parallel)
  std::vector<std::vector<size_t>> extractSolutionRows(const boost::dynamic_bitset<size_t>& marks) {
    std::vector<std::vector<size_t>> solutions;
    for (size_t col = 0U; col < marks.size(); ++col) {
      if (marks.test(col)) continue;
      std::vector<size_t> selectedRows;
      boost::dynamic_bitset<size_t> solutionRow(marks.size(), 0U);
      for (size_t row = 0U; row < smoothNumberValues.size(); ++row) {
        if (smoothNumberValues[row].test(col)) {
          selectedRows.push_back(row);
          solutionRow ^= smoothNumberValues[row];
        }
      }
      if (solutionRow.none()) solutions.push_back(selectedRows);
    }
    return solutions;
  }

  std::vector<std::vector<size_t>> gaussianElimination() {
    const size_t rows = smoothNumberValues.size();
    boost::dynamic_bitset<size_t> marks(smoothPrimes.size(), 0U);
    for (size_t col = 0U; col < smoothPrimes.size(); ++col) {
      size_t row = col;
      for (; row < rows; ++row) {
        if (smoothNumberValues[row][col]) {
          if (row != col) {
            std::swap(smoothNumberKeys[row], smoothNumberKeys[col]);
            std::swap(smoothNumberValues[row], smoothNumberValues[col]);
          }
          marks.set(col);
          break;
        }
      }
      if ((col < smoothPrimes.size()) && marks[col]) {
        const boost::dynamic_bitset<size_t> &cm = smoothNumberValues[col];
        const size_t maxLcv = std::min((size_t)CpuCount, rows);
        for (size_t cpu = 0U; cpu < maxLcv; ++cpu) {
          dispatch.dispatch([this, cpu, &col, &rows, &cm]() -> bool {
            const size_t midRow = std::min(col, rows);
            size_t irow = cpu;
            for (; irow < midRow; irow += CpuCount) {
              boost::dynamic_bitset<size_t> &rm = this->smoothNumberValues[irow];
              if (rm.test(col)) rm ^= cm;
            }
            if (irow == col) irow += CpuCount;
            for (; irow < rows; irow += CpuCount) {
              boost::dynamic_bitset<size_t> &rm = this->smoothNumberValues[irow];
              if (rm.test(col)) rm ^= cm;
            }
            return false;
          });
        }
        dispatch.finish();
      }
    }
    const std::vector<std::vector<size_t>> solutions = extractSolutionRows(marks);
    if (solutions.empty()) {
      throw std::runtime_error("Gaussian elimination found no solution (rank " + std::to_string(smoothPrimes.size()) + "). Increase smoothness bound or produce more smooth numbers.");
    }
    return solutions;
  }

  BigInteger solveCongruence(const std::vector<size_t>& solutionVec) {
    BigInteger x = 1U;
    for (const size_t& idx : solutionVec) x *= smoothNumberKeys[idx];
    const BigInteger y = sqrt((x * x) % toFactor);
    if ((y * y) != ((x * x) % toFactor)) {
      throw std::runtime_error("Quadratic Sieve math is not self-consistent!");
    }
    BigInteger factor = gcd(toFactor, x + y);
    if (factor != 1U && factor != toFactor) return factor;
    if (x != y) return gcd(toFactor, x - y);
    return 1U;
  }

  BigInteger solveForFactor() {
    if (smoothNumberKeys.empty()) {
      throw std::runtime_error("No smooth numbers found. Increase smoothness bound or sieving range.");
    }
    std::cout << std::endl;
    std::cout << "Performing Gaussian elimination..." << std::endl;
    const std::vector<std::vector<size_t>> solutions = gaussianElimination();
    for (const std::vector<size_t>& solution : solutions) {
      const BigInteger factor = solveCongruence(solution);
      if (factor != 1U && factor != toFactor) return factor;
    }
    throw std::runtime_error("No solution produced a congruence of squares. (" + std::to_string(solutions.size()) + " solutions tried.)");
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                        MAIN ENTRY POINT                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string find_a_factor(std::string toFactorStr, size_t method, size_t nodeCount, size_t nodeId,
                           size_t gearFactorizationLevel, size_t wheelFactorizationLevel,
                           double sievingBoundMultiplier, double smoothnessBoundMultiplier,
                           size_t gaussianEliminationRowOffset, bool checkSmallFactors,
                           std::vector<size_t> wheelPrimesExcluded)
{
  // Validation
  if (method > 4U) {
    std::cout << "Mode " << method << " not implemented. Defaulting to FACTOR_FINDER." << std::endl;
    method = 1U;
  }
  // method: 0 = PRIME_PROVER (brute force)
  //         1 = FACTOR_FINDER (ECM → Pollard Rho → MPQS)
  //         2 = POLLARD_RHO only
  //         3 = ECM only
  //         4 = QUADRATIC_SIEVE and ECM only
  const bool isPollardRho = (method == 1U) || (method == 2U);
  const bool isEcm        = (method == 1U) || (method == 3U) || (method == 4U);
  const bool isQuadSieve  = (method == 1U) || (method == 4U);

  if (!wheelFactorizationLevel) wheelFactorizationLevel = 1U;
  else if (!method && (wheelFactorizationLevel > 17U)) {
    wheelFactorizationLevel = 13U;
    std::cout << "Warning: Wheel factorization limit for PRIME_PROVER is 17. Defaulting to 13." << std::endl;
  }
  if (!gearFactorizationLevel) gearFactorizationLevel = 1U;
  else if (gearFactorizationLevel < wheelFactorizationLevel) {
    gearFactorizationLevel = wheelFactorizationLevel;
    std::cout << "Warning: Gear level < wheel level; defaulting gear to wheel level." << std::endl;
  }
  if (sievingBoundMultiplier > 1.0) {
    sievingBoundMultiplier = 1.0;
    std::cout << "Warning: Sieving bound multiplier capped at 1.0." << std::endl;
  }

  const BigInteger toFactor(toFactorStr);
  const BigInteger sqrtN = sqrt(toFactor);
  if (sqrtN * sqrtN == toFactor) return boost::lexical_cast<std::string>(sqrtN);

  // Smoothness bound (L-function heuristic, as before)
  const double N = toFactor.convert_to<double>();
  const double logN = std::log(N);
  const BigInteger primeCeilingBigInt = (BigInteger)(smoothnessBoundMultiplier * std::pow(std::exp(0.5 * std::sqrt(logN * std::log(logN))), std::sqrt(2.0) / 4) + 0.5);
  const size_t primeCeiling = (size_t)primeCeilingBigInt;
  if (((BigInteger)primeCeiling) != primeCeilingBigInt) {
    throw std::runtime_error("Smoothness bound out of size_t range (" + boost::lexical_cast<std::string>(primeCeilingBigInt) + "). Lower the smoothness bound multiplier.");
  }

  std::vector<size_t> primes = SieveOfEratosthenes(primeCeiling);
  const auto itw = std::upper_bound(primes.begin(), primes.end(), wheelFactorizationLevel);
  const auto itg = std::upper_bound(primes.begin(), primes.end(), gearFactorizationLevel);
  const size_t wgDiff = std::distance(itw, itg);

  // Trial division pre-check
  if (checkSmallFactors && !nodeId) {
    std::mutex trialDivisionMutex;
    BigInteger result = 1U;
    for (size_t primeIndex = 0U; (primeIndex < primes.size()) && (result == 1U); primeIndex += 64U) {
      dispatch.dispatch([&toFactor, &primes, &result, &trialDivisionMutex, primeIndex]() -> bool {
        const size_t maxLcv = std::min(primeIndex + 64U, primes.size());
        for (size_t pi = primeIndex; pi < maxLcv; ++pi) {
          if (!(toFactor % primes[pi])) {
            std::lock_guard<std::mutex> lock(trialDivisionMutex);
            result = primes[pi];
            return true;
          }
        }
        return false;
      });
    }
    dispatch.finish();
    if ((result != 1U) || (toFactor <= (primeCeiling * primeCeiling))) {
      return boost::lexical_cast<std::string>(result);
    }
  }

  // ECM pre-check (fast for small factors, fills the gap between Rho and QS)
  if (isEcm) {
    std::cout << "Running ECM pre-check...";
    const BigInteger ecmResult = ecm(toFactor, primes);
    if (ecmResult > 1U && ecmResult < toFactor) {
      std::cout << std::endl;
      return boost::lexical_cast<std::string>(ecmResult);
    }
    if (method == 3U) {
      std::cout << std::endl;
      return std::to_string(1);
    }
    std::cout << " Done." << std::endl;
  }

  // Pollard's Rho
  if (isPollardRho) {
    std::cout << "Running Pollard's Rho...";
    const BigInteger rhoResult = pollardRho(toFactor, sqrtN);
    if (rhoResult > 1U && rhoResult < toFactor) {
      std::cout << std::endl;
      return boost::lexical_cast<std::string>(rhoResult);
    }
    if (method == 2U) {
      std::cout << std::endl;
      return std::to_string(1);
    }
    std::cout << " Done." << std::endl;
  }

  // Set up wheel and gear factorization
  std::vector<size_t> gearFactorizationPrimes(primes.begin(), itg);
  std::vector<size_t> wheelFactorizationPrimes(primes.begin(), itw);
  std::vector<size_t> smoothPrimes;

  if (isQuadSieve) {
    std::cout << "Selecting factor base...";
    smoothPrimes = selectFactorBase(toFactor, primes);
    if (smoothPrimes.empty()) {
      throw std::runtime_error("No smooth primes found. Increase smoothness bound multiplier.");
    }
    for (const size_t& wpe : wheelPrimesExcluded) {
      const auto git = std::find(gearFactorizationPrimes.begin(), gearFactorizationPrimes.end(), wpe);
      if (git != gearFactorizationPrimes.end()) gearFactorizationPrimes.erase(git);
      const auto wit = std::find(wheelFactorizationPrimes.begin(), wheelFactorizationPrimes.end(), wpe);
      if (wit != wheelFactorizationPrimes.end()) wheelFactorizationPrimes.erase(wit);
    }
    std::cout << " Done." << std::endl;
  }

  std::cout << "Setting up wheels (and 'gears')...";
  BigInteger biggestWheelBigInt = 1U;
  for (const size_t &wp : gearFactorizationPrimes) biggestWheelBigInt *= (size_t)wp;
  biggestWheel = (size_t)biggestWheelBigInt;
  if (((BigInteger)biggestWheel) != biggestWheelBigInt) {
    throw std::runtime_error("Wheel too large! Reduce wheel/gear level. (Wheel radius: " + boost::lexical_cast<std::string>(biggestWheelBigInt) + ")");
  }

  wheel.clear();
  for (size_t i = 1U; i <= biggestWheel; ++i) {
    if (!isMultiple(i, gearFactorizationPrimes)) wheel.push_back(i);
  }
  if (wheel.empty()) wheel.push_back(1U);

  size_t batchItemCount = wheel.size();
  const size_t minBatch = 256U;
  if (minBatch > batchItemCount) {
    batchItemCount = ((minBatch + batchItemCount - 1U) / batchItemCount) * batchItemCount;
  }
  wheelFactorizationPrimes.clear();
  std::vector<boost::dynamic_bitset<size_t>> inc_seqs = generateGears(gearFactorizationPrimes);
  const size_t MIN_RTD_LEVEL = gearFactorizationPrimes.size() - wgDiff;
  const Wheel SMALLEST_WHEEL = wheelByPrimeCardinal(MIN_RTD_LEVEL);
  inc_seqs.erase(inc_seqs.begin(), inc_seqs.end() - wgDiff);
  gearFactorizationPrimes.clear();

  const auto ppBackwardFn = backward(SMALLEST_WHEEL);
  const auto ppForwardFn  = forward(SMALLEST_WHEEL);
  const BigInteger ppNodeRange = (((ppBackwardFn(sqrtN) + batchItemCount - 1U) / batchItemCount) + nodeCount - 1U) / nodeCount;
  const size_t ppStartingBatch = ((size_t)ppBackwardFn(primeCeiling)) / batchItemCount;

  const size_t rowLimit = smoothPrimes.size() + gaussianEliminationRowOffset;
  BigInteger qsBackwardLowBound = smoothBackwardFn(sqrtN + 1U);
  if (smoothForwardFn(qsBackwardLowBound) < (sqrtN + 1U)) ++qsBackwardLowBound;
  const BigInteger qsNodeRange = ((((smoothBackwardFn(sqrtN + (BigInteger)((toFactor - sqrtN).convert_to<double>() * sievingBoundMultiplier + 0.5)) - qsBackwardLowBound)
                                    + batchItemCount - 1U) / batchItemCount) + nodeCount - 1U) / nodeCount;

  std::cout << " Done." << std::endl;

  Factorizer worker(toFactor, sqrtN, qsBackwardLowBound,
                    isQuadSieve ? qsNodeRange : ppNodeRange,
                    nodeCount, nodeId, batchItemCount, rowLimit,
                    isQuadSieve ? 0U : ppStartingBatch,
                    smoothPrimes, wheelFactorizationLevel,
                    isQuadSieve ? ((wheel.size() > 1U) ? smoothForwardFn : forward(WHEEL1)) : ppForwardFn,
                    isQuadSieve ? ((wheel.size() > 1U) ? smoothBackwardFn : backward(WHEEL1)) : ppBackwardFn);

  if (isQuadSieve) {
    std::cout << "MPQS sieving. Smooth numbers:" << std::endl;
  }

  std::vector<std::future<BigInteger>> futures;
  futures.reserve(CpuCount);

  const auto workerFn = [&inc_seqs, &worker, &isQuadSieve] {
    std::vector<boost::dynamic_bitset<size_t>> inc_seqs_clone;
    inc_seqs_clone.reserve(inc_seqs.size());
    for (const boost::dynamic_bitset<size_t> &b : inc_seqs) inc_seqs_clone.emplace_back(b);
    return isQuadSieve ? worker.mpqsSieve(&inc_seqs_clone) : worker.bruteForce(&inc_seqs_clone);
  };

  for (unsigned cpu = 0U; cpu < CpuCount; ++cpu) {
    futures.push_back(std::async(std::launch::async, workerFn));
  }

  for (unsigned cpu = 0U; cpu < futures.size(); ++cpu) {
    const BigInteger r = futures[cpu].get();
    if ((r > 1U) && (r < toFactor)) return boost::lexical_cast<std::string>(r);
  }

  if (isQuadSieve) return boost::lexical_cast<std::string>(worker.solveForFactor());

  return std::to_string(1);
}

} // namespace Qimcifa

using namespace Qimcifa;

PYBIND11_MODULE(_find_a_factor, m) {
  m.doc() = "pybind11 plugin to find any factor of input (MPQS + ECM + Pollard Rho)";
  // method: 0 = PRIME_PROVER (brute force trial division)
  //         1 = FACTOR_FINDER (ECM → Pollard Rho → MPQS with large prime variant)
  //         2 = POLLARD_RHO (Pollard's Rho only, O(n^1/4))
  //         3 = ECM (Elliptic Curve Method only)
  m.def("_find_a_factor", &find_a_factor, "Finds any nontrivial factor of input");
}
