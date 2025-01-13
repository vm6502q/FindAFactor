///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2025. All rights reserved.
//
// "A quantum-inspired Monte Carlo integer factoring algorithm"
//
// This library was originally called ["Qimcifa"](https://github.com/vm6502q/qimcifa) and demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has
// since been developed into a general factoring algorithm and tool.
//
// `FindAFactor` uses heavily wheel-factorized brute-force "exhaust" numbers as "smooth" inputs to Quadratic Sieve, widely regarded as the asymptotically second fastest algorithm
// class known for cryptographically relevant semiprime factoring. `FindAFactor` is C++ based, with `pybind11`, which tends to make it faster than pure Python approaches. For the
// quick-and-dirty application of finding _any single_ nontrivial factor, something like at least 80% of positive integers will factorize in a fraction of a second, but the most
// interesting cases to consider are semiprime numbers, for which `FindAFactor` should be about as asymptotically competitive as similar Quadratic Sieve implementations.
//
// Our original contribution to Quadratic Sieve seems to be wheel factorization to 13 or 17 and maybe the idea of using the "exhaust" of a brute-force search for smooth number
// inputs for Quadratic Sieve. For wheel factorization (or "gear factorization"), we collect a short list of the first primes and remove all of their multiples from a "brute-force"
// guessing range by mapping a dense contiguous integer set, to a set without these multiples, relying on both a traditional "wheel," up to a middle prime number (of `11`), and a
// "gear-box" that stores increment values per prime according to the principles of wheel factorization, but operating semi-independently, to reduce space of storing the full
// wheel.
//
// Beyond this, we gain a functional advantage of a square-root over a more naive approach, by setting the brute force guessing range only between the highest prime in wheel
// factorization and the (modular) square root of the number to factor: if the number is semiprime, there is exactly one correct answer in this range, but including both factors in
// the range to search would cost us the square root advantage.
//
// Factoring this way is surprisingly easy to distribute: basically 0 network communication is needed to coordinate an arbitrarily high amount of parallelism to factor a single
// number. Each brute-force trial division instance is effectively 100% independent of all others (i.e. entirely "embarrassingly parallel"), and these guesses can seed independent
// Gaussian elimination matrices, so `FindAFactor` offers an extremely simply interface that allows work to be split between an arbitrarily high number of nodes with absolutely no
// network communication at all. In terms of incentives of those running different, cooperating nodes in the context of this specific number of integer factoring, all one
// ultimately cares about is knowing the correct factorization answer _by any means._ For pratical applications, there is no point at all in factoring a number whose factors are
// already known. When a hypothetical answer is forwarded to the (0-communication) "network" of collaborating nodes, _it is trivial to check whether the answer is correct_ (such as
// by simply entering the multiplication and equality check with the original number into a Python shell console)! Hence, collaborating node operators only need to trust that all
// participants in the "network" are actually performing their alloted segment of guesses and would actually communicate the correct answer to the entire group of collaborating
// nodes if any specific invidual happened to find the answer, but any purported answer is still trivial to verify.
//
//**Special thanks to OpenAI GPT "Elara," for indicated region of contributed code!**
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or
// https://www.gnu.org/licenses/lgpl-3.0.en.html for details.

#include "dispatchqueue.hpp"

#include <algorithm>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <stdlib.h>
#include <string>

#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Qimcifa {

typedef boost::multiprecision::cpp_int BigInteger;

const unsigned CpuCount = std::thread::hardware_concurrency();
DispatchQueue dispatch(CpuCount);

enum Wheel { ERROR = 0, WHEEL1 = 1, WHEEL2 = 2, WHEEL3 = 6, WHEEL5 = 30, WHEEL7 = 210, WHEEL11 = 2310 };

Wheel wheelByPrimeCardinal(int i) {
  switch (i) {
  case 0:
    return WHEEL1;
  case 1:
    return WHEEL2;
  case 2:
    return WHEEL3;
  case 3:
    return WHEEL5;
  case 4:
    return WHEEL7;
  case 5:
    return WHEEL11;
  default:
    return ERROR;
  }
}

// See https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
BigInteger ipow(BigInteger base, unsigned exp) {
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
  // Otherwise, find b = sqrt(b^2).
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
  // Otherwise, find b = sqrt(b^2).
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

// We are multiplying out the first distinct primes, below.

// Make this NOT a multiple of 2.
inline size_t forward2(const size_t &p) { return (p << 1U) | 1U; }

inline size_t backward2(const size_t &p) { return (size_t)(p >> 1U); }

// Make this NOT a multiple of 2 or 3.
inline size_t forward3(const size_t &p) { return (p << 1U) + (~(~p | 1U)) - 1U; }

inline size_t backward3(const size_t &n) { return (size_t)((~(~n | 1U)) / 3U) + 1U; }

constexpr unsigned char wheel5[8U] = {1U, 7U, 11U, 13U, 17U, 19U, 23U, 29U};

// Make this NOT a multiple of 2, 3, or 5.
size_t forward5(const size_t &p) { return wheel5[p & 7U] + (p >> 3U) * 30U; }

size_t backward5(const size_t &n) { return std::distance(wheel5, std::lower_bound(wheel5, wheel5 + 8U, (size_t)(n % 30U))) + 8U * (size_t)(n / 30U) + 1U; }

constexpr unsigned char wheel7[48U] = {1U,   11U,  13U,  17U,  19U,  23U,  29U,  31U,  37U,  41U,  43U,  47U,  53U,  59U,  61U,  67U,
                                       71U,  73U,  79U,  83U,  89U,  97U,  101U, 103U, 107U, 109U, 113U, 121U, 127U, 131U, 137U, 139U,
                                       143U, 149U, 151U, 157U, 163U, 167U, 169U, 173U, 179U, 181U, 187U, 191U, 193U, 197U, 199U, 209U};

// Make this NOT a multiple of 2, 3, 5, or 7.
size_t forward7(const size_t &p) { return wheel7[p % 48U] + (p / 48U) * 210U; }

size_t backward7(const size_t &n) { return std::distance(wheel7, std::lower_bound(wheel7, wheel7 + 48U, (size_t)(n % 210U))) + 48U * (size_t)(n / 210U) + 1U; }

constexpr unsigned short wheel11[480U] = {
    1U,    13U,   17U,   19U,   23U,   29U,   31U,   37U,   41U,   43U,   47U,   53U,   59U,   61U,   67U,   71U,   73U,   79U,   83U,   89U,   97U,   101U,  103U,  107U,
    109U,  113U,  127U,  131U,  137U,  139U,  149U,  151U,  157U,  163U,  167U,  169U,  173U,  179U,  181U,  191U,  193U,  197U,  199U,  211U,  221U,  223U,  227U,  229U,
    233U,  239U,  241U,  247U,  251U,  257U,  263U,  269U,  271U,  277U,  281U,  283U,  289U,  293U,  299U,  307U,  311U,  313U,  317U,  323U,  331U,  337U,  347U,  349U,
    353U,  359U,  361U,  367U,  373U,  377U,  379U,  383U,  389U,  391U,  397U,  401U,  403U,  409U,  419U,  421U,  431U,  433U,  437U,  439U,  443U,  449U,  457U,  461U,
    463U,  467U,  479U,  481U,  487U,  491U,  493U,  499U,  503U,  509U,  521U,  523U,  527U,  529U,  533U,  541U,  547U,  551U,  557U,  559U,  563U,  569U,  571U,  577U,
    587U,  589U,  593U,  599U,  601U,  607U,  611U,  613U,  617U,  619U,  629U,  631U,  641U,  643U,  647U,  653U,  659U,  661U,  667U,  673U,  677U,  683U,  689U,  691U,
    697U,  701U,  703U,  709U,  713U,  719U,  727U,  731U,  733U,  739U,  743U,  751U,  757U,  761U,  767U,  769U,  773U,  779U,  787U,  793U,  797U,  799U,  809U,  811U,
    817U,  821U,  823U,  827U,  829U,  839U,  841U,  851U,  853U,  857U,  859U,  863U,  871U,  877U,  881U,  883U,  887U,  893U,  899U,  901U,  907U,  911U,  919U,  923U,
    929U,  937U,  941U,  943U,  947U,  949U,  953U,  961U,  967U,  971U,  977U,  983U,  989U,  991U,  997U,  1003U, 1007U, 1009U, 1013U, 1019U, 1021U, 1027U, 1031U, 1033U,
    1037U, 1039U, 1049U, 1051U, 1061U, 1063U, 1069U, 1073U, 1079U, 1081U, 1087U, 1091U, 1093U, 1097U, 1103U, 1109U, 1117U, 1121U, 1123U, 1129U, 1139U, 1147U, 1151U, 1153U,
    1157U, 1159U, 1163U, 1171U, 1181U, 1187U, 1189U, 1193U, 1201U, 1207U, 1213U, 1217U, 1219U, 1223U, 1229U, 1231U, 1237U, 1241U, 1247U, 1249U, 1259U, 1261U, 1271U, 1273U,
    1277U, 1279U, 1283U, 1289U, 1291U, 1297U, 1301U, 1303U, 1307U, 1313U, 1319U, 1321U, 1327U, 1333U, 1339U, 1343U, 1349U, 1357U, 1361U, 1363U, 1367U, 1369U, 1373U, 1381U,
    1387U, 1391U, 1399U, 1403U, 1409U, 1411U, 1417U, 1423U, 1427U, 1429U, 1433U, 1439U, 1447U, 1451U, 1453U, 1457U, 1459U, 1469U, 1471U, 1481U, 1483U, 1487U, 1489U, 1493U,
    1499U, 1501U, 1511U, 1513U, 1517U, 1523U, 1531U, 1537U, 1541U, 1543U, 1549U, 1553U, 1559U, 1567U, 1571U, 1577U, 1579U, 1583U, 1591U, 1597U, 1601U, 1607U, 1609U, 1613U,
    1619U, 1621U, 1627U, 1633U, 1637U, 1643U, 1649U, 1651U, 1657U, 1663U, 1667U, 1669U, 1679U, 1681U, 1691U, 1693U, 1697U, 1699U, 1703U, 1709U, 1711U, 1717U, 1721U, 1723U,
    1733U, 1739U, 1741U, 1747U, 1751U, 1753U, 1759U, 1763U, 1769U, 1777U, 1781U, 1783U, 1787U, 1789U, 1801U, 1807U, 1811U, 1817U, 1819U, 1823U, 1829U, 1831U, 1843U, 1847U,
    1849U, 1853U, 1861U, 1867U, 1871U, 1873U, 1877U, 1879U, 1889U, 1891U, 1901U, 1907U, 1909U, 1913U, 1919U, 1921U, 1927U, 1931U, 1933U, 1937U, 1943U, 1949U, 1951U, 1957U,
    1961U, 1963U, 1973U, 1979U, 1987U, 1993U, 1997U, 1999U, 2003U, 2011U, 2017U, 2021U, 2027U, 2029U, 2033U, 2039U, 2041U, 2047U, 2053U, 2059U, 2063U, 2069U, 2071U, 2077U,
    2081U, 2083U, 2087U, 2089U, 2099U, 2111U, 2113U, 2117U, 2119U, 2129U, 2131U, 2137U, 2141U, 2143U, 2147U, 2153U, 2159U, 2161U, 2171U, 2173U, 2179U, 2183U, 2197U, 2201U,
    2203U, 2207U, 2209U, 2213U, 2221U, 2227U, 2231U, 2237U, 2239U, 2243U, 2249U, 2251U, 2257U, 2263U, 2267U, 2269U, 2273U, 2279U, 2281U, 2287U, 2291U, 2293U, 2297U, 2309U};

// Make this NOT a multiple of 2, 3, 5, 7, or 11.
size_t forward11(const size_t &p) { return wheel11[p % 480U] + (p / 480U) * 2310U; }

size_t backward11(const size_t &n) { return std::distance(wheel11, std::lower_bound(wheel11, wheel11 + 480U, (size_t)(n % 2310U))) + 480U * (size_t)(n / 2310U) + 1U; }

inline BigInteger _forward2(const BigInteger &p) { return (p << 1U) | 1U; }

inline BigInteger _backward2(const BigInteger &n) { return n >> 1U; }

inline BigInteger _forward3(const BigInteger &p) { return (p << 1U) + (~(~p | 1U)) - 1U; }

inline BigInteger _backward3(const BigInteger &n) { return ((~(~n | 1U)) / 3U) + 1U; }

BigInteger _forward5(const BigInteger &p) { return wheel5[(size_t)(p & 7U)] + (p >> 3U) * 30U; }

BigInteger _backward5(const BigInteger &n) { return std::distance(wheel5, std::lower_bound(wheel5, wheel5 + 8U, (size_t)(n % 30U))) + 8U * (n / 30U) + 1U; }

BigInteger _forward7(const BigInteger &p) { return wheel7[(size_t)(p % 48U)] + (p / 48U) * 210U; }

BigInteger _backward7(const BigInteger &n) { return std::distance(wheel7, std::lower_bound(wheel7, wheel7 + 48U, n % 210U)) + 48U * (n / 210U) + 1U; }

BigInteger _forward11(const BigInteger &p) { return wheel11[(size_t)(p % 480U)] + (p / 480U) * 2310U; }

BigInteger _backward11(const BigInteger &n) { return std::distance(wheel11, std::lower_bound(wheel11, wheel11 + 480U, (size_t)(n % 2310U))) + 480U * (n / 2310U) + 1U; }

typedef BigInteger (*ForwardFn)(const BigInteger &);
inline ForwardFn forward(const Wheel &w) {
  switch (w) {
  case WHEEL2:
    return _forward2;
  case WHEEL3:
    return _forward3;
  case WHEEL5:
    return _forward5;
  case WHEEL7:
    return _forward7;
  case WHEEL11:
    return _forward11;
  case WHEEL1:
  default:
    return [](const BigInteger &n) -> BigInteger { return n; };
  }
}

inline ForwardFn backward(const Wheel &w) {
  switch (w) {
  case WHEEL2:
    return _backward2;
  case WHEEL3:
    return _backward3;
  case WHEEL5:
    return _backward5;
  case WHEEL7:
    return _backward7;
  case WHEEL11:
    return _backward11;
  case WHEEL1:
  default:
    return [](const BigInteger &n) -> BigInteger { return n; };
  }
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

  knownPrimes.reserve(std::expint(log((double)n)) - std::expint(log(2)));

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

boost::dynamic_bitset<size_t> wheel_inc(std::vector<size_t> primes) {
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

std::vector<boost::dynamic_bitset<size_t>> wheel_gen(const std::vector<size_t> &primes) {
  std::vector<boost::dynamic_bitset<size_t>> output;
  std::vector<size_t> wheelPrimes;
  for (const size_t &p : primes) {
    wheelPrimes.push_back(p);
    output.push_back(wheel_inc(wheelPrimes));
  }

  return output;
}

size_t GetWheelIncrement(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
  size_t wheelIncrement = 0U;
  bool is_wheel_multiple = false;
  do {
    for (size_t i = 0U; i < inc_seqs->size(); ++i) {
      boost::dynamic_bitset<size_t> &wheel = (*inc_seqs)[i];
      is_wheel_multiple = wheel.test(0U);
      wheel >>= 1U;
      if (is_wheel_multiple) {
        wheel[wheel.size() - 1U] = true;
        break;
      }
    }
    ++wheelIncrement;
  } while (is_wheel_multiple);

  return wheelIncrement;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  WRITTEN WITH ELARA (GPT) BELOW                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Utility to perform modular exponentiation
inline BigInteger modExp(BigInteger base, BigInteger exp, const BigInteger &mod) {
  BigInteger result = 1U;
  while (exp) {
    if (exp & 1U) {
      result = (result * base) % mod;
    }
    base = (base * base) % mod;
    exp >>= 1U;
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  WRITTEN WITH ELARA (GPT) ABOVE                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Factorizer {
  std::mutex batchMutex;
  std::mutex smoothNumberMapMutex;
  std::default_random_engine rng;
  std::mt19937_64 gen;
  std::uniform_int_distribution<size_t> dis;
  BigInteger toFactorSqr;
  BigInteger toFactor;
  BigInteger toFactorSqrt;
  BigInteger batchRange;
  BigInteger batchNumber;
  BigInteger batchOffset;
  BigInteger batchTotal;
  BigInteger wheelRadius;
  size_t wheelEntryCount;
  size_t smoothPartsLimit;
  size_t rowOffset;
  bool isIncomplete;
  std::vector<size_t> primes;
  ForwardFn forwardFn;
  std::vector<BigInteger> smoothNumberKeys;
  std::vector<boost::dynamic_bitset<size_t>> smoothNumberValues;

  Factorizer(const BigInteger &tfsqr, const BigInteger &tf, const BigInteger &tfsqrt, const BigInteger &range, size_t nodeCount, size_t nodeId, size_t w, size_t spl,
             const std::vector<size_t> &p, ForwardFn fn)
    : rng({}), gen(rng()), dis(0U, p.size() - 1U), toFactorSqr(tfsqr), toFactor(tf), toFactorSqrt(tfsqrt), batchRange(range), batchNumber(0U), batchOffset(nodeId * range), batchTotal(nodeCount * range),
    wheelRadius(1U), wheelEntryCount(w), smoothPartsLimit(spl), rowOffset(p.size()), isIncomplete(true), primes(p), forwardFn(fn)
  {
    for (size_t i = 0U; i < primes.size(); ++i) {
      const size_t& p = primes[i];
      wheelRadius *= p;
      smoothNumberKeys.push_back(p);
      smoothNumberValues.emplace_back(primes.size(), 0);
      smoothNumberValues.back()[i] = true;
    }
  }

  BigInteger getNextAltBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);

    if (batchNumber >= batchRange) {
      isIncomplete = false;
    }

    const BigInteger halfIndex = batchOffset + (batchNumber++ >> 1U) + 1U;

    return ((batchNumber & 1U) ? batchTotal - halfIndex : halfIndex);
  }

  BigInteger bruteForce(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs) {
    // Up to wheel factorization, try all batches up to the square root of toFactor.
    for (BigInteger batchNum = getNextAltBatch(); isIncomplete; batchNum = getNextAltBatch()) {
      const BigInteger batchStart = batchNum * wheelEntryCount;
      const BigInteger batchEnd = batchStart + wheelEntryCount;
      for (BigInteger p = batchStart; p < batchEnd;) {
        const BigInteger n = forwardFn(p);
        if (!(toFactor % n) && (n != 1U) && (n != toFactor)) {
          isIncomplete = false;
          return n;
        }
        p += GetWheelIncrement(inc_seqs);
      }
    }

    return 1U;
  }

  BigInteger smoothCongruences(std::vector<boost::dynamic_bitset<size_t>> *inc_seqs, std::vector<BigInteger> *semiSmoothParts, bool isGaussElim) {
    // Up to wheel factorization, try all batches up to the square root of toFactor.
    // Since the largest prime factors of these numbers is relatively small,
    // use the "exhaust" of brute force to produce smooth numbers for Quadratic Sieve.
    for (BigInteger batchNum = getNextAltBatch(); isIncomplete; batchNum = getNextAltBatch()) {
      const BigInteger batchStart = batchNum * wheelEntryCount;
      const BigInteger batchEnd = batchStart + wheelEntryCount;
      for (BigInteger p = batchStart; p < batchEnd;) {
        // Brute-force check if the sequential number is a factor.
        const BigInteger n = forwardFn(p);
        // If so, terminate this node and return the answer.
        if (!(toFactor % n) && (n != 1U) && (n != toFactor)) {
          isIncomplete = false;
          return n;
        }
        // Use the "exhaust" to produce smoother numbers.
        semiSmoothParts->push_back(n);
        // Skip increments on the "wheels" (or "gears").
        p += GetWheelIncrement(inc_seqs);
      }

      // Batch this work, to reduce contention.
      if (semiSmoothParts->size() >= smoothPartsLimit) {
        makeSmoothNumbers(semiSmoothParts, isGaussElim);

        return 1U;
      }
    }

    return 1U;
  }

  // Compute the prime factorization modulo 2
  boost::dynamic_bitset<size_t> factorizationVector(BigInteger num) {
    boost::dynamic_bitset<size_t> vec(primes.size(), 0);
    while (true) {
      BigInteger factor = gcd(num, wheelRadius);
      if (factor == 1U) {
        break;
      }
      num /= factor;
      // Remove smooth primes from factor
      for (size_t pi = 0U; pi < primes.size(); ++pi) {
        const size_t& p = primes[pi];
        if (factor % p) {
          continue;
        }
        factor /= p;
        vec.flip(pi);
        if (factor == 1U) {
          break;
        }
      }
      if (num == 1U) {
        return vec;
      }
    }
    if (num != 1U) {
      return boost::dynamic_bitset<size_t>();
    }

    return vec;
  }

  void makeSmoothNumbers(std::vector<BigInteger> *semiSmoothParts, bool isGaussElim) {
    // Factorize all "smooth parts."
    std::vector<BigInteger> smoothParts;
    std::map<BigInteger, boost::dynamic_bitset<size_t>> smoothPartsMap;
    for (const BigInteger &n : (*semiSmoothParts)) {
      const boost::dynamic_bitset<size_t> fv = factorizationVector(n);
      if (fv.size()) {
        smoothPartsMap[n] = fv;
        smoothParts.push_back(n);
      }
    }
    // We can clear the thread's buffer vector.
    semiSmoothParts->clear();

    // This is the only nondeterminism in the algorithm.
    std::shuffle(smoothParts.begin(), smoothParts.end(), rng);

    const BigInteger limit = isGaussElim ? toFactor : toFactorSqrt;

    // Now that smooth parts have been shuffled, just multiply down the list until they are larger than square root of toFactor.
    BigInteger smoothNumber = 1U;
    boost::dynamic_bitset<size_t> fv(primes.size(), 0);
    for (size_t spi = 0U; spi < smoothParts.size(); ++spi) {
      const BigInteger &sp = smoothParts[spi];
      // This multiplies together the factorizations of the smooth parts
      // (producing the overall factorization of their multiplication)
      fv ^= smoothPartsMap[sp];
      smoothNumber *= sp;
      // Check if the number is big enough
      if (smoothNumber <= limit) {
        continue;
      }
      if (true) {
        std::lock_guard<std::mutex> lock(smoothNumberMapMutex);
        smoothNumberValues.emplace_back(fv);
        smoothNumberKeys.push_back(smoothNumber);
      }
      // Reset "smoothNumber" and its factorization vector.
      smoothNumber = 1U;
      fv = boost::dynamic_bitset<size_t>(primes.size(), 0);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                   WRITTEN WITH ELARA (GPT) BELOW                                                      //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Perform Gaussian elimination on a binary matrix
  void gaussianElimination() {
    const unsigned cpuCount = CpuCount;
    auto mColIt = smoothNumberValues.begin();
    auto nColIt = smoothNumberKeys.begin();
    const size_t rows = smoothNumberValues.size();
    for (size_t col = 0U; col < primes.size(); ++col) {
      auto mRowIt = mColIt;
      auto nRowIt = nColIt;

      int64_t pivot = -1;
      for (size_t row = col; row < rows; ++row) {
        if ((*mRowIt)[col]) {
          // Swapping matrix rows corresponds
          // with swapping factorized numbers.
          if (row != col) {
            std::swap(*mColIt, *mRowIt);
            std::swap(*nColIt, *nRowIt);
          }
          pivot = row;
          break;
        }
        ++nRowIt;
        ++mRowIt;
      }

      if (pivot != -1) {
        const boost::dynamic_bitset<size_t> &cm = *mColIt;
        const BigInteger &cn = *nColIt;
        mRowIt = smoothNumberValues.begin();
        nRowIt = smoothNumberKeys.begin();
        for (unsigned cpu = 0U; (cpu < CpuCount) && (cpu < rows); ++cpu) {
          dispatch.dispatch([cpu, &cpuCount, &col, &rows, &cm, &cn, nRowIt, mRowIt]() -> bool {
            auto mrIt = mRowIt;
            auto nrIt = nRowIt;
            for (size_t row = cpu; ; row += cpuCount) {
              boost::dynamic_bitset<size_t> &rm = *mrIt;
              BigInteger &rn = *nrIt;
              if ((row != col) && rm[col]) {
                // XOR-ing factorization rows
                // is like multiplying the numbers.
                rm ^= cm;
                rn *= cn;
              }
              if ((row + cpuCount) >= rows) {
                return false;
              }
              std::advance(nrIt, cpuCount);
              std::advance(mrIt, cpuCount);
            }

            return false;
          });
          ++mRowIt;
          ++nRowIt;
        }
        dispatch.finish();
      }

      ++mColIt;
      ++nColIt;
    }
  }

  BigInteger checkPerfectSquare(BigInteger perfectSquare) {
    // Compute x and y
    const BigInteger x = perfectSquare % toFactor;
    const BigInteger y = modExp(x, toFactor >> 1U, toFactor);

    // Check congruence of squares
    BigInteger factor = gcd(toFactor, x + y);
    if ((factor != 1U) && (factor != toFactor)) {
      return factor;
    }

    if (x == y) {
      return 1U;
    }

    // Try x - y as well
    factor = gcd(toFactor, x - y);
    if ((factor != 1U) && (factor != toFactor)) {
      return factor;
    }

    return 1U;
  }

  // Find duplicate rows
  BigInteger findDuplicateRows(const BigInteger &target) {
    // Check for linear dependencies and find a congruence of squares
    std::mutex rowMutex;
    BigInteger result = 1U;
    std::set<size_t> toStrike;
    auto iIt = smoothNumberValues.begin();
    const size_t rowCount = smoothNumberValues.size();
    const size_t rowCountMin1 = rowCount - 1U;
    for (size_t i = primes.size(); (i < rowCountMin1) && (result == 1U); ++i) {
      dispatch.dispatch([this, &target, i, iIt, &rowCount, &result, &rowMutex, &toStrike]() -> bool {
        boost::dynamic_bitset<size_t> &iRow = *iIt;
        const BigInteger& iInt = this->smoothNumberKeys[i];

        const size_t startJ = std::max(this->rowOffset, i + 1U);
        auto jIt = this->smoothNumberValues.begin();
        std::advance(jIt, (startJ - 1U));
        for (size_t j = startJ; j < rowCount; ++j) {
          ++jIt;

          boost::dynamic_bitset<size_t> &jRow = *jIt;
          if (iRow != jRow) {
            continue;
          }

          const BigInteger& jInt = this->smoothNumberKeys[j];
          if (iInt < jInt) {
            std::lock_guard<std::mutex> lock(rowMutex);
            toStrike.insert(j);
          } else {
            std::lock_guard<std::mutex> lock(rowMutex);
            toStrike.insert(i);
          }

          const BigInteger factor = checkPerfectSquare(this->smoothNumberKeys[i]);
          if ((factor != 1U) && (factor != target)) {
            std::lock_guard<std::mutex> lock(rowMutex);
            result = factor;

            return true;
          }
        }

        return false;
      });
      ++iIt;
    }
    dispatch.finish();

    if (result != 1U) {
      return result;
    }

    // These numbers have been tried already:
    for (const size_t& i : toStrike) {
      smoothNumberKeys.erase(smoothNumberKeys.begin() + i);
      smoothNumberValues.erase(smoothNumberValues.begin() + i);
    }

    rowOffset = smoothNumberKeys.size();

    return 1U; // No factor found
  }

  // Use Gaussian elimination
  BigInteger findFactor(const BigInteger &target) {
    // Gaussian elimination multiplies these numbers
    // with small primes, to produce squares
    gaussianElimination();

    // Check for linear dependencies and find a congruence of squares
    std::mutex rowMutex;
    BigInteger result = 1U;
    const size_t rowCount = smoothNumberKeys.size();
    for (size_t i = primes.size(); (i < rowCount) && (result == 1U); ++i) {
      dispatch.dispatch([this, &target, i, &result, &rowMutex]() -> bool {
        const BigInteger factor = checkPerfectSquare(this->smoothNumberKeys[i]);

        if ((factor != 1U) && (factor != target)) {
          std::lock_guard<std::mutex> lock(rowMutex);
          result = factor;

          return true;
        }

        return false;
      });
    }
    dispatch.finish();

    if (result != 1U) {
      return result;
    }

    // These numbers have been tried already:
    smoothNumberKeys.resize(primes.size());
    smoothNumberValues.resize(primes.size());

    return 1U; // No factor found
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                   WRITTEN WITH ELARA (GPT) ABOVE                                                      //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};

std::string find_a_factor(const std::string &toFactorStr, const bool &isConOfSqr, const bool &isGaussElim, const size_t &nodeCount, const size_t &nodeId,
                          size_t trialDivisionLevel, size_t gearFactorizationLevel, size_t wheelFactorizationLevel, double smoothnessBoundMultiplier, double batchSizeMultiplier) {
  // (At least) level 11 wheel factorization is baked into basic functions.
  if (!wheelFactorizationLevel) {
    wheelFactorizationLevel = 1U;
  } else if (wheelFactorizationLevel > 11U) {
    wheelFactorizationLevel = 11U;
    std::cout << "Warning: Wheel factorization limit is 11. (Parameter will be ignored and default to 11.)" << std::endl;
  }
  if (!gearFactorizationLevel) {
    gearFactorizationLevel = 1U;
  } else if (gearFactorizationLevel < wheelFactorizationLevel) {
    gearFactorizationLevel = wheelFactorizationLevel;
    std::cout << "Warning: Gear factorization level must be at least as high as wheel level. (Parameter will be ignored and default to wheel level.)" << std::endl;
  }

  // Convert from string.
  const BigInteger toFactor(toFactorStr);

  // The largest possible discrete factor of "toFactor" is its square root (as with any integer).
  const BigInteger fullMaxBase = sqrt(toFactor);
  if (fullMaxBase * fullMaxBase == toFactor) {
    return boost::lexical_cast<std::string>(fullMaxBase);
  }

  // We only need to try trial division about as high as would be necessary for 4096 bits of semiprime.
  const size_t primeCeiling = (trialDivisionLevel < fullMaxBase) ? trialDivisionLevel : (size_t)fullMaxBase;
  BigInteger result = 1U;
  // This uses very little memory and time, to find primes.
  std::vector<size_t> primes = SieveOfEratosthenes(primeCeiling);
  // "it" is the end-of-list iterator for a list up-to-and-including wheelFactorizationLevel.
  const auto itw = std::upper_bound(primes.begin(), primes.end(), wheelFactorizationLevel);
  const auto itg = std::upper_bound(primes.begin(), primes.end(), gearFactorizationLevel);
  const size_t wgDiff = std::distance(itw, itg);

  // This is simply trial division up to the ceiling.
  std::mutex trialDivisionMutex;
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

  // Set up wheel factorization (or "gear" factorization)
  std::vector<size_t> gearFactorizationPrimes(primes.begin(), itg);
  std::vector<size_t> wheelFactorizationPrimes(primes.begin(), itw);
  // Keep as many "smooth" primes as bits in number to factor.
  const size_t toFactorBits = (size_t)log2(toFactor);
  size_t smoothPrimeCount = (size_t)(smoothnessBoundMultiplier * toFactorBits);
  if (!smoothPrimeCount) {
    smoothPrimeCount = 1U;
    std::cout << "Warning: smoothness bound multiplier would retain no primes, but it must retain at least 1. (Defaulting to retaining 1 prime.)" << std::endl;
  }
  // Primes are only present in range above wheel factorization level
  primes.erase(primes.begin(), itg);
  const size_t maxPrimeCount = std::min(primes.size(), smoothPrimeCount);
  std::vector<size_t> smoothPrimes;
  for (size_t primeId = 0U; (primeId < primes.size()) && (smoothPrimes.size() < maxPrimeCount); ++primeId) {
    const size_t p = primes[primeId];
    const size_t residue = (size_t)(toFactor % p);
    const size_t sr = _sqrt(residue);
    if ((sr * sr) == residue) {
      smoothPrimes.push_back(p);
    }
  }
  if (isConOfSqr && (smoothPrimes.size() < maxPrimeCount)) {
    std::cout << "Warning: Factor base truncated to " << smoothPrimes.size() << " factors. If you don't want to truncate, set the trial division level option higher." << std::endl;
  }
  // From 1, this is a period for wheel factorization
  size_t biggestWheel = 1ULL;
  for (const size_t &wp : gearFactorizationPrimes) {
    biggestWheel *= (size_t)wp;
  }
  // Wheel entry count per largest "gear" scales our brute-force range.
  size_t wheelEntryCount = 0U;
  for (size_t i = 0U; i < biggestWheel; ++i) {
    if (!isMultiple(i, wheelFactorizationPrimes)) {
      ++wheelEntryCount;
    }
  }
  wheelFactorizationPrimes.clear();
  // These are "gears," for wheel factorization (with a "wheel" already in place up to 11).
  std::vector<boost::dynamic_bitset<size_t>> inc_seqs = wheel_gen(gearFactorizationPrimes);
  // We're done with the lowest primes.
  const size_t MIN_RTD_LEVEL = gearFactorizationPrimes.size() - wgDiff;
  const Wheel SMALLEST_WHEEL = wheelByPrimeCardinal(MIN_RTD_LEVEL);
  // Skip multiples removed by wheel factorization.
  inc_seqs.erase(inc_seqs.begin(), inc_seqs.end() - wgDiff);
  gearFactorizationPrimes.clear();

  // Range per parallel node
  const BigInteger nodeRange = (((backward(SMALLEST_WHEEL)(fullMaxBase) + nodeCount - 1U) / nodeCount) + wheelEntryCount - 1U) / wheelEntryCount;
  // This manages the work of all threads.
  Factorizer worker(toFactor * toFactor, toFactor, fullMaxBase,
                    nodeRange, nodeCount, nodeId,
                    wheelEntryCount, (size_t)((wheelEntryCount << 1U) * batchSizeMultiplier),
                    smoothPrimes, forward(SMALLEST_WHEEL));

  if (!isConOfSqr) {
    const auto workerFn = [&inc_seqs, &worker] {
      // inc_seq needs to be independent per thread.
      std::vector<boost::dynamic_bitset<size_t>> inc_seqs_clone;
      inc_seqs_clone.reserve(inc_seqs.size());
      for (const boost::dynamic_bitset<size_t> &b : inc_seqs) {
        inc_seqs_clone.emplace_back(b);
      }

      // "Brute force" includes extensive wheel multiplication and can be faster.
      return worker.bruteForce(&inc_seqs_clone);
    };

    std::vector<std::future<BigInteger>> futures;
    futures.reserve(CpuCount);

    for (unsigned cpu = 0U; cpu < CpuCount; ++cpu) {
      futures.push_back(std::async(std::launch::async, workerFn));
    }

    for (unsigned cpu = 0U; cpu < futures.size(); ++cpu) {
      const BigInteger r = futures[cpu].get();
      if ((r > result) && (r != toFactor)) {
        result = r;
      }
    }

    return boost::lexical_cast<std::string>(result);
  }

  const auto smoothNumberFn = [&inc_seqs, &wheelEntryCount, &batchSizeMultiplier, &worker, &isGaussElim] {
    // inc_seq needs to be independent per thread.
    std::vector<boost::dynamic_bitset<size_t>> inc_seqs_clone;
    inc_seqs_clone.reserve(inc_seqs.size());
    for (const boost::dynamic_bitset<size_t> &b : inc_seqs) {
      inc_seqs_clone.emplace_back(b);
    }

    // Different collections per thread;
    std::vector<BigInteger> semiSmoothParts;
    semiSmoothParts.reserve((size_t)((wheelEntryCount << 1U) * batchSizeMultiplier));

    // While brute-forcing, use the "exhaust" to feed "smooth" number generation and check conguence of squares.
    return worker.smoothCongruences(&inc_seqs_clone, &semiSmoothParts, isGaussElim);
  };

  std::vector<std::future<BigInteger>> futures;
  futures.reserve(CpuCount);

  do {
    for (unsigned cpu = 0U; cpu < CpuCount; ++cpu) {
      futures.push_back(std::async(std::launch::async, smoothNumberFn));
    }

    for (unsigned cpu = 0U; cpu < futures.size(); ++cpu) {
      const BigInteger r = futures[cpu].get();
      if ((r > result) && (r != toFactor)) {
        result = r;
      }
    }

    if ((result != 1U) && (result != toFactor)) {
      return boost::lexical_cast<std::string>(result);
    }

    futures.clear();

    // This next section is for (Quadratic Sieve) Gaussian elimination.
    result = isGaussElim ? worker.findFactor(toFactor) : worker.findDuplicateRows(toFactor);
  } while ((result == 1U) || (result == toFactor));

  return boost::lexical_cast<std::string>(result);
}
} // namespace Qimcifa

using namespace Qimcifa;

PYBIND11_MODULE(_find_a_factor, m) {
  m.doc() = "pybind11 plugin to find any factor of input";
  m.def("_find_a_factor", &find_a_factor, "Finds any nontrivial factor of input (or returns 1 if prime)");
}