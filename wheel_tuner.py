# Choose a wheel factorization ceiling:
wheel_factorization_level = 50
# Grab a sample of smooth numbers from `FACTOR_FINDER` mode, and paste them here.
smooth_nums = [42571, 42576, 42583, 42589, 42591, 42592, 42593, 42595, 42598, 42601, 42612, 42615, 42617, 42634, 42635, 42639, 42641, 42849, 43079, 42851, 43094, 42852, 43105, 43109, 42864, 43114, 42871, 43133, 43150, 42909, 43168, 42929, 42931, 42936, 43187, 42937, 42942, 43209, 42962, 42978, 42985, 42986, 42650, 43247, 42677, 43279, 43280, 43285, 42994, 43366, 43314, 43013, 43397, 43594, 43409, 43596, 43846, 43855, 42687, 43059, 43864, 44109, 43671, 43423, 42689,]

# Primes up to 1000 are likely enough to choose from in general (but you can add more).
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]
small_primes = [p for p in primes if p < wheel_factorization_level]
smooth_num_count = len(smooth_nums)

for p in small_primes:
    mult_count = 0
    for x in smooth_nums:
        if (x % p) == 0:
            mult_count += 1

    print("Prime = " + str(p) + ", quality = " + (str((1 / p) / ((mult_count / smooth_num_count))) if mult_count else "MAX"))

print("Any low quality score (particularly less than 1.0, but even higher) might be worth excluding.")
print("MAX quality scores are the best to include.")