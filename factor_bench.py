import sys
import time

from FindAFactor import find_a_factor

argv_len = len(sys.argv)

if argv_len < 2:
    print("Enter a number to factor on the command line when calling " + sys.argv[0])

to_factor = int(sys.argv[1])
node_count = 1
node_id = 0
wheel_factorization_level = 17

if argv_len > 3:
    node_count=int(sys.argv[2])
    node_id=int(sys.argv[3])
if argv_len > 4:
    wheel_factorization_level=int(sys.argv[4])

start = time.perf_counter()
print(find_a_factor(
    to_factor,
    wheel_factorization_level=wheel_factorization_level,
    node_count=node_count,
    node_id=node_id
))
print(time.perf_counter() - start)
