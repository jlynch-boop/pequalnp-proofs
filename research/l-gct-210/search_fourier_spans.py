import argparse
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product

import numpy as np


PARTITIONS = [(5,), (4, 1), (3, 2), (3, 1, 1), (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)]
PERMUTATION_CELLS = np.array([[5*i + p[i] for i in range(5)] for p in permutations(range(5))])
PRIME = 61


def cycle_vectors(partition):
    generator = next(g for g in range(2, PRIME)
                     if all(pow(g, (PRIME-1)//factor, PRIME) != 1 for factor in (2, 3, 5)))
    start = 0
    for length in partition:
        root = pow(generator, (PRIME-1)//length, PRIME)
        for frequency in range(length):
            vector = [0] * 5
            vector[start:start+length] = [pow(root, frequency*k, PRIME) for k in range(length)]
            yield vector
        start += length


@cache
def positive_simplex(size):
    return np.array([xs for xs in product(range(1, 6), repeat=size) if sum(xs) <= 5],
                    dtype=np.int64)


def search(row_partition, column_partition, target):
    modes = np.array([np.outer(x, y).reshape(25) % PRIME
                      for x in cycle_vectors(row_partition) for y in cycle_vectors(column_partition)])

    @cache
    def passes(indices):
        directions = modes[list(indices)]
        inputs = (1 + positive_simplex(len(indices)) @ directions) % PRIME
        values = np.prod(inputs[:, PERMUTATION_CELLS], axis=-1, dtype=np.int64).sum(axis=1) % PRIME
        return bool(np.all(values == 120 % PRIME))

    vertices = [i for i in range(25) if passes((i,))]
    edges = {i: set() for i in vertices}
    for i, j in combinations(vertices, 2):
        if passes((i, j)):
            edges[i].add(j)

    # All remembered incompatibilities are exact failures of the permanent
    # identity at an integer linear combination, reduced modulo PRIME.
    largest = ()

    def extend(prefix, remaining):
        nonlocal largest
        if len(prefix) > len(largest):
            largest = prefix
        if len(prefix) == target:
            return True
        for offset, i in enumerate(remaining):
            if len(prefix) + len(remaining) - offset <= len(largest):
                break
            if all(passes((*group, i)) for count in range(2, min(4, len(prefix)) + 1)
                   for group in combinations(prefix, count)):
                tail = [j for j in remaining[offset+1:] if j in edges[i]]
                if extend((*prefix, i), tail):
                    return True
        return False

    extend((), vertices)
    return largest, len(vertices), passes.cache_info().misses


def main():
    parser = argparse.ArgumentParser(description='Search coordinate spans of cycle Fourier modes.')
    parser.add_argument('--target', type=int, default=9)
    args = parser.parse_args()
    if not 1 <= args.target <= 25:
        parser.error('target must be between 1 and 25')
    best = 0
    checked = 0
    for rows, columns in combinations_with_replacement(PARTITIONS, 2):
        largest, vertices, checks = search(rows, columns, args.target)
        best = max(best, len(largest))
        checked += checks
        print(f'{rows} x {columns}: {vertices} modes pass individually; '
              f'largest passing span {len(largest)}; {checks} tested supports.', flush=True)
        if len(largest) == args.target:
            print(f'Modular candidate (requires an independent complex-field check): {largest}')
            return
    print(f'Across all 28 partition pairs, maximum modular dimension {best}; {checked} tested supports.')


if __name__ == '__main__':
    main()
