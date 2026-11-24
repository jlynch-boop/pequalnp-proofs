from fractions import Fraction
from itertools import permutations, product
from math import prod


SIZE = 5
CORE = SIZE - 1
COORDINATES = range(SIZE)
MINOR_PERMUTATIONS = tuple(permutations(range(CORE)))
MATRIX_PERMUTATIONS = tuple(permutations(COORDINATES))


def star_direction(row, column):
    matrix = [[0] * SIZE for _ in COORDINATES]
    matrix[-1][:-1] = [*row, -sum(row)]
    for index, entry in enumerate(column):
        matrix[index][-1] = entry
    matrix[-1][-1] = -sum(column)
    return matrix


def gradient(matrix):
    result = []
    for row, column in product(COORDINATES, repeat=2):
        minor = [[value for j, value in enumerate(line) if j != column]
                 for i, line in enumerate(matrix) if i != row]
        result.append(sum(minor[0][p[0]] * minor[1][p[1]]
                          * minor[2][p[2]] * minor[3][p[3]]
                          for p in MINOR_PERMUTATIONS))
    return result


def permanent_along(direction, scalar):
    return sum(prod(1 + scalar * direction[i][p[i]] for i in COORDINATES)
               for p in MATRIX_PERMUTATIONS)


def exact_rank(rows):
    basis = {}
    for row in rows:
        reduced = list(map(Fraction, row))
        for index, pivot in sorted(basis.items()):
            factor = reduced[index]
            reduced = [entry - factor * component for entry, component in zip(reduced, pivot)]
        for index, entry in enumerate(reduced):
            if entry:
                basis[index] = [value / entry for value in reduced]
                break
    return len(basis)


def standard_samples(dimension):
    return [tuple(0 for _ in range(dimension)), *[
        tuple(int(i == j) for i in range(dimension)) for j in range(dimension)
    ]]


def flatten(matrix):
    return [value for row in matrix for value in row]


def main():
    rows = standard_samples(CORE - 1)
    columns = standard_samples(CORE)
    gradients = [gradient([[1 + value for value in line]
                           for line in star_direction(row, column)])
                 for row, column in product(rows, columns)]
    directions = [star_direction(row, columns[0]) for row in rows[1:]]
    directions.extend(star_direction(rows[0], column) for column in columns[1:])
    extra = [[1] * CORE + [int(i == CORE) * -20] for i in COORDINATES]
    directions.append(extra)
    tangent_basis = [flatten(direction) for direction in directions]
    assert exact_rank(gradients) == 17
    assert exact_rank(tangent_basis) == 8
    assert all(sum(x * y for x, y in zip(row, vector)) == 0
               for row, vector in product(gradients, tangent_basis))
    # Both sides have degree at most five, so six distinct values certify identity.
    assert all(permanent_along(extra, s) == 120 * (1+s)**4 * (1-4*s) for s in range(6))
    print('Gradient rank 17. The 8-dimensional tangent kernel has the stated basis.')
    print('The extra direction has permanent polynomial 120(1+z)^4(1-4z).')


if __name__ == '__main__':
    main()
