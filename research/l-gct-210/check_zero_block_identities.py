from fractions import Fraction
from itertools import combinations, permutations, product
from math import prod


PERMUTATIONS_3 = tuple(permutations(range(3)))


def per3(matrix):
    return sum(prod(matrix[row][p[row]] for row in range(3)) for p in PERMUTATIONS_3)


def check_projection_charts():
    # Every formula is degree at most three in each parameter separately.
    for r, s, z0, z1, z2, z3, z4 in product(range(4), repeat=7):
        product_chart = [[s*z0 + r*z2 - r*s*z4, -r*z3, z0],
                         [-s*z1, 0, z1], [z2, z3, z4]]
        hirzebruch_chart = [[-r*z2, r*z3 - s*z0, r*z4 + s*z1],
                            [0, z0, z1], [z2, z3, z4]]
        assert per3(product_chart) == per3(hirzebruch_chart) == 0
        assert hirzebruch_chart[0][1] + s*z0 - r*z3 == 0
        assert hirzebruch_chart[0][2] - s*z1 - r*z4 == 0

    # Substitute the printed plus sign instead, with r=1 and s=0.
    assert per3([[-1, -1, 0], [0, 0, 1], [1, 1, 0]]) == -2


def cofactors(direction):
    left = [[1, 1, 1], [1, 1, 1]] + [[1 + x for x in row] for row in direction]
    result = [[0] * 5 for _ in range(5)]
    for i, j in combinations(range(5), 2):
        value = per3([row for r, row in enumerate(left) if r not in (i, j)])
        result[i][j] = result[j][i] = value
    return result


def two_subpermanent_sum(x, y):
    return sum(x) * sum(y) - sum(a*b for a, b in zip(x, y))


def check_radical_coefficient():
    # Bottom row h=3 in the 5x3 left block, with p=2 and q=4.
    coordinate_rows = [[0] * 8] + [[int(i == j) for i in range(8)] for j in range(8)]
    probes = coordinate_rows + [[a+b for a, b in zip(x, y)]
                                for x, y in combinations(coordinate_rows[1:], 2)]
    base = cofactors([[0, 0, 0] for _ in range(3)])
    for d in probes:
        direction = [d[:3], [d[3], d[4], -d[3]-d[4]], d[5:]]
        positive = cofactors(direction)
        negative = cofactors([[-x for x in row] for row in direction])
        quadratic = [[Fraction(a+b, 2)-c for a, b, c in zip(pos, neg, zero)]
                     for pos, neg, zero in zip(positive, negative, base)]
        q_p = two_subpermanent_sum(direction[1], direction[2])
        q_q = two_subpermanent_sum(direction[0], direction[1])
        for A, B in [(1, 0), (0, 1)]:
            level = -Fraction(A+B, 3)
            x = [level, level, A, level, B]
            for y in [(1, -1, 0, 0, 0), (1, 0, 0, -1, 0)]:
                actual = sum(x[i]*quadratic[i][j]*y[j] for i, j in product(range(5), repeat=2))
                expected = Fraction(y[3]*((B-2*A)*q_p + (A-2*B)*q_q), 3)
                assert actual == expected


def check_exceptional_cofactors():
    # Every entry is degree at most three in each parameter separately.
    group = {0, 1, 3}
    for a, b, c, d, e in product(range(4), repeat=5):
        matrix = cofactors([[a, a, b], [c, -c, 0], [d, d, e]])
        u, v, s = 2*a+b, 2*d+e, 2*(a*d+a*e+b*d)
        for i, j in combinations(range(5), 2):
            if i in group and j in group:
                expected = 6+2*(u+v)+s
            elif i == 2 and j == 4:
                expected = 6
            elif i == 2 or j == 2:
                expected = 6+2*v
            else:
                expected = 6+2*u
            assert matrix[i][j] == expected
    assert 2*(1*1 + 1*0 + 0*1) != 2*(2*1+0)*(2*1+0)


if __name__ == '__main__':
    check_projection_charts()
    check_radical_coefficient()
    check_exceptional_cofactors()
    print('Exact checks passed for the projection charts, radical coefficient, and exceptional cofactors.')
