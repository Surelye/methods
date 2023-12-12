import random as rnd
from math import isqrt
from egcd import egcd


TRIES = 10000


def is_B_smooth(b, factor_base):
    len_factor_base = len(factor_base)
    canonical_decomposition = [0] * len_factor_base
    for i in range(len_factor_base):
        while b % factor_base[i] == 0:
            b /= factor_base[i]
            canonical_decomposition[i] += 1
    return canonical_decomposition if b == 1 else False


def read_params():
    with open("params", 'r') as f:
        return list(map(lambda num: int(num), f.readline().split()))


def multiply_matrices(a, b, p):
    result, b_n, b_m = list(), len(b), len(b[0])
    for i in range(len(a)):
        vec = list()
        for j in range(b_m):
            elt = 0
            for k in range(b_n):
                elt = (elt + a[i][k] * b[k][j]) % p
            vec.append(elt)
        result.append(vec)
    return result


def modified_Jordan(equations, p):
    n = len(equations)
    for i in range(n):
        for j in range(i + 1, n):
            rs = [equations[j][i], (-1 * equations[i][i]) % p]
            egcd_res = egcd(equations[i][i], equations[j][i])
            xy = [egcd_res[1] % p, egcd_res[2] % p]
            a = [xy, rs]
            b = [equations[i], equations[j]]
            res = multiply_matrices(a, b, p)
            equations[i] = res[0]
            equations[j] = res[1]
        a_ii_egcd = egcd(equations[i][i], p)
        if a_ii_egcd[0] != 1:
            return list()
        a_ii_inv = a_ii_egcd[1] % p
        for k in range(len(equations[i])):
            equations[i][k] = (equations[i][k] * a_ii_inv) % p
        for j in range(i):
            a_ji = equations[j][i]
            for k in range(len(equations[j])):
                equations[j][k] = (equations[j][k] - equations[i][k] * a_ji) % p
    return [equations[i][-1] for i in range(n)]


def isprime(p):
    for i in range(3, isqrt(p) + 1):
        if p % i == 0:
            return False
    return True


def index_method():
    p, a, h, B = read_params()
    factor_base = [2] + [prime for prime in range(3, B + 1, 2) if isprime(prime)]
    equations, x = list(), list()
    len_factor_base, len_equations = len(factor_base), len(equations)

    attempt = 0
    while len_equations != len_factor_base and attempt < TRIES:
        m = rnd.randint(0, p - 2)
        b = pow(a, m, p)
        canonical_decomposition = is_B_smooth(b, factor_base)
        if canonical_decomposition:
            canonical_decomposition.append(m)
            equations.append(canonical_decomposition)
            len_equations += 1
        else:
            continue
        x = modified_Jordan(equations, p - 1)
        if len(x) == 0:
            equations = list()
            len_equations = 0
            attempt += 1
            continue

    if len(x) == 0:
        return list()
    while True:
        m = rnd.randint(0, p - 2)
        b = (pow(a, m, p) * h) % p
        ri = is_B_smooth(b, factor_base)
        if not ri:
            continue
        res = 0
        for i in range(len(ri)):
            res = res + (ri[i] * x[i]) % (p - 1)
        res = (res - m) % (p - 1)
        return res


def main():
    res = index_method()
    f = open("result", "w")
    f.write(str(res))
    f.close()


if __name__ == '__main__':
    main()
