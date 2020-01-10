import numpy as np


def _legendre_polynomial(k):
    """Returns k'th Legendre Polynomial using their recursive definition"""
    if k == 0:
        return lambda x: 1
    elif k == 1:
        return lambda x: x

    def lk(x):
        x_2 = 1
        x_1 = x
        for l in range(2, k + 1):
            x_0 = ((2 * l - 1) * x * x_1 - (l - 1) * x_2) / l
            x_2 = x_1
            x_1 = x_0
        return x_0

    return lk


def _legendre_polynomial_deriv1(k):
    """Returns the first derivative of k'th Legendre Polynomial"""
    if k == 0:
        return lambda x: 0
    elif k == 1:
        return lambda x: 1

    def lk1(x):
        x1_2 = 0
        x1_1 = 1
        x_1 = x
        for l in range(2, k + 1):
            x1_0 = ((2 * l - 1) * (x * x1_1 + x_1) - (l - 1) * x1_2) / l
            x1_2 = x1_1
            x1_1 = x1_0
            x_1 = _legendre_polynomial(l)(x)
        return x1_0

    return lk1


def _legendre_roots(n):
    """Returns the n+1 roots of the n+1'th Legendre Polynomial as a list"""
    roots = []
    for k in range(n + 1):
        x0 = np.cos((4 * k + 3) * np.pi / (4 * n + 6))
        xk = _newton_iteration(
            _legendre_polynomial(n + 1), _legendre_polynomial_deriv1(n + 1), x0
        )
        roots.append(xk)

    return roots


def _newton_iteration(f, f_deriv, x0):
    """Does Newton using f and its first derivative, starting at x0"""
    for _ in range(10000):
        # prevent lack of sufficient convergence
        xk = x0 - (f(x0) / f_deriv(x0))
        if abs(xk - x0) < 1e-15:
            return xk
        x0 = xk
    return xk


def _int_weights(roots):
    """Calculate the Integration weights using Legendre Polynomial"""
    n = len(roots) - 1
    return [
        2 * (1 - xk ** 2) / ((n + 1) ** 2 * _legendre_polynomial(n)(xk) ** 2)
        for xk in roots
    ]


def gauss_quadratur(n):
    """Approximate the integral of f on the interval [-1,1] using
    Gauss-Quadrature and (n+1) 'Stützstellen'"""
    roots = _legendre_roots(n)
    weights = _int_weights(roots)

    def cached(f):
        total = 0
        for (root, weight) in zip(roots, weights):
            total += f(root) * weight
        return total
        # return sum((f(root) * weight for (root, weight) in zip(roots, weights)))

    return cached
