import functools
import numpy as np


@functools.lru_cache(maxsize=None)
def _legendre_polynomial(k):
    """Returns k'th Legendre Polynomial using their recursive definition"""
    # assert k >= 0, f"Error! Called _legendre_polynomial with {k} < 0"
    if k == 0:
        return lambda x: 1
    elif k == 1:
        return lambda x: x
    return (lambda x: ((2 * k - 1) * x * _legendre_polynomial(k - 1)(x) -
                       (k - 1) * _legendre_polynomial(k - 2)(x)) / k)


@functools.lru_cache(maxsize=None)
def _legendre_polynomial_deriv1(k):
    """Returns the first derivative of k'th Legendre Polynomial"""
    if k == 0:
        return lambda x: 0
    elif k == 1:
        return lambda x: 1
    return (lambda x: ((2 * k - 1) * (x * _legendre_polynomial_deriv1(k - 1)
                                      (x) + _legendre_polynomial(k - 1)(x)) -
                       (k - 1) * _legendre_polynomial_deriv1(k - 2)(x)) / k)


@functools.lru_cache(maxsize=None)
def _legendre_roots(n, eps=1e-16):
    """Returns the n+1 roots of the n+1'th Legendre Polynomial as a list"""
    roots = []
    for k in range(n + 1):
        x0 = np.cos((4 * k + 3) * np.pi / (4 * n + 6))
        xk = _newton_iteration(_legendre_polynomial(n + 1),
                               _legendre_polynomial_deriv1(n + 1), x0, eps)
        roots.append(xk)

    return roots


def _newton_iteration(f, f_deriv, x0, eps):
    """Does Newton using f and its first derivative, starting at x0"""
    for _ in range(10000):
        # prevent lack of sufficient convergence
        xk = x0 - (f(x0) / f_deriv(x0))
        if abs(xk - x0) < eps:
            return xk
        x0 = xk
    return xk


def _int_weights(roots):
    """Calculate the Integration weights using Legendre Polynomial"""
    n = len(roots) - 1
    return [
        2 * (1 - xk**2) / ((n + 1)**2 * _legendre_polynomial(n)(xk)**2)
        for xk in roots
    ]


def gauss_quadratur(n):
    """Approximate the integral of f on the interval [-1,1] using
    Gauss-Quadrature and (n+1) 'Stützstellen'"""
    roots = _legendre_roots(n)
    weights = _int_weights(roots)

    def intern(f):
        return sum((f(roots[k]) * weights[k] for k in range(0, n + 1)))

    return intern
