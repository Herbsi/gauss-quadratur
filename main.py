#!/usr/bin/env python3
import math
import numpy as np
import legendre
from polynomial import Poly


def main():
    f = lambda x: math.log(x + 2)
    actual_value = math.log(27) - 2
    integral_values = []
    errors = []
    for n in [2, 4, 8, 16]:
        int_value = legendre.gauss_quadratur(f, n)
        error = abs(int_value - actual_value)
        print(f"n = {n}, int = {int_value}, err={error}")


if __name__ == "__main__":
    main()
