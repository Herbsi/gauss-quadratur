#!/usr/bin/env python3
import numpy as np
import legendre


def main():
    f = lambda x: np.log(x + 2)
    actual_value = -2 + np.log(27)
    integral_values = []
    errors = []
    for n in [2, 4, 8, 16]:
        int_value = legendre.gauss_quadratur(f, n)
        error = abs(int_value - actual_value)
        integral_values.append(int_value)
        errors.append(error)
    for (i, n) in enumerate([2, 4, 8, 16]):
        # TODO nicer formatting
        print(f"n = {n}, int = {integral_values[i]}, err = {errors[i]}\n")


if __name__ == "__main__":
    main()
