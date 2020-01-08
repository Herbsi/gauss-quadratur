class Poly:
    def __init__(self, coeffs):
        self.coefficients = coeffs
        self.degree = len(coeffs) - 1

    def nth_coefficient(self, n):
        try:
            return self.coefficients[n]
        except:
            return 0

    def __str__(self):
        pprint = f"{self.coefficients[0]} "
        for (k, ak) in zip(range(1, 1 + self.degree), self.coefficients[1:]):
            if ak != 0:
                if ak != 1:
                    pprint += f"+ {ak}*x^{k} "
                else:
                    pprint += f"+ x^{k} "
        return pprint

    def __add__(self, other):
        return Poly(
            [
                self.nth_coefficient(n) + other.nth_coefficient(n)
                for n in range(1 + max(self.degree, other.degree))
            ]
        )

    def __sub__(self, other):
        return Poly(
            [
                self.nth_coefficient(n) - other.nth_coefficient(n)
                for n in range(1 + max(self.degree, other.degree))
            ]
        )

    def __mul__(self, other):
        def kth_coeff(k):
            return sum(
                [
                    self.nth_coefficient(l) * other.nth_coefficient(k - l)
                    for l in range(k + 1)
                ]
            )

        return Poly([kth_coeff(k) for k in range(1 + self.degree + other.degree)])

    def derivative(self):
        return Poly(
            [a * k for (k, a) in zip(range(1, self.degree + 1), self.coefficients[1:])]
        )

    def eval(self, x):
        total = 0
        for ak in self.coefficients[::-1]:
            total = ak + x * total
        return total
