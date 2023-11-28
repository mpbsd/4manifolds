#!/usr/bin/env python3
# -*- coding: utf-8 -*- #


import sympy as sp


# pseudo_riemannian_metric class {{{
class pseudo_riemannian_metric:
    def __init__(self, diag, coord):
        self.x = coord
        self.g = sp.diag(diag[0], diag[1], diag[2], diag[3])
        self.g_inv = self.g**-1

    def christoffel_symbol(self, i, j, k):
        return sum(
            [
                0.5
                * (
                    + self.g[i, l].diff(self.x[j], 1)
                    + self.g[j, l].diff(self.x[i], 1)
                    - self.g[i, j].diff(self.x[k], 1)
                )
                * self.g_inv[l, k]
                for l in range(4)
            ]
        )

    # R_{ijk}^{l}
    def curvature_tensor(self, i, j, k, l):
        return (
            + self.christoffel_symbol(i, k, l).diff(self.x[j], 1)
            - self.christoffel_symbol(i, j, l).diff(self.x[k], 1)
        ) + sum(
            [
                + self.christoffel_symbol(i, k, m)
                * self.christoffel_symbol(m, j, l)
                - self.christoffel_symbol(i, j, m)
                * self.christoffel_symbol(m, k, l)
                for m in range(4)
            ]
        )

    # R_{ijkl}
    def curvature_tensor_lower_ind(self, i, j, k, l):
        return sum(
            [
                self.curvature_tensor(i, j, k, m) * self.g[m, l]
                for m in range(4)
            ]
        )

    def sectional_curvature(self, i, j):
        return (
            self.curvature_tensor_lower_ind(i, j, i, j)
            / (self.g[i, i] * self.g[j, j])
            if i != j
            else None
        )

    def ricci_tensor(self, i, j):
        return sum([self.curvature_tensor(i, k, j, k) for k in range(4)])

    def scalar_curvature(self):
        return sum(
            [
                sum(
                    [
                        self.ricci_tensor(i, j) * self.g_inv[j, i]
                        for j in range(4)
                    ]
                )
                for i in range(4)
            ]
        )


# }}}


def main():
    x = sp.Array([sp.symbols("x_%d" % i) for i in range(4)])

    # sigma = [1, 2, 3, 0]

    # u_0 = sp.Function("u_0")(x[sigma[0]])
    # u_1 = sp.Function("u_1")(x[sigma[1]])
    # u_2 = sp.Function("u_2")(x[sigma[2]])
    # u_3 = sp.Function("u_3")(x[sigma[3]])

    u_0 = 0
    u_1 = sp.ln(sp.cosh(x[2]))
    u_2 = 0
    u_3 = sp.ln(sp.cosh(x[0]))

    diag = [
        +sp.exp(2 * u_0),
        +sp.exp(2 * u_1),
        +sp.exp(2 * u_2),
        +sp.exp(2 * u_3),
    ]

    diagonal_metric = pseudo_riemannian_metric(diag, x)

    with open("ricci_tensor.tex", "w") as log_file:
        for i in range(4):
            for j in range(i, 4):
                print(
                    "%d, %d: %s"
                    % (i, j, sp.latex(diagonal_metric.ricci_tensor(i, j))),
                    file=log_file,
                )
        print("---", file=log_file)
        for i in range(3):
            for j in range(i + 1, 4):
                print(
                    "%d, %d: %s"
                    % (i, j, sp.latex(diagonal_metric.sectional_curvature(i, j))),
                    file=log_file,
                )


if __name__ == "__main__":
    main()
