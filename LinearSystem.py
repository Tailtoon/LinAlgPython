from Matrix import Matrix


class LinearSystem:
    def __init__(self, a: Matrix, b: Matrix):
        self.a = a
        self.b = b

    def solve_gaussian(self):
        augmented = Matrix([[0 for t in range(self.a.m + self.b.m)] for k in range(self.a.n)],
                           self.a.n, self.a.m + self.b.m)
        for i in range(augmented.n):
            for j in range(augmented.m):
                if j > augmented.m - 2:
                    augmented[i][j] = self.b[i][0]
                else:
                    augmented[i][j] = self.a[i][j]
        rank_before = augmented.rank()
        augmented = augmented.upper_triangular(clear_zero_rows=True)
        if rank_before != augmented.n:
            raise Exception("Solve error: Kronecker–Capelli (Rouché–Capelli) rank(self) != rank(upper_triangular)")
        if rank_before == augmented.m - 1:  # rank(self) = n - only one solution
            x = Matrix([[0] for k in range(augmented.m)], augmented.m - 1, 1)
            for i in range(x.n):
                curi = augmented.n - i - 1
                curj = augmented.m - 1
                s = 0
                for k in range(i):
                    s += augmented[curi][curj - k - 1] * x[curj - k - 1][0]
                x[curi][0] = augmented[curi][curj] - s
            return x
        else:  # more than one solution
            pass

    def solve_gauss_jordan(self):
        augmented = Matrix([[0 for t in range(self.a.m + self.b.m)] for k in range(self.a.n)],
                           self.a.n, self.a.m + self.b.m)
        for i in range(augmented.n):
            for j in range(augmented.m):
                if j > augmented.m - 2:
                    augmented[i][j] = self.b[i][0]
                else:
                    augmented[i][j] = self.a[i][j]
        rank_before = augmented.rank()
        augmented = augmented.gauss_jordan()
        if rank_before != augmented.n:
            raise Exception("Solve error: Kronecker–Capelli (Rouché–Capelli) rank(self) != rank(upper_triangular)")
        if rank_before == augmented.m - 1:  # rank(self) = n - only one solution
            x = Matrix([[0] for k in range(augmented.m)], augmented.m - 1, 1)
            for i in range(x.n):
                curi = augmented.n - i - 1
                curj = augmented.m - 1
                x[curi][0] = augmented[curi][curj]
            return x
        else:  # more than one solution
            pass

    def solve_cramer(self):
        det_a = self.a.det()
        x = Matrix([[0] for k in range(self.a.m)], self.a.m, 1)
        for i in range(x.n):
            replaced_matrix = self.a.copy()
            replaced_matrix[i] = self.b
            pass

    def solve_inverse(self):
        inverted = self.a.inverse_adjugate()
        x = inverted * self.b
        return x


if __name__ == "__main__":
    a = Matrix([[1, -2, 3], [2, -3, -4], [2, -5, 1]], 3, 3)
    b = Matrix([[2], [-5], [-2]], 3, 1)
    ls = LinearSystem(a, b)
    print(ls.solve_gaussian())
    print(ls.solve_gauss_jordan())
    print(ls.solve_inverse())
