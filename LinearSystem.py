from Matrix import Matrix


class LinearSystem:
    def __init__(self, a, b):
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

        augmented = augmented.upper_triangular(False)
        x = Matrix([[0 for t in range(1)] for k in range(self.a.n)],
                           self.a.n, self.a.m + self.b.m)
    def solve_gauss_jordan(self):
        pass

    def solve_cramer(self):
        pass

    def solve_inverse(self):
        pass


if __name__ == "__main__":
    a = Matrix([[1, -2, 3], [2, -3, -4], [2, -5, 1]], 3, 3)
    b = Matrix([[2], [-5], [-2]], 3, 1)
    ls = LinearSystem(a, b)
    ls.solve_gaussian()
