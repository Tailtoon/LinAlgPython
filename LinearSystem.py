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
        rank_before = self.a.rank()
        augmented = augmented.upper_triangular(clear_zero_rows=True)
        if rank_before != augmented.n:
            raise Exception("Solve error: Kronecker–Capelli (Rouché–Capelli) rank(self) != rank(upper_triangular)")
        if rank_before == self.a.m:  # rank(self) = a.m - only one solution
            x = Matrix([[0] for k in range(augmented.m - 1)], augmented.m - 1, 1)
            for i in range(x.n):
                curi = augmented.n - i - 1
                curj = augmented.m - 1
                s = 0
                for k in range(i):
                    s += augmented[curi][curj - k - 1] * x[curj - k - 1][0]
                x[curi][0] = augmented[curi][curj] - s
            return x
        else:  # more than one solution
            # It will be easier to get a general solution by doing gauss-jordan solving method
            for curij in range(augmented.n - 1, -1, -1):
                for i in range(1, curij + 1):  # Get zeros above 1 (skipped if k == range minimum)
                    muler = augmented[curij - i][curij]
                    for j in range(augmented.m):
                        augmented[curij - i][j] -= muler * augmented[curij][j]
            x = Matrix([[""] for k in range(self.a.m)], self.a.m, 1)
            c_name_count = 1
            for i in range(x.n - 1, -1, -1):
                if i > augmented.n - 1:
                    x[i][0] += "C{}".format(c_name_count)
                    c_name_count += 1
                else:
                    for k in range(self.a.m - rank_before):
                        x[i][0] += "({} ".format(-augmented[i][augmented.n + k]) + x[augmented.n + k][0] + ") + "
                    x[i][0] += "(" + str(augmented[i][augmented.m - 1]) + ")"
            return x

    def solve_gauss_jordan(self):
        augmented = Matrix([[0 for t in range(self.a.m + self.b.m)] for k in range(self.a.n)],
                           self.a.n, self.a.m + self.b.m)
        for i in range(augmented.n):
            for j in range(augmented.m):
                if j > augmented.m - 2:
                    augmented[i][j] = self.b[i][0]
                else:
                    augmented[i][j] = self.a[i][j]
        rank_before = self.a.rank()
        augmented = augmented.gauss_jordan()
        if rank_before != augmented.n:
            raise Exception("Solve error: Kronecker–Capelli (Rouché–Capelli) rank(self) != rank(upper_triangular)")
        if rank_before == self.a.m:  # rank(self) = m - only one solution
            x = Matrix([[0] for k in range(augmented.m - 1)], augmented.m - 1, 1)
            for i in range(x.n):
                curi = augmented.n - i - 1
                curj = augmented.m - 1
                x[curi][0] = augmented[curi][curj]
            return x
        else:  # more than one solution
            x = Matrix([[""] for k in range(self.a.m)], self.a.m, 1)
            c_name_count = 1
            for i in range(x.n - 1, -1, -1):
                if i > augmented.n - 1:
                    x[i][0] += "C{}".format(c_name_count)
                    c_name_count += 1
                else:
                    for k in range(self.a.m - rank_before):
                        x[i][0] += "({} ".format(-augmented[i][augmented.n + k]) + x[augmented.n + k][0] + ") + "
                    x[i][0] += "(" + str(augmented[i][augmented.m - 1]) + ")"
            return x

    def solve_cramer(self):
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
        if rank_before != self.a.m:
            raise Exception("Solve error: Cramer's rule works only for systems with unique solution")
        det_a = self.a.det()
        x = Matrix([[0] for k in range(self.a.m)], self.a.m, 1)
        for i in range(x.n):
            replaced_matrix = self.a.copy()
            replaced_matrix.replace_col(i, self.b)
            x[i][0] = replaced_matrix.det() / det_a
        return x

    def solve_inverse(self):
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
        if rank_before != self.a.m:
            raise Exception("Solve error: Inverse solve method works only for systems with unique solution")
        inverted = self.a.inverse_adjugate()
        x = inverted * self.b
        return x

    def solve_LU(self):
        LU = self.a.LU_decomposition(compact_return=True)
        y = Matrix([[0] for k in range(self.a.m)], self.a.m, 1)
        for i in range(y.n):
            s = 0
            for k in range(i):
                s += LU[i][k] * y[k][0]
            y[i][0] = self.b[i][0] - s
        x = Matrix([[0] for k in range(self.a.m)], self.a.m, 1)
        for i in range(x.n - 1, -1, -1):
            s = 0
            for k in range(i, x.n):
                s += LU[i][k] * x[k][0]
            x[i][0] = (y[i][0] - s) / LU[i][i]  # Zero division error excluded by LU_decomposition
        return x


if __name__ == "__main__":
    a = Matrix([[1, -2, 3], [2, -3, -4], [2, -5, 1]], 3, 3)
    b = Matrix([[2], [-5], [-2]], 3, 1)
    print(a)
    print(b)
    print("det(A) = ", a.det())
    print("Ранг A = ", a.rank())
    L, U = a.LU_decomposition()
    print(L)
    print(U)
    print(a.LU_decomposition(compact_return=True))
    print(L * U)
    augmented = Matrix([[0 for t in range(a.m + b.m)] for k in range(a.n)],
                       a.n, a.m + b.m)
    for i in range(augmented.n):
        for j in range(augmented.m):
            if j > augmented.m - 2:
                augmented[i][j] = b[i][0]
            else:
                augmented[i][j] = a[i][j]
    augmented = augmented.gauss_jordan()
    print("Ранг A|B = ", augmented.n)
    print(a.inverse_adjugate())
    ls = LinearSystem(a, b)
    print(ls.solve_gaussian())
    print(ls.solve_gauss_jordan())
    print(ls.solve_inverse())
    print(ls.solve_cramer())
    print(ls.solve_LU())
    a = Matrix([[1, 2, -1, 2, 2], [1, 0, -3, 6, 1], [0, 4, 2, -1, 1]], 3, 5)
    print(a)
    b = Matrix([[2], [4], [-2]], 3, 1)
    print(b)
    print("Ранг A = ", a.rank())
    augmented = Matrix([[0 for t in range(a.m + b.m)] for k in range(a.n)],
                       a.n, a.m + b.m)
    for i in range(augmented.n):
        for j in range(augmented.m):
            if j > augmented.m - 2:
                augmented[i][j] = b[i][0]
            else:
                augmented[i][j] = a[i][j]
    augmented = augmented.gauss_jordan()
    print("Ранг A|B = ", augmented.n)
    ls = LinearSystem(a, b)
    print(ls.solve_gaussian())
    print(a * Matrix([[1], [0], [-1], [0], [0]], 5, 1))
