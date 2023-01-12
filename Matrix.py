import random


class Matrix:
    PRECISION = 1e-10

    def __init__(self, mat, n, m):
        self.n = n
        self.m = m
        self.matrix = mat

    def __add__(self, other):  # self + other
        if self.n != other.n or self.m != other.m:
            raise Exception("Add error: Matrix sizes aren't equal")
        res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = self[i][j] + other[i][j]
        return res

    def __sub__(self, other):  # self - other
        if self.n != other.n or self.m != other.m:
            raise Exception("Sub error: Matrix sizes aren't equal")
        res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = self[i][j] - other[i][j]
        return res

    def __mul__(self, other):  # self x other
        if isinstance(other, float) or isinstance(other, int):
            res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
            for i in range(self.n):
                for j in range(self.m):
                    res[i][j] = self[i][j] * other
            return res
        else:
            return self.dot(other)

    def __truediv__(self, other):  # self / other (only for numbers)
        if not isinstance(other, float) and not isinstance(other, int):
            raise Exception("Div error: " + type(other) + " isn't a number")
        res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = self[i][j] / other
        return res

    def __eq__(self, other):  # self == other
        if (self.n != other.n) or (self.m != other.m):
            return False
        for i in range(self.n):
            for j in range(self.m):
                if self[i][j] != other[i][j]:
                    return False
        return True

    def __ne__(self, other):  # self != other
        return not self.__eq__(other)

    def __getitem__(self, item):
        return self.matrix[item]

    def __str__(self):
        s = ""
        for i in range(self.n):
            s += self[i].__str__() + "\n"
        return s

    def dot(self, other):  # self x other_matrix
        if self.m != other.n:
            raise Exception("Dot error: Number of columns isn't equal to number of rows")
        res = Matrix([[0 for j in range(other.m)] for i in range(self.n)], self.n, other.m)
        for i in range(self.n):
            for j in range(other.m):
                res[i][j] = 0
                for k in range(self.m):
                    res[i][j] += self[i][k] * other[k][j]
        return res

    def transpose(self):  # creates new matrix without overwriting an existing one
        res = Matrix([[0 for j in range(self.n)] for i in range(self.m)], self.m, self.n)
        for i in range(self.m):
            for j in range(self.n):
                res[i][j] = self[j][i]
        return res

    def det(self):  # Calculate determinant of Matrix
        if self.n != self.m:
            raise Exception("Det error: not square matrix")
        return Matrix._det(self.matrix)

    @staticmethod
    def _det(matrix):  # Calculate determinant of matrix
        n = len(matrix)
        if n == 1:
            return matrix[0][0]
        s = 0
        for j in range(n):
            s += matrix[0][j] * Matrix._minor(matrix, 0, j) * ((-1) ** (2 + j))
        return s

    @staticmethod
    def _minor(matrix, i, j):  # Calculate minor of matrix without i row and j column
        minor = Matrix._del_iRow_jCol(matrix, i, j)
        return Matrix._det(minor)

    @staticmethod
    def _del_iRow_jCol(matrix, i, j):
        n = len(matrix)
        minor = [[0 for t in range(n - 1)] for k in range(n - 1)]
        skipi = 0
        for k in range(n - 1):
            if k == i:
                skipi = 1
            skipj = 0
            for t in range(n - 1):
                if t == j:
                    skipj = 1
                minor[k][t] = matrix[k + skipi][t + skipj]
        return minor

    def is_invertible(self):
        return False if abs(self.det()) < self.PRECISION or (self.n != self.m) else True

    def inverse_adjugate(self):
        if not self.is_invertible():
            raise Exception("Inverse error: Matrix isn't invertible")
        transposed = self.transpose()
        cofactor_matrix = Matrix([[0 for t in range(self.n)] for k in range(self.n)], self.n, self.n)
        for i in range(self.n):
            for j in range(self.n):
                cofactor_matrix[i][j] = Matrix._minor(transposed.matrix, i, j) * ((-1) ** (2 + i + j))
        res = cofactor_matrix * (1 / self.det())
        return res

    def upper_triangular(self, clear_zero_rows=True):
        uppmatrix = Matrix([[0 for t in range(self.m)] for k in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                uppmatrix[i][j] = self[i][j]
        for curij in range(uppmatrix.n):
            if abs(uppmatrix[curij][curij]) < self.PRECISION:  # Get nonzero element at top
                for i in range(uppmatrix.n - curij):
                    if abs(uppmatrix[curij + i][curij]) > self.PRECISION:
                        uppmatrix.matrix[curij], uppmatrix.matrix[curij + i] = \
                            uppmatrix.matrix[curij + i], uppmatrix.matrix[curij]
                        break
            if abs(uppmatrix[curij][curij]) < self.PRECISION:  # If we didn't find element nonzero - continue
                continue
            diver = uppmatrix[curij][curij]
            for j in range(self.m):  # Get 1 element at top
                uppmatrix[curij][j] /= diver

            for i in range(1, uppmatrix.n - curij):  # Get zeros below 1 (skipped if k == range maximum)
                muler = uppmatrix[curij + i][curij]
                for j in range(uppmatrix.m):
                    uppmatrix[curij + i][j] -= muler * uppmatrix[curij][j]
        if clear_zero_rows:
            delete_list = [False for t in range(uppmatrix.n)]
            for i in range(uppmatrix.n):
                delete = True
                for j in range(uppmatrix.m):
                    if uppmatrix[i][j] > self.PRECISION:
                        delete = False
                        break
                delete_list[i] = delete

            deleted_cnt = 0
            for i in range(len(delete_list)):
                if delete_list[i]:
                    uppmatrix.matrix.pop(i - deleted_cnt)
                    uppmatrix.n -= 1
                    deleted_cnt += 1
        return uppmatrix

    def gauss_jordan(self):
        res = self.upper_triangular()
        for curij in range(res.n - 1, -1, -1):
            for i in range(1, curij + 1):  # Get zeros above 1 (skipped if k == range minimum)
                muler = res[curij - i][curij]
                for j in range(res.m):
                    res[curij - i][j] -= muler * res[curij][j]
        return res

    def inverse_gaussian(self):
        if not self.is_invertible():
            raise Exception("Inverse error: Matrix isn't invertible")
        attached_matrix = Matrix([[0 for t in range(self.n * 2)] for k in range(self.n)], self.n, self.n * 2)
        for i in range(self.n):
            for j in range(self.n):
                attached_matrix[i][j] = self[i][j]
        for i in range(self.n):
            attached_matrix[i][i + self.m] = 1

        attached_matrix = attached_matrix.gauss_jordan()

        res = Matrix([[0 for t in range(self.n)] for k in range(self.n)], self.n, self.n)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = attached_matrix[i][j + self.n]
        return res

    def rank(self):
        row_echelon = self.upper_triangular()
        return row_echelon.n


if __name__ == '__main__':
    a = Matrix([[random.randint(0, 5) for j in range(4)] for i in range(3)], 3, 4)
    print("Матрица A")
    print(a)
    b = a.transpose()
    print("Матрица B")
    print(b)
    c = a * b
    print("Матрица C = A x B")
    print(c)
    print("Определитель матрицы C = ", c.det())
    print(c.is_invertible())
    print(c.inverse_adjugate())
    print(c.inverse_gaussian())
    print(c.rank())
    # d = Matrix([[2, 5, 7], [6, 3, 4], [5, -2, -3]], 3, 3)
    # print(d)
    # print(d.inverse_adjugate())
    # gg = Matrix([[3, 1, -1, -2, 8], [7, 1, -2, -1, 12], [11, 1, -3, 0, 16], [2, 2, -1, -5, 12]], 4, 5)
    # print(gg)
    # gg1 = gg.upper_triangular()
    # print(gg1)
    # gg2 = gg.upper_triangular(False)
    # print(gg2)
