class Matrix:
    PRECISION = 1e-10

    def __init__(self, mat: list, n: int, m: int):
        self.n = n
        self.m = m
        self.matrix = mat

    def __add__(self, other):  # self + other
        if self.n != other.n or self.m != other.m:
            raise Exception("Add error: Matrix sizes aren't equal")
        res = Matrix([[self[i][j] + other[i][j] for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        return res

    def __sub__(self, other):  # self - other
        if self.n != other.n or self.m != other.m:
            raise Exception("Sub error: Matrix sizes aren't equal")
        res = Matrix([[self[i][j] - other[i][j] for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        return res

    def __mul__(self, other):  # self x other
        if isinstance(other, float) or isinstance(other, int):
            res = Matrix([[self[i][j] * other for j in range(self.m)] for i in range(self.n)], self.n, self.m)
            return res
        else:
            return self.dot(other)

    def __truediv__(self, other):  # self / other (only for numbers)
        if not isinstance(other, float) and not isinstance(other, int):
            raise Exception("Div error: " + type(other) + " isn't a number")
        res = Matrix([[self[i][j] / other for j in range(self.m)] for i in range(self.n)], self.n, self.m)
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

    def copy(self):
        return Matrix([[self[i][j] for j in range(self.m)] for i in range(self.n)], self.n, self.m)

    def replace_row(self, i: int, value):
        if isinstance(value, list):
            if len(list) != self.m:
                raise Exception("Replace_row error: invalid len(list) = {}".format(len(value)))
            self.matrix[i] = value
        elif isinstance(value, Matrix):
            if value.n != 1 or value.m != self.m:
                raise Exception("Replace_row error: invalid value Matrix size = ({}, {})".format(value.n, value.m))
            for j in range(self.m):
                self[i][j] = value[0][j]
        else:
            raise Exception("Replace_row error: invalid value class " + type(value))

    def replace_col(self, j: int, value):
        if isinstance(value, list):
            if len(list) != self.n:
                raise Exception("Replace_col error: invalid len(list) = {}".format(len(value)))
            for i in range(self.n):
                self[i][j] = value[i]
        elif isinstance(value, Matrix):
            if value.m != 1 or value.n != self.n:
                raise Exception("Replace_col error: invalid value Matrix size = ({}, {})".format(value.n, value.m))
            for i in range(self.n):
                self[i][j] = value[i][0]
        else:
            raise Exception("Replace_col error: invalid value class " + type(value))

    def clear_zero_rows(self):
        delete_list = [False for t in range(self.n)]
        for i in range(self.n):
            delete = True
            for j in range(self.m):
                if self[i][j] > self.PRECISION:
                    delete = False
                    break
            delete_list[i] = delete

        deleted_cnt = 0
        for i in range(len(delete_list)):
            if delete_list[i]:
                self.matrix.pop(i - deleted_cnt)
                self.n -= 1
                deleted_cnt += 1
        return self

    @staticmethod
    def zeros(n: int, m: int):
        return Matrix([[0 for j in range(m)] for i in range(n)], n, m)

    def dot(self, other):  # self x other_matrix
        if self.m != other.n:
            raise Exception("Dot error: Number of columns isn't equal to number of rows")
        res = Matrix.zeros(self.n, other.m)
        for i in range(self.n):
            for j in range(other.m):
                for k in range(self.m):
                    res[i][j] += self[i][k] * other[k][j]
        return res

    def transpose(self):  # creates new matrix without overwriting an existing one
        res = Matrix.zeros(self.m, self.n)
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
        if n == 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        if n == 3:
            return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] \
                + matrix[1][0] * matrix[2][1] * matrix[0][2] - matrix[0][2] * matrix[1][1] * matrix[2][0] \
                - matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][1] * matrix[0][0]
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
        return False if (self.n != self.m) or abs(self.det()) < self.PRECISION else True

    def inverse_adjugate(self, check_is_invertible=True):
        if check_is_invertible:
            if not self.is_invertible():
                raise Exception("Inverse error: Matrix isn't invertible")
        transposed = self.transpose()
        cofactor_matrix = Matrix.zeros(self.n, self.n)
        for i in range(self.n):
            for j in range(self.n):
                cofactor_matrix[i][j] = Matrix._minor(transposed.matrix, i, j) * ((-1) ** (2 + i + j))
        res = cofactor_matrix * (1 / self.det())
        return res

    def upper_triangular(self, clear_zero_rows=True):
        uppmatrix = self.copy()
        for curij in range(min(uppmatrix.n, uppmatrix.m)):
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
            uppmatrix.clear_zero_rows()
        return uppmatrix

    def gauss_jordan(self):  # TODO: zero rows reduction param
        res = self.upper_triangular()
        for curij in range(res.n - 1, -1, -1):
            for i in range(1, curij + 1):  # Get zeros above 1 (skipped if k == range minimum)
                muler = res[curij - i][curij]
                for j in range(res.m):
                    res[curij - i][j] -= muler * res[curij][j]
        return res

    def inverse_gaussian(self, check_is_invertible=True):
        if check_is_invertible:
            if not self.is_invertible():
                raise Exception("Inverse error: Matrix isn't invertible")
        attached_matrix = Matrix.zeros(self.n, self.n * 2)
        for i in range(self.n):
            for j in range(self.n):
                attached_matrix[i][j] = self[i][j]
        for i in range(self.n):
            attached_matrix[i][i + self.m] = 1

        attached_matrix = attached_matrix.gauss_jordan()

        res = Matrix.zeros(self.n, self.n)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = attached_matrix[i][j + self.n]
        return res

    def rank(self):
        row_echelon = self.upper_triangular()
        return row_echelon.n

    def LU_decomposition(self, compact_return=False):
        if self.n != self.m:
            raise Exception("LU_decomposition error: Matrix isn't square")
        LU = Matrix.zeros(self.n, self.n)
        for i in range(self.n):
            for j in range(i):
                if abs(LU[j][j]) < self.PRECISION:  # U[j][j] == 0
                    raise Exception("LU_decomposition error: U[{}][{}] = 0", j, j)
                s = 0
                for k in range(j):
                    s += LU[i][k] * LU[k][j]
                LU[i][j] = (self[i][j] - s) / LU[j][j]
            for j in range(i, self.n):
                s = 0
                for k in range(i):
                    s += LU[i][k] * LU[k][j]
                LU[i][j] = self[i][j] - s
        if compact_return:
            return LU
        L = Matrix.zeros(self.n, self.n)
        U = Matrix.zeros(self.n, self.n)
        for i in range(self.n):
            for j in range(i):
                L[i][j] = LU[i][j]
            L[i][i] = 1
            for j in range(i, self.n):
                U[i][j] = LU[i][j]
        return L, U


if __name__ == '__main__':
    a = Matrix([[5, 0, 1, 0], [4, 4, 2, 4], [1, 2, 2, 1]], 3, 4)
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
    print("Обратная матрица методом присоединенной матрицы")
    adjugate = c.inverse_adjugate(check_is_invertible=False)
    print(adjugate)
    print("Обратная матрица методом Гаусса-Жордана")
    gaussian = c.inverse_gaussian(check_is_invertible=False)
    print(gaussian)
    print("rank(C) = ", c.rank())
