import random


class Matrix:
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
                res.matrix[i][j] = self.matrix[i][j] + other.matrix[i][j]
        return res

    def __sub__(self, other):  # self - other
        if self.n != other.n or self.m != other.m:
            raise Exception("Sub error: Matrix sizes aren't equal")
        res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res.matrix[i][j] = self.matrix[i][j] - other.matrix[i][j]
        return res

    def __mul__(self, other):  # self x other
        if isinstance(other, float) or isinstance(other, int):
            res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
            for i in range(self.n):
                for j in range(self.m):
                    res.matrix[i][j] = self.matrix[i][j] * other
            return res
        else:
            return self.dot(other)

    def __truediv__(self, other):  # self / other (only for numbers)
        if not isinstance(other, float) and not isinstance(other, int):
            raise Exception("Div error: " + type(other) + " isn't a number")
        res = Matrix([[0 for j in range(self.m)] for i in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res.matrix[i][j] = self.matrix[i][j] / other
        return res

    def __eq__(self, other):  # self == other
        if (self.n != other.n) or (self.m != other.m):
            return False
        for i in range(self.n):
            for j in range(self.m):
                if self.matrix[i][j] != other.matrix[i][j]:
                    return False
        return True

    def __ne__(self, other):  # self != other
        return not self.__eq__(other)

    def __str__(self):
        s = ""
        for i in range(self.n):
            s += self.matrix[i].__str__() + "\n"
        return s

    def dot(self, other):  # self x other_matrix
        if self.m != other.n:
            raise Exception("Dot error: Number of columns isn't equal to number of rows")
        res = Matrix([[0 for j in range(other.m)] for i in range(self.n)], self.n, other.m)
        for i in range(self.n):
            for j in range(other.m):
                res.matrix[i][j] = 0
                for k in range(self.m):
                    res.matrix[i][j] += self.matrix[i][k] * other.matrix[k][j]
        return res

    def transpose(self):  # creates new matrix without overwriting an existing one
        res = Matrix([[0 for j in range(self.n)] for i in range(self.m)], self.m, self.n)
        for i in range(self.m):
            for j in range(self.n):
                res.matrix[i][j] = self.matrix[j][i]
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
        minor = Matrix._deliRowjCol(matrix, i, j)
        return Matrix._det(minor)

    @staticmethod
    def _deliRowjCol(matrix, i, j):
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
        return False if self.det() == 0 or (self.n != self.m) else True

    def inverse_adjugate(self):
        if not self.is_invertible():
            raise Exception("Inverse error: Matrix isn't invertible")
        transposed = self.transpose()
        cofactor_matrix = Matrix([[0 for t in range(self.n)] for k in range(self.n)], self.n, self.n)
        for i in range(self.n):
            for j in range(self.n):
                cofactor_matrix.matrix[i][j] = Matrix._minor(transposed.matrix, i, j) * ((-1) ** (2 + i + j))
        res = cofactor_matrix * (1 / self.det())
        return res

    def upper_triangular(self):
        uppmatrix = Matrix([[0 for t in range(self.m)] for k in range(self.n)], self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                uppmatrix.matrix[i][j] = self.matrix[i][j]
        for curij in range(uppmatrix.n):
            if uppmatrix.matrix[curij][curij] == 0:  # Get nonzero element at top
                for i in range(uppmatrix.n - curij):
                    if uppmatrix.matrix[curij + i][curij] != 0:
                        uppmatrix.matrix[curij], uppmatrix.matrix[curij + i] = \
                            uppmatrix.matrix[curij + i], uppmatrix.matrix[curij]
                        break

            diver = uppmatrix.matrix[curij][curij]
            for j in range(self.m):  # Get 1 element at top
                uppmatrix.matrix[curij][j] /= diver

            for i in range(1, uppmatrix.n - curij):  # Get zeros below 1 (skipped if k == range maximum)
                muler = uppmatrix.matrix[curij + i][curij]
                for j in range(uppmatrix.m):
                    uppmatrix.matrix[curij + i][j] -= muler * uppmatrix.matrix[curij][j]
        return uppmatrix

    def inverse_gaussian(self):
        if not self.is_invertible():
            raise Exception("Inverse error: Matrix isn't invertible")
        attached_matrix = Matrix([[0 for t in range(self.n * 2)] for k in range(self.n)], self.n, self.n * 2)
        for i in range(self.n):
            for j in range(self.n):
                attached_matrix.matrix[i][j] = self.matrix[i][j]
        for i in range(self.n):
            attached_matrix.matrix[i][i + self.m] = 1
        attached_matrix = attached_matrix.upper_triangular()
        #print(attached_matrix)
        # for i in range(attached_matrix.n): TODO: make zero rows and column check
        #     for j in range(attached_matrix.m):
        #         if attached_matrix.matrix[i][j] != 0:
        #             break
        #         if j == attached_matrix.m - 1
        for curij in range(attached_matrix.n - 1, -1, -1):
            for i in range(1, curij + 1):  # Get zeros above 1 (skipped if k == range minimum)
                muler = attached_matrix.matrix[curij - i][curij]
                for j in range(attached_matrix.m):
                    attached_matrix.matrix[curij - i][j] -= muler * attached_matrix.matrix[curij][j]
        #print(attached_matrix)
        res = Matrix([[0 for t in range(self.n)] for k in range(self.n)], self.n, self.n)
        for i in range(self.n):
            for j in range(self.m):
                res.matrix[i][j] = attached_matrix.matrix[i][j + self.n]
        return res


if __name__ == '__main__':
    a = Matrix([[random.randint(0, 5) for j in range(5)] for i in range(4)], 4, 5)
    print(a)
    b = a.transpose()
    print(b)
    c = a * b
    print(c)
    print(c.det())
    print(c.is_invertible())
    print(c.inverse_adjugate())
    print(c.inverse_gaussian())
    d = Matrix([[2, 5, 7], [6, 3, 4], [5, -2, -3]], 3, 3)
    print(d)
    print(d.inverse_adjugate())
