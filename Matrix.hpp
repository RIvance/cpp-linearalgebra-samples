/**
 *
 * This file is for educational purposes only,
 * the matrix operations are not optimized in any way,
 * please do not use it in a real production environment.
 *
 */

#pragma once

#include <array>
#include <string>
#include <complex>
#include <utility>
#include <algorithm>
#include <functional>

using usize = std::size_t;
using String = std::string;

template <typename Type, usize size>
using Array = std::array<Type, size>;

template <typename Type>
using InitList = std::initializer_list<Type>;

template <typename Signature>
using Function = std::function<Signature>;

String str(const char* s)    { return s; }
String str(const String & s) { return s; }

template <typename Type>
String str(const Type & value)
{
    return std::to_string(value);
}

template <typename FirstType, typename ... Types>
String str(const FirstType & first, const Types & ... args)
{
    return str(first) + std::move(str(args...));
}

String operator * (const String & s, usize nTimes)
{
    String result;
    for (usize i = 0; i < nTimes; i++) {
        result += s;
    }
    return result;
}

class Exception : std::exception
{
  private:

    [[nodiscard]]
    const char *what() const noexcept override
      { return message.data(); }

  protected:

    String message;

  public:

    Exception() : std::exception()
      { /* return */ }

    explicit Exception(String message) : message(std::move(message))
      { /* return */ }
};

class MatrixIndexException : Exception
{
  private:

    explicit MatrixIndexException(String message)
        : Exception(std::move(message))
      { /* return */ }

  public:

    MatrixIndexException() = delete;

    static MatrixIndexException indexOutOfBound(usize rows, usize cols, usize row, usize col)
    {
        return MatrixIndexException(str(
            "The matrix is ", rows, " x ", cols, ".",
            " Index (", row, ", ", col, ") out of bound"
        ));
    }

    static MatrixIndexException rowIndexOutOfBound(usize rows, usize cols, usize row)
    {
        return MatrixIndexException(str(
            "The matrix is ", rows, " x ", cols, ".",
            " Index (row = ", row, ") out of bound"
        ));
    }

    static MatrixIndexException colIndexOutOfBound(usize rows, usize cols, usize col)
    {
        return MatrixIndexException(str(
            "The matrix is ", rows, " x ", cols, ".",
            " Index (col = ", col, ") out of bound"
        ));
    }
};

/**
 * NxM matrix
 * @tparam nRows N
 * @tparam nCols M
 * @tparam Type (Type, +, *)
 */
template <usize nRows, usize nCols, typename Type>
class MatrixData
{
  private:

    Array<Type, nRows * nCols> data = Array<Type, nRows * nCols>();

  public:

    MatrixData() = default;

    Type & operator [] (usize index)
    {
        if (index >= data.size()) {
            throw Exception("Index out of bound");
        }
        return data[index];
    }

    Type & operator () (usize row, usize col = 1)
    {
        usize index = (row - 1) * nCols + (col - 1);
        if (col > nCols || row > nRows || col <= 0 || row <= 0) {
            throw MatrixIndexException::indexOutOfBound(nRows, nCols, row, col);
        }
        return data[index];
    }

    Type operator () (usize row, usize col = 1) const
    {
        usize index = (row - 1) * nCols + (col - 1);
        if (col > nCols || row > nRows || col <= 0 || row <= 0) {
            throw MatrixIndexException::indexOutOfBound(nRows, nCols, row, col);
        }
        return data[index];
    }

    MatrixData & operator = (const Array<Type, nRows * nCols> & newData)
    {
        this->data = newData;
        return *this;
    }
};

template <usize nRows, usize nCols = nRows, typename Type = double>
class Matrix
{
  private:

    MatrixData<nRows, nCols, Type> self = MatrixData<nRows, nCols, Type>();

    template <usize rhsCols>
    using RightMulMatrix = Matrix<nCols, rhsCols, Type>;

    template <usize rhsCols>
    using RightMulResult = Matrix<nRows, rhsCols, Type>;

    explicit Matrix(const MatrixData<nRows, nCols, Type> & data)
      : self(data) { /* return */ }

  public:

    struct Shape
    {
        const usize rows = nRows;
        const usize cols = nCols;
    };

    const usize rows = nRows;
    const usize cols = nCols;

    static const bool isSquare = (nCols == nRows);

    template <typename NewType>
    using SameShape = Matrix<nRows, nCols, NewType>;
    using This = SameShape<Type>;
    using Empty = Matrix<0, 0, Type>;

    using RowVector = Matrix<nRows, 1, Type>;
    using ColVector = Matrix<nCols, 1, Type>;

    Matrix() = default;

    Matrix(InitList<Type> values)
    {
        for (auto iter = values.begin(); iter != values.end(); iter++) {
            self[iter - values.begin()] = *iter;
        }
    }

    explicit Matrix(const Array<Type, nRows * nCols> & data)
    {
        this->self = data;
    }

    static This identity()
    {
        if (!isSquare) {
            throw Exception("Identity matrix must be a square matrix");
        }

        This matrix;
        for (usize i = 1; i < nRows; i++) {
            matrix(i, i) = 1;
        }
        return matrix;
    }

    static This generate(Function<Type(usize, usize)> function)
    {
        This matrix;
        for (usize i = 1; i <= matrix.rows; i++) {
            for (usize j = 1; j <= matrix.cols; j++) {
                matrix(i, j) = function(i, j);
            }
        }
        return matrix;
    }

    /**
     * Get element from matrix
     * @param row 1..nDims, inclusive both
     * @param col 1..nCols, inclusive both
     * @return element at (row, col)
     */

    Type & operator () (usize row, usize col = 1)
    {
        return self(row, col);
    }

    Type operator () (usize row, usize col = 1) const
    {
        return self(row, col);
    }

    Shape shape() const
    {
        return Shape();
    }

    /**
     * Vector length
     * @return
     */
    Type length() const
    {
        if (nCols != 1) {
            throw Exception("Length is only defined for a vector");
        }

        Type squareLength = 0;
        for (usize i = 1; i <= nRows; i++) {
            squareLength += self(i) * self(i);
        }

        return std::sqrt(self);
    }

    void forEach(Function<void(Type)> function) const
    {
        for (usize i = 1; i <= nRows; i++) {
            for (usize j = 1; j <= nCols; j++) {
                function(self(i, j));
            }
        }
    }

    template <typename NewType>
    Matrix<nRows, nCols, NewType> map(Function<NewType(Type)> function) const
    {
        Matrix<nRows, nCols, NewType> result;
        for (usize i = 1; i <= nRows; i++) {
            for (usize j = 1; j <= nCols; j++) {
                result(i, j) = function(self(i, j));
            }
        }
        return result;
    }

    using Comparator = Function<Type(const Type &, const Type &)>;

    Type max(Comparator comparator = std::max) const
    {
        Type result = self(1, 1);
        for (usize i = 1; i <= rows; i++) {
            for (usize j = 1; j <= cols; j++) {
                result = comparator(result, self(i, j));
            }
        }
        return result;
    }

    Type min(Comparator comparator = std::min) const
    {
        return max(comparator);
    }

    String toString(
      #ifdef __clang__
        Function<String(Type)> formatFunc
      #else
        Function<String(Type)> formatFunc = str<Type>
      #endif
    ) const {
        if (rows == 0 || cols == 0) {
            return "";
        }
        String result;

        SameShape<String> stringMatrix = this->map<String>(
            [&formatFunc] (Type value) { return formatFunc(value); }
        );

        usize maxLength = stringMatrix.max([](const String & s1, const String & s2) {
            return s1.length() > s2.length() ? s1 : s2;
        }).length();

        for (usize i = 1; i <= stringMatrix.rows; i++) {
            result += (rows == 1 || (i != 1 && i != rows) ? "|" : i == 1 ? "/" : "\\");
            for (usize j = 1; j <= stringMatrix.cols; j++) {
                String elementString = stringMatrix(i, j);
                result += str(" ") * (maxLength - elementString.length() + 1) + elementString;
            }
            result += " ";
            result += (rows == 1 || (i != 1 && i != rows) ? "|\n" : i == 1 ? "\\\n" : "/\n");
        }

        return result;
    }

    RowVector row(usize index) const
    {
        if (index > rows) {
            throw MatrixIndexException::rowIndexOutOfBound(rows, cols, index);
        }
        RowVector vector;
        for (usize i = 0; i < cols; i++) {
            vector(i) = self(index, i);
        }
        return vector;
    }

    ColVector col(usize index) const
    {
        if (index > cols) {
            throw MatrixIndexException::colIndexOutOfBound(rows, cols, index);
        }
        ColVector vector;
        for (usize i = 0; i < rows; i++) {
            vector(i) = self(i, index);
        }
        return vector;
    }

    void setRow(usize row, const RowVector & rowVector)
    {
        for (usize i = 1; i <= cols; i++) {
            self(row, i) = rowVector(i);
        }
    }

    void setCol(usize col, const ColVector & colVector)
    {
        for (usize i = 1; i <= rows; i++) {
            self(i, col) = colVector(i);
        }
    }

    void swapRow(usize row1, usize row2)
    {
        for (usize i = 1; i <= cols; i++) {
            std::swap(self(row1, i), self(row2, i));
        }
    }

    void swapCol(usize col1, usize col2)
    {
        for (usize i = 1; i <= rows; i++) {
            std::swap(self(i, col1), self(i, col2));
        }
    }

    inline RowVector operator [] (usize index) const
    {
        return row(index);
    }

    Matrix<nRows - 1, nRows == 1 ? 0 : nCols, Type>
    dropRow(usize index)
    {
        Matrix<nRows - 1, nRows == 1 ? 0 : nCols, Type> result;
        for (usize i = 1; i <= result.rows; i++) {
            for (usize j = 1; j <= result.cols; j++) {
                result(i, j) = self(i < index ? i : i + 1, j);
            }
        }
        return result;
    }

    Matrix<nCols == 1 ? 0 : nRows, nCols - 1, Type>
    dropCol(usize index)
    {
        Matrix<nCols == 1 ? 0 : nRows, nCols - 1, Type> result;
        for (usize i = 1; i <= result.rows; i++) {
            for (usize j = 1; j <= result.cols; j++) {
                result(i, j) = self(i, j < index ? j : j + 1);
            }
        }
        return result;
    }

    /**
     * NxM matrix mul MxP matrix -> NxP matrix
     * @tparam nRows   N
     * @tparam nCols   M
     * @tparam rhsCols P
     * @param rhs MxP matrix
     * @return NxP matrix
     */

    template <usize rhsCols>
    RightMulResult<rhsCols> operator * (const RightMulMatrix<rhsCols> & rhs) const
    {
        RightMulResult<rhsCols> result;
        for (usize i = 1; i <= result.rows; i++) {
            for (usize j = 1; j <= result.cols; j++) {
                // (AB)_{ij} = sum_{k=1}^{M} A_{ik} B_{kj}
                Type sum = 0; // (Type, +, *) must have zero element
                for (usize k = 1; k <= cols; k++) {
                    sum += self(i, k) * rhs(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    /**
     * NxM matrix transpose -> MxN matrix
     * @return MxN transpose matrix
     */

    using Transpose = Matrix<nCols, nRows, Type>;

    Transpose transpose() const
    {
        Transpose result;
        for (usize i = 1; i <= rows; i++) {
            for (usize j = 1; j <= cols; j++) {
                result(j, i) = self(i, j);
            }
        }
        return result;
    }

    inline Transpose operator ~ () const
    {
        return transpose();
    }

    This adjugate() const
    {
        return This::generate([this] (usize row, usize col) {
            return this->cofactor(row, col);
        });
    }

    [[nodiscard]]
    bool isInvertible() const
    {
        return isSquare && this->determinant() != 0;
    }

    This inverse() const
    {
        if (nCols != nRows) {
            throw Exception("There is no inverse matrix for a non-square matrix");
        }

        Type det = this->determinant();
        if (det == 0) {
            throw Exception("Determinant equals to zero, no inverse matrix");
        }

        return (1 / det) * this->adjugate();
    }

    /**
     * (N-1)x(M-1) submatrix that removed n-th column and m-th row
     * @param row n
     * @param col m
     * @return (N-1)x(M-1) submatrix
     */

    using Minor = Matrix<
        nRows == 1 ? 1 : nRows - 1,
        nCols == 1 ? 1 : nCols - 1,
    Type>;

    Minor minor(usize row, usize col) const
    {
        if (cols == 1 && rows == 1) {
            return Minor::identity();
        }

        Minor result;
        for (usize i = 1; i <= result.rows; i++) {
            for (usize j = 1; j <= result.cols; j++) {
                // skip `col` and `row`
                result(i, j) = self(i < row ? i : i + 1, j < col ? j : j + 1);
            }
        }
        return result;
    }

    /**
     * determinant of (i, j) minor multiplied by (-1)^{i+j}
     * @param row i
     * @param col j
     * @return cofactor of
     */

    Type cofactor(usize row, usize col) const
    {
        return ((row + col) % 2 == 0 ? 1 : -1) * minor(row, col).determinant();
    }

    /**
     * simplify to a triangular matrix by gauss elimination
     * @return triangular matrix
     */

    This triangular() const
    {
        This result(self);
        /**
         * @var transferRow: which row to eliminate to 1
         */
        for (usize transferRow = 1; transferRow <= rows; transferRow++) {

            // find the longest non-zero row
            usize firstNonZeroCol = 1;
            for (NULL; firstNonZeroCol <= cols; firstNonZeroCol++) {
                usize row = transferRow;
                bool foundNonZeroElement = false;
                for (NULL; row <= rows; row++) {
                    if (result(row, firstNonZeroCol) != 0) {
                        foundNonZeroElement = true;
                        break;
                    }
                }
                if (foundNonZeroElement) {
                    result.swapRow(transferRow, row);
                    break;
                }
            }

            // zero matrix or the following rows are all zeros
            if (firstNonZeroCol > cols) {
                return result;
            }

            double factor1 = result(transferRow, firstNonZeroCol);
            for (usize col = firstNonZeroCol; col <= cols; col++) {
                result(transferRow, col) /= factor1;
            }

            for (usize i = 1; i <= rows; i++) {
                if (i == transferRow) continue;
                double factor2 = result(i, firstNonZeroCol);
                for (usize j = firstNonZeroCol; j <= cols; j++) {
                    result(i, j) = result(i, j) - result(transferRow, j) * factor2;
                }
            }
        }

        return result;
    }

    template <usize rowBegin, usize rowEnd, usize colBegin, usize colEnd>
    Matrix<rowEnd - rowBegin + 1, colEnd - colBegin + 1, Type> submatrix() const
    {
        Matrix<rowEnd - rowBegin + 1, colEnd - colBegin + 1, Type> result;
        for (usize i = rowBegin; i <= rowEnd; i++) {
            for (usize j = colBegin; j <= colEnd; j++) {
                result(i - rowBegin + 1, j - colBegin + 1) = self(i, j);
            }
        }
        return result;
    }

    Type determinant() const
    {
        if (cols != rows) {
            throw Exception("Determinant is not defined for a non-square matrix");
        }

        if (rows == 2) { // ac - bd
            return self(1, 1) * self(2, 2) - self(1, 2) * self(2, 1);
        } else if (rows == 1) {
            return self(1, 1);
        }

        Type result = 0;
        for (usize i = 1; i <= rows; i++) {
            for (usize j = 1; j <= cols; j++) {
                result += self(i, j) * cofactor(i, j);
            }
        }
        return result;
    }

    [[nodiscard]]
    usize rank() const
    {
        usize result = rows;
        for (usize i = 1; i <= rows; i++) {
            bool isAllElementZero = true;
            for (usize j = 1; j <= cols; j++) {
                if (self(i, j) != 0) isAllElementZero = false;
            }
            if (isAllElementZero) result -= 1;
        }
        return result;
    }
};

template <usize nRows, usize nCols, typename Type, typename Scalar>
Matrix<nRows, nCols, Type> operator * (Scalar scalar, const Matrix<nRows, nCols, Type> & matrix)
{
    return (matrix.template map<Type>) (
        [&scalar] (Type value) { return (Type) (scalar * value); }
    );
}
