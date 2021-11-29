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
#include <utility>
#include <functional>

using usize = std::size_t;
using String = std::string;


template <typename Type, usize size>
using Array = std::array<Type, size>;

template <typename Signature>
using Function = std::function<Signature>;

auto str(const char* s)
{
    return String(s);
}

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

    Array<Type, nRows * nCols> data;

  public:

    Type & operator [] (usize index)
    {
        if (index >= data.size()) {
            throw Exception("Index out of bound");
        }
        return data[index];
    }

    Type & operator () (usize row, usize col)
    {
        usize index = (row - 1) * nCols + (col - 1);
        if (col > nCols || row > nRows || col <= 0 || row <= 0) {
            throw MatrixIndexException::indexOutOfBound(nRows, nCols, row, col);
        }
        return data[index];
    }

    Type operator () (usize row, usize col) const
    {
        usize index = (row - 1) * nCols + (col - 1);
        if (col > nCols || row > nRows || col <= 0 || row <= 0) {
            throw MatrixIndexException::indexOutOfBound(nRows, nCols, row, col);
        }
        return data[index];
    }
};

template <usize nRows, usize nCols, typename Type = double>
class Matrix
{
  private:

    MatrixData<nRows, nCols, Type> self;

    template <usize rhsCols>
    using RightMulMatrix = Matrix<nCols, rhsCols, Type>;

    template <usize rhsCols>
    using RightMulResult = Matrix<nRows, rhsCols, Type>;

  public:

    struct Shape
    {
        const usize rows = nRows;
        const usize cols = nCols;
    };

    const usize rows = nRows;
    const usize cols = nCols;

    template <typename NewType>
    using SameShape = Matrix<nRows, nCols, NewType>;

    using RowVector = Matrix<nRows, 1, Type>;
    using ColVector = Matrix<nCols, 1, Type>;

    Matrix() = default;

    /**
     * Get element from matrix
     * @param row 1..nRows, inclusive both
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
        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                result = comparator(result, self(i, j));
            }
        }
        return result;
    }

    Type min(Comparator comparator = std::min) const
    {
        return max(comparator);
    }

    String toString(Function<String(Type)> formatFunc = str<Type>) const
    {
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

    String toLatexString()
    {
        // TODO
    }

    RowVector row(usize index) const
    {
        if (index > rows) {
            throw MatrixIndexException::rowIndexOutOfBound(rows, cols, index);
        }
        RowVector vector;
        for (int i = 0; i < cols; i++) {
            vector(i) = self(index, i);
        }
    }

    ColVector col(usize index) const
    {
        if (index > cols) {
            throw MatrixIndexException::colIndexOutOfBound(rows, cols, index);
        }
        ColVector vector;
        for (int i = 0; i < rows; i++) {
            vector(i) = self(i, index);
        }
    }

    inline RowVector operator [] (usize index) const
    {
        return row(index);
    }

    Matrix<nRows - 1, nCols, Type> dropRow(usize index)
    {
        Matrix<nRows - 1, nCols, Type> result;
        for (int i = 1; i <= result.rows; i++) {
            for (int j = 1; j <= result.cols; j++) {
                result(i, j) = self(i < index ? i : i + 1, j);
            }
        }
        return result;
    }

    Matrix<nRows, nCols - 1, Type> dropCol(usize index)
    {
        Matrix<nRows, nCols - 1, Type> result;
        for (int i = 1; i <= result.rows; i++) {
            for (int j = 1; j <= result.cols; j++) {
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

    /**
     * (N-1)x(M-1) submatrix that removed n-th column and m-th row
     * @param row n
     * @param col m
     * @return (N-1)x(M-1) submatrix
     */

    using Cofactor = Matrix<nRows - 1, nCols - 1, Type>;

    Cofactor cofactor(usize row, usize col) const
    {
        Cofactor result;
        for (usize i = 1; i <= result.rows; i++) {
            for (usize j = 1; j <= result.cols; j++) {
                // skip `col` and `row`
                result(i, j) = self(i < row ? i : i + 1, j < col ? j : j + 1);
            }
        }
        return result;
    }

    template <usize rowBegin, usize rowEnd, usize colBegin, usize colEnd>
    Matrix<rowEnd - rowBegin + 1, colEnd - colBegin + 1, Type> submatrix() const
    {
        // TODO
    }

    Type determinant() const
    {
        // TODO
    }
};

template <usize nRows, typename Type = double>
using Vector = Matrix<nRows, 1, Type>;