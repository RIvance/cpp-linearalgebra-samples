#pragma once

#include "Matrix.hpp"
#include "Vectors.hpp"

#include <valarray>
#include <iostream>

ulong gcd(ulong a, ulong b)
{
    if (b != 0) return gcd(b, a % b);
    else return a;
}

class Fraction
{
  public:
    uint numerator;
    uint denominator;
    bool isNegative;

    [[nodiscard]]
    String toString() const
    {
        return (
            numerator == 0 ? "0" :
            denominator == 1 ? str(isNegative ? "-" : "", numerator) :
            str(isNegative ? "-" : "", numerator, "/", denominator)
        );
    }

    [[nodiscard]]
    String toLatex() const
    {
        return (
            numerator == 0 ? "0" :
            denominator == 1 ? str(isNegative ? "-" : "", numerator) :
            str(isNegative ? "-\\frac{" : "\\frac{", numerator, "}{", denominator, "}")
        );
    }

    static Fraction from(double value, usize cycles = 10, double precision = 5e-4)
    {
        Fraction fraction {};
        fraction.isNegative = value < 0;
        value = value >= 0 ? value : -value;

        double newNumber, wholePart;
        double decimalPart = value - (uint) value;

        std::valarray<double> vec1 = { double((uint) value), 1};
        std::valarray<double> vec2 = { 1, 0 };


        for (int i = 0; decimalPart > precision & i < cycles; i++) {
            newNumber = 1 / decimalPart;
            wholePart = (uint) newNumber;

            auto tmp = vec1;
            vec1 = wholePart * vec1 + vec2;
            vec2 = tmp;

            decimalPart = newNumber - wholePart;
        }

        fraction.numerator   = (uint) vec1[0];
        fraction.denominator = (uint) vec1[1];

        return fraction;
    }

    static String doubleToLatexFrac(double value)
    {
        return Fraction::from(value).toLatex();
    }

    static String doubleToStringFrac(double value)
    {
        return Fraction::from(value).toString();
    }
};

namespace MatrixUtils
{
    template <usize nRows, usize nCols, typename ValueType>
    void printValueMatrix(const Matrix<nRows, nCols, ValueType> & matrix)
    {
        std::cout << matrix.toString(Fraction::doubleToStringFrac);
    }

    template <usize nRows, usize nCols, typename Type>
    String matrixToLatex(
        const Matrix<nRows, nCols, Type> & matrix,
        Function<String(Type)> formatFunc = Fraction::doubleToLatexFrac
    ) {
        String result = str("\\left[\\begin{array}{", str("c") * matrix.cols, "}");

        for (usize i = 1; i <= nRows; i++) {
            for (usize j = 1; j <= nCols; j++) {
                result += formatFunc(matrix(i, j)) + " & ";
            }
            result += i == nRows ? "" : "\\\\ ";
        }

        return str(result, "\\end{array}\\right]");
    }
}