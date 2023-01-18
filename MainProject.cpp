#define _USE_MATH_DEFINES
#define PRINT true

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>
#include <complex>
#include <utility>
#include <algorithm>
#include "except.h"

typedef float fp_t;

using namespace std;

// ================================================================== //

template<typename fp_t>
inline bool isZero(const fp_t& x)
{
    return FP_ZERO == fpclassify(x);
}

template<typename fp_t>
inline fp_t fms(fp_t a, fp_t b, fp_t c, fp_t d)
{
    fp_t cd = -c * d;

    return fma(a, b, cd) - fma(c, d, cd);
}

template <typename fp_t>
inline int sgn(fp_t x)
{
    return (fp_t(0) < x) - (x < fp_t(0));
}

template<typename fp_t>
vector<complex<fp_t>> cubeRoot(const complex<fp_t>& z)
{
    fp_t r = abs(z);
    fp_t phi = arg(z);

    fp_t three = 3.L;

    return
    {
        complex<fp_t>(polar(cbrt(r), phi / three)),
        complex<fp_t>(polar(cbrt(r), fma(static_cast<fp_t>(2.L), static_cast<fp_t>(M_PI), phi) / three)),
        complex<fp_t>(polar(cbrt(r), fma(static_cast<fp_t>(4.L), static_cast<fp_t>(M_PI), phi) / three)) // *исправлено
    };
}

template<typename fp_t>
complex<fp_t> fmac(complex<fp_t> x, complex<fp_t> y, complex<fp_t> z)
{
    fp_t r = fma(-x.imag(), y.imag(), fma(x.real(), y.real(), z.real()));
    fp_t i = fma(y.real(), x.imag(), fma(x.real(), y.imag(), z.imag()));

    return complex<fp_t>(r, i);
}

template<typename fp_t>
complex<fp_t> fmsc(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d)
{
    complex<fp_t> cd = -c * d;

    return fmac(a, b, cd) - fmac(c, d, cd);
}

template<typename fp_t>
complex<fp_t> epsilonC(complex<fp_t> x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() > abs(x.imag()) ? complex<fp_t>(x.real(), 0) : x;
}

template<typename fp_t>
complex<fp_t> cubicPolynomial(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d, complex<fp_t> x)
{
    return fmac(fmac(fmac(a, x, b), x, c), x, d);
}

template<typename fp_t>
vector<complex<fp_t>> rootSearch(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d, vector<complex<fp_t>>& probRoots)
{
    // подставляем полученные корни в исходное уравнение
    complex<fp_t> checkC = cubicPolynomial(a, b, c, d, probRoots[0]);
    complex<fp_t> checkC_ = cubicPolynomial(a, b, c, d, probRoots[3]);

    vector<complex<fp_t>> roots(3);

    vector<pair<complex<fp_t>, complex<fp_t>>> temp(6);

    int count = 0;

    for (size_t i = 0; i < 6; ++i)
    {
        if (!isnan(probRoots[i].real()))
        {
            temp[count] = pair<complex<fp_t>, complex<fp_t>>(epsilonC(cubicPolynomial(a, b, c, d, probRoots[i])), probRoots[i]);
            ++count;
        }
    }

    if (count == 6)
    {
        sort(temp.begin(), temp.end(), [](const pair<complex<fp_t>, complex<fp_t>>& left, const pair<complex<fp_t>, complex<fp_t>>& right)
            {
                return abs(left.first) < abs(right.first);
            });
    }

    roots =
    {
        epsilonC(temp[0].second),
        epsilonC(temp[1].second),
        epsilonC(temp[2].second)
    };

    return roots;
}

// ================================================================== //

template<typename fp_t>
unsigned int solveCubic(fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // нормировка коэффициентов
    if (isZero(a) || isinf(b /= a))
        return 0;
    if (isinf(c /= a))
        return 0;
    if (isinf(d /= a))
        return 0;
    a = 1;
    // числовые константы
    const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);
    const fp_t ONE_TWENTY_SEVEN = static_cast<fp_t>(1.0L / 27.0L);

    complex<fp_t> ONE_HALF_C(0.5L, 0);
    complex<fp_t> ONE_THIRD_C(ONE_THIRD, 0);

    // коэффициенты полинома в комплексной обертке
    complex<fp_t> A_C(1, 0);
    complex<fp_t> B_C(b, 0);
    complex<fp_t> C_C(c, 0);
    complex<fp_t> D_C(d, 0);

    // кол-во вещественных корней
    unsigned numOfRoots = 0;

    // расчетные коэффициенты
    complex<fp_t> C1 = fms(static_cast<fp_t>(3), c, b, b);
    complex<fp_t> C2 = fms(b, fms(static_cast<fp_t>(2) * b, b, static_cast<fp_t>(9), c), static_cast<fp_t>(-27), d);

    complex<fp_t> CCC = fmsc(ONE_HALF_C, sqrt(fmsc(static_cast<fp_t>(4) * C1 * C1, C1, -C2, C2)), ONE_HALF_C, C2);
    complex<fp_t> CCC_ = fmsc(ONE_HALF_C, sqrt(fmsc(static_cast<fp_t>(4) * C1 * C1, C1, -C2, C2)), -ONE_HALF_C, C2);

    vector<complex<fp_t>> C = cubeRoot(CCC);
    vector<complex<fp_t>> C_ = cubeRoot(CCC_);

    // возможные корни уравнение, которые следует уточнить
    vector<complex<fp_t>> probRoots(6);
    probRoots =
    {
        fmac(ONE_THIRD_C, -C1 / C[0], fmsc(ONE_THIRD_C, C[0], ONE_THIRD_C, B_C)),
        fmac(ONE_THIRD_C, -C1 / C[1], fmsc(ONE_THIRD_C, C[1], ONE_THIRD_C, B_C)),
        fmac(ONE_THIRD_C, -C1 / C[2], fmsc(ONE_THIRD_C, C[2], ONE_THIRD_C, B_C)),
        fmac(ONE_THIRD_C, -C1 / C_[0], fmsc(ONE_THIRD_C, C_[0], ONE_THIRD_C, B_C)),
        fmac(ONE_THIRD_C, -C1 / C_[1], fmsc(ONE_THIRD_C, C_[1], ONE_THIRD_C, B_C)),
        fmac(ONE_THIRD_C, -C1 / C_[2], fmsc(ONE_THIRD_C, C_[2], ONE_THIRD_C, B_C)),
    };

    // ищем истинные корни и присваиваем их вектору
    vector<complex<fp_t>> rootsC = rootSearch(A_C, B_C, C_C, D_C, probRoots);

    // проверка на наличие мнимой части
    if (isZero(rootsC[0].imag()) && !isnan(rootsC[0].real()) && !isnan(rootsC[0].imag()))
    {
        roots[0] = rootsC[0].real();
        ++numOfRoots;
    }
    if (isZero(rootsC[1].imag()) && !isnan(rootsC[1].real()) && !isnan(rootsC[1].imag()))
    {
        roots[numOfRoots] = rootsC[1].real();
        ++numOfRoots;
    }
    if (isZero(rootsC[2].imag()) && !isnan(rootsC[2].real()) && !isnan(rootsC[2].imag()))
    {
        roots[numOfRoots] = rootsC[2].real();
        ++numOfRoots;
    }

    return numOfRoots;
}

// ================================================================== //

template <typename fp_t>
void testCubicAdv(const int testCount, const fp_t dist)
{
    int P = 3; // power, total number of tests
    fp_t low = -1, high = 1; // [low, high], max distance between clustered roots
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0;
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;
    int count = 0;

    vector<fp_t> coefficients(P + 1);
    vector<fp_t> trueRoots(P);

    for (size_t i = 0; i < testCount; ++i)
    {
        vector<fp_t> foundRoots(P);

        generate_polynomial<fp_t>(P, 0, P, 0, dist, low, high, trueRoots, coefficients);
        numOfFoundRoots = solveCubic<fp_t>(coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        compare_roots<fp_t>(numOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError);

        if (isinf(absMaxError))
            cantFind += 1;
        else
        {
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;

            count += relMaxError > 1 ? 1 : 0;
        }
    }

    if (PRINT)
    {
        cout << "\n\n\t\t\tCUBIC TEST RESULTS\n\n";
        cout << "Max distance: " << dist << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "Mean absMaxError = " << absErrors / (testCount - cantFind) << endl;
        cout << "Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: " << maxAbsAllofTest << endl;
        cout << "Mean RelMaxError = " << relError / (testCount - cantFind) << endl;
        cout << "Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: " << maxRelAllofTest << endl;
        cout << "RelMaxError > 1: " << count << " times" << endl;
    }
}

int main()
{
    testCubicAdv(1000000, static_cast<fp_t>(1e-05));

    system("pause");
    return 0;
}
