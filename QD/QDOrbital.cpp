/* 
 * File:   QDOrbital.cpp
 * Author: Sigve
 * 
 * Created on May 10, 2012, 8:48 AM
 */

#include "QDOrbital.h"

////////////////////////////////////////////////////////////////////////////////

QDOrbital::QDOrbital() {
}

////////////////////////////////////////////////////////////////////////////////

QDOrbital::QDOrbital(const QDOrbital& orig) {
}

////////////////////////////////////////////////////////////////////////////////

QDOrbital::~QDOrbital() {
}

////////////////////////////////////////////////////////////////////////////////

QDOrbital::QDOrbital(int dim, double alpha, double w) : w(w), Orbital(dim, alpha) {
    sqrtWAlpha = sqrt(w * alpha);
    wAlpha = w*alpha;
    sqrtW = sqrt(w);
    sqrtWDividedAlpha = sqrt(w / alpha);
    gradient = zeros(1, dim);
#if 0 
    cout
            << "\talpha = " << alpha
            << "\t w = " << w
            << "\t sqrtWAlpha = " << sqrtWAlpha
            << "\t wAlpha = " << wAlpha
            << "\t sqrtW = " << sqrtW
            << "\t sqrtWDividedAlpha = " << sqrtWDividedAlpha
            << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////

double QDOrbital::evaluate(const rowvec &r, const int nx, const int ny) {
    double Hx = H->polynomial(nx, sqrtWAlpha * r(0));
    double Hy = H->polynomial(ny, sqrtWAlpha * r(1));

    return Hx * Hy * exp(evaluateExp(r));
}

////////////////////////////////////////////////////////////////////////////////

double QDOrbital::evaluateExp(const rowvec &r) {
    return -0.5 * wAlpha * dot(r, r);
}

////////////////////////////////////////////////////////////////////////////////

double QDOrbital::evaluateLaplacian(const rowvec &r, const int nx, const int ny) {
    double sum = 0;

    double Hx = H->polynomial(nx, sqrtWAlpha * r(0));
    double Hy = H->polynomial(ny, sqrtWAlpha * r(1));

    if (nx > 1)
        sum += 4.0 * nx * (nx - 1) * H->polynomial(nx - 2, sqrtWAlpha * r(0)) / Hx;

    if (ny > 1)
        sum += 4.0 * ny * (ny - 1) * H->polynomial(ny - 2, sqrtWAlpha * r(1)) / Hy;

    if (nx > 0)
        sum -= 4.0 * sqrtWAlpha * nx * r(0) * H->polynomial(nx - 1, sqrtWAlpha * r(0)) / Hx;

    if (ny > 0)
        sum -= 4 * sqrtWAlpha * ny * r(1) * H->polynomial(ny - 1, sqrtWAlpha * r(1)) / Hy;

    sum += wAlpha * dot(r, r) - 2;
    sum *= wAlpha;

    return sum * evaluate(r, nx, ny);
}

////////////////////////////////////////////////////////////////////////////////

rowvec QDOrbital::getGradient(const rowvec &r, const int nx, const int ny) {
    for (int i = 0; i < dim; i++) {
        gradient(i) = 0;
    }

    double Hx = H->polynomial(nx, sqrtWAlpha * r(0));
    double Hy = H->polynomial(ny, sqrtWAlpha * r(1));

    // x-component.
    if (nx > 0)
        gradient(0) = 2 * nx * sqrtWAlpha * H->polynomial(nx - 1, sqrtWAlpha * r(0)) / Hx;

    // y-component.
    if (ny > 0)
        gradient(1) = 2 * ny * sqrtWAlpha * H->polynomial(ny - 1, sqrtWAlpha * r(1)) / Hy;

    gradient -= wAlpha*r;

    return gradient * evaluate(r, nx, ny);
}

////////////////////////////////////////////////////////////////////////////////

double QDOrbital::variationalDerivative(const rowvec &r, const int nx, const int ny) {
    double derivative;
    double r_sq = dot(r, r);

    double Hx = H->polynomial(nx, sqrtWAlpha * r(0));
    double Hy = H->polynomial(ny, sqrtWAlpha * r(1));

    // Fist part of the derivative.
    derivative = -0.5 * w * r_sq;

    if (nx > 0)
        derivative += r(0) * nx * sqrtWDividedAlpha * H->polynomial(nx - 1, sqrtWAlpha * r(0)) / Hx;

    if (ny > 0)
        derivative += r(1) * ny * sqrtWDividedAlpha * H->polynomial(ny - 1, sqrtWAlpha * r(1)) / Hy;

    return derivative * evaluate(r, nx, ny);
}

////////////////////////////////////////////////////////////////////////////////

void QDOrbital::setNewAlpha(double newAlpha) {
    alpha = newAlpha;
    sqrtWAlpha = sqrt(w * newAlpha);
    wAlpha = w*newAlpha;
    sqrtWDividedAlpha = sqrt(w / newAlpha);
};