/* 
 * File:   Orbital.cpp
 * Author: zigg
 * 
 * Created on May 10, 2012, 8:43 AM
 */

#include "Orbital.h"

Orbital::Orbital() {
}

Orbital::Orbital(int dim, double alpha) : dim(dim), alpha(alpha) {
    H = new Hermite();
#if 0
    cout
            << "\t dim = " << this->dim 
            << "\t alpha = " << this->alpha 
            << endl;
#endif
}

Orbital::Orbital(const Orbital& orig) {
}

Orbital::~Orbital() {
    delete H;
}

