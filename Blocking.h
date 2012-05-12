/* 
 * File:   Blocking.h
 * Author: sigve
 *
 * Created on 11. mai 2012, 16:06
 */

#ifndef BLOCKING_H
#define	BLOCKING_H

#include <armadillo>
using namespace arma;

class Blocking {
public:
    Blocking();
    Blocking(int);
    Blocking(const Blocking& orig);
    virtual ~Blocking();
    vec block(vec, int, int);
private:

};

#endif	/* BLOCKING_H */

