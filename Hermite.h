/* 
 * File:   Hermite.h
 * Author: zigg
 *
 * Created on May 10, 2012, 9:04 AM
 */

#ifndef HERMITE_H
#define	HERMITE_H

class Hermite {
public:
    Hermite();
    Hermite(const Hermite& orig);
    virtual ~Hermite();
    double polynomial(const int, const double);
private:

};

#endif	/* HERMITE_H */

