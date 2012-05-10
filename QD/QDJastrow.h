/* 
 * File:   QDJastrow.h
 * Author: zigg
 *
 * Created on May 10, 2012, 9:33 AM
 */

#ifndef QDJASTROW_H
#define	QDJASTROW_H

#include "../Jastrow.h"

class QDJastrow : public Jastrow {
public:
    QDJastrow();
    QDJastrow(int, int, double);
    QDJastrow(const QDJastrow& orig);
    virtual ~QDJastrow();

    virtual double getVariationalGradient(const mat &);
    virtual double evaluate(const mat &);
    virtual void computeGradient(const mat &, int);
    virtual double getLaplacian(const mat &, int);
    double getRatio(const mat &, const mat &);

    virtual void setNewBeta(double beta) {
        this->beta = beta;
    }
private:
    mat a;
};

#endif	/* QDJASTROW_H */

