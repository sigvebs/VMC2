/* 
 * File:   ComputeOneWF.h
 * Author: zigg
 *
 * Created on May 9, 2012, 10:51 PM
 */

#ifndef COMPUTEONEWF_H
#define	COMPUTEONEWF_H

class ComputeOneWF {
public:
    ComputeOneWF();
    ComputeOneWF(const ComputeOneWF& orig);
    virtual ~ComputeOneWF();
private:
    double E;
    double Esq;
    double accepted;
};

#endif	/* COMPUTEONEWF_H */

