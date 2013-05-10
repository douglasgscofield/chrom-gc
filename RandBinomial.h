#ifndef BINOM_H
#define BINOM_H
#include "RandUniform.h"

#include <iostream>
#include <cassert>
#include <cmath>

class RandBinomial {
    /* long ignbin(long n,double pp)
                    GENerate BINomial random deviate Function
     Generates a single random deviate from a binomial
     distribution whose number of trials is N and whose
     probability of an event in each trial is P.
                              Arguments
     n  --> The number of trials in the binomial distribution
            from which a random deviate is to be generated.
     p  --> The probability of an event in each trial of the
            binomial distribution from which a random deviate
            is to be generated.
     draw <-- A random deviate yielding the number of events
                from N independent trials, each of which has
                a probability of event P.
                              Method
     This is algorithm BTPE from:
         Kachitvichyanukul, V. and Schmeiser, B. W.
         Binomial Random Variate Generation.
         Communications of the ACM, 31, 2
         (February, 1988) 216.
         */
    public:
        RandBinomial();
    private:
        RandUniform uniform;
        double psave;
        long nsave;
        long ignbin, i, ix, ix1, k, m, mp, T1;
        double al, alv, amaxp, c, f, f1, f2, ffm, fm, g, p, p1, 
               p2, p3, p4, q, qn, r, u, v, w, w2, x, x1, x2, xl,
               xll, xlr, xm, xnp, xnpq, xr, ynorm, z, z2;
    public:
        long draw(const long n, const double pp);
    private:
        template<class T> inline const T ABS(const T a) const {
            return((a) > ((T)0) ? (a) : -(a));
        }
        template<class T> inline const T MAX(const T& a, const T& b) const {
            return((b) > (a) ? (b) : (a)); 
        }
        template<class T> inline const T MIN(const T& a, const T& b) const {
            return((b) < (a) ? (b) : (a)); 
        }
};

inline RandBinomial::RandBinomial() 
    : psave(-1.0), nsave(-1)
{ 
    /* empty */
}

inline long RandBinomial::draw(const long n, const double pp) 
{
    if(pp != psave) goto S10;
    if(n != nsave) goto S20;
    if(xnp < 30.0) goto S150;
    goto S30;
S10:
/*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
*/
    psave = pp;
    p = MIN(psave,1.0-psave);
    q = 1.0-p;
S20:
    xnp = n*p;
    nsave = n;
    if(xnp < 30.0) goto S140;
    ffm = xnp+p;
    m = (long) ffm;
    fm = m;
    xnpq = xnp*q;
    p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
    xm = fm+0.5;
    xl = xm-p1;
    xr = xm+p1;
    c = 0.134+20.5/(15.3+fm);
    al = (ffm-xl)/(ffm-xl*p);
    xll = al*(1.0+0.5*al);
    al = (xr-ffm)/(xr*q);
    xlr = al*(1.0+0.5*al);
    p2 = p1*(1.0+c+c);
    p3 = p2+c/xll;
    p4 = p3+c/xlr;
S30:
/*
*****GENERATE VARIATE
*/
    u = uniform.draw()*p4;
    v = uniform.draw();
/*
     TRIANGULAR REGION
*/
    if(u > p1) goto S40;
    ix = (long) (xm-p1*v+u);
    goto S170;
S40:
/*
     PARALLELOGRAM REGION
*/
    if(u > p2) goto S50;
    x = xl+(u-p1)/c;
    v = v*c+1.0-ABS(xm-x)/p1;
    if(v > 1.0 || v <= 0.0) goto S30;
    ix = (long) x;
    goto S70;
S50:
/*
     LEFT TAIL
*/
    if(u > p3) goto S60;
    ix = (long) (xl+log(v)/xll);
    if(ix < 0) goto S30;
    v *= ((u-p2)*xll);
    goto S70;
S60:
/*
     RIGHT TAIL
*/
    ix = (long) (xr-log(v)/xlr);
    if(ix > n) goto S30;
    v *= ((u-p3)*xlr);
S70:
/*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
    k = ABS(ix-m);
    if(k > 20 && k < xnpq/2-1) goto S130;
/*
     EXPLICIT EVALUATION
*/
    f = 1.0;
    r = p/q;
    g = (n+1)*r;
    T1 = m-ix;
    if(T1 < 0) goto S80;
    else if(T1 == 0) goto S120;
    else  goto S100;
S80:
    mp = m+1;
    for(i=mp; i<=ix; i++) f *= (g/i-r);
    goto S120;
S100:
    ix1 = ix+1;
    for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
    if(v <= f) goto S170;
    goto S30;
S130:
/*
     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
*/
    amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
    ynorm = -(k*k/(2.0*xnpq));
    alv = log(v);
    if(alv < ynorm-amaxp) goto S170;
    if(alv > ynorm+amaxp) goto S30;
/*
     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
     THE FINAL ACCEPTANCE/REJECTION TEST
*/
    x1 = ix+1.0;
    f1 = fm+1.0;
    z = n+1.0-fm;
    w = n-ix+1.0;
    z2 = z*z;
    x2 = x1*x1;
    f2 = f1*f1;
    w2 = w*w;
    if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
      (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
    goto S30;
S140:
/*
     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
*/
    qn = pow(q,(double)n);
    r = p/q;
    g = r*(n+1);
S150:
    ix = 0;
    f = qn;
    u = uniform.draw();
S160:
    if(u < f) goto S170;
    if(ix > 110) goto S150;
    u -= f;
    ix += 1;
    f *= (g/ix-r);
    goto S160;
S170:
    if(psave > 0.5) ix = n-ix;
    ignbin = ix;
    return ignbin;
}

#endif // BINOM_H

