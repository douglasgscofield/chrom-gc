#ifndef RANDUNIFORM_H
#define RANDUNIFORM_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <ctime>

class RandUniform {
    public:
        RandUniform(bool random_seed = false) : test(false) {
            init(random_seed);
        };
    private:
        bool test;
        double u[98], c, cd, cm;
        int i97, j97;
    public:
        void   init(const bool random_seed = false,
                    const int ij = 1802, const int kl = 9373);
        double draw();
};

inline void RandUniform::init(const bool random_seed, const int ij,
    const int kl)
{
    int iij = ij; 
    int i, j, k, l, ii, jj, m;
    double s, t;
    
    // NOTE: The seed variables can have values between 
    // 0 <= IJ <= 31328 and 0 <= KL <= 30081
    // int iij = 1802; // int kl = 9373;
    if (iij < 0 || iij > 31328 || kl < 0 || kl > 30081) {
        std::cerr << "RandUniform::init: range ij [0, 31328]" << std::endl
            << "RandUniform::init: range kl 0, 30081]" << std::endl;
    }
    assert(iij >= 0 && iij <= 31328 && kl >= 0 && kl <= 30081);
    if (random_seed == true) {
        iij = time(NULL) % 31329;  // system clock seed
    }
    i = (iij/177)%177 + 2;
    j = iij%177 + 2;
    k = (kl/169)%178 + 1;
    l = kl%169;
    for (ii=1; ii<=97; ii++) {
        s = 0.0;
        t = 0.5;
        for (jj=1; jj<=24; jj++) {
            m = (((i*j)%179)*k) % 179;
            i = j;
            j = k;
            k = m;
            l = (53*l + 1) % 169;
            if ((l*m)%64 >= 32) s += t;
            t *= 0.5;
        }
        u[ii] = s;
    }
    c = 362436.0 / 16777216.0;
    cd = 7654321.0 / 16777216.0;
    cm = 16777213.0 / 16777216.0;
    i97 = 97;
    j97 = 33;
    test = true;
}

inline double RandUniform::draw()
{
    double uni;

    if (test == false) {
        std::cerr << "RandUniform::draw: Call init() first." << std::endl;
    }
    assert(test == true);
	uni = u[i97] - u[j97];
	if (uni < 0.0) uni += 1.0;
	u[i97] = uni;
	i97--;
	if (i97==0) i97 = 97;
	j97--;
	if (j97==0) j97 = 97;
	c -= cd;
	if (c<0.0) c += cm;
	uni -= c;
	if (uni<0.0) uni += 1.0;
	return uni;
}

#endif // RANDUNIFORM_H

