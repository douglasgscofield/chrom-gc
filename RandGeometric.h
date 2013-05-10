#ifndef __RANDGEOMETRIC_H__
#define __RANDGEOMETRIC_H__

// class RandGeometric
//
// Creates a geometrically-distributed random deviate.  The geometric
// distribution is parameterized by a single parameter, p, which is the
// probability of a success in a trial.  This function returns the number of
// trials k until a success.  The probability for each k is p*[(1-p)^(k-1)],
// and a uniform random deviate is used to determine when exactly the success
// occurs; k is then calculated and returned.

// Routines here are derived from the (G)nu (S)cientific (L)ibrary, version
// 1.8.  Routine-specific copyright notices from the GSL are provided within
// the class definition, though the paragraphs common to each copyright notice
// and to GPL-covered code in general are reproduced below.

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <cmath>
#include <ctime>
#include <cassert>
#include "RandUniform.h"

class RandGeometric {
    private:
        RandUniform    unif;
        double         prob;  // the probability until the first success
        double         log_1_minus_prob;  // log(1.0 - prob)
        bool           seed_set;

    public:
        //RandGeometric(const double p = (1.0 - 0.99717), const int seed = 1)
        RandGeometric(const double p = (1.0 - 0.9), const int seed = -100)
            : prob(p), seed_set(false)
        {
            init(prob, seed);
        };

        void           init(const double p, const int seed = 1)
        {
            seed_set = false;
            prob = p;
            if (prob != 1.0) { log_1_minus_prob = log(1.0 - prob); }
            unif.init(seed);
            seed_set = true;
        };

        long           draw()
        {
            assert(seed_set == true);
            if (prob == 1.0) { return(0); }
            long k = static_cast<long>(log(unif.draw()) / log_1_minus_prob);
            return(k);
        };
};

#endif // __RANDGEOMETRIC_H__


