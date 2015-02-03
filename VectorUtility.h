#ifndef __VECTORUTILITY_H__
#define __VECTORUTILITY_H__

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <iostream>

namespace VectorUtility {

    typedef long Scalar;

    template<class T> inline const T        Max(const T& a, const T& b);
    template<class T> inline const T        Max(const std::vector<T>& Vec);
    template<class T> inline const T        Abs(const T a);
    template<class T> inline const T        Min(const T& a, const T& b);
    template<class T> inline const T        Min(const std::vector<T>& Vec);
    template<class T> inline const double   Sum(const std::vector<T>& Vec);
    template<class T> inline const double   SumSquares(const std::vector<T>& Vec);
    template<class T> inline const double   Mean(const std::vector<T>& Vec);
    template<class T> inline const double   Var(const std::vector<T>& Vec, 
                                                bool sample = true);
    template<class T> inline const double   Var2(const std::vector<T>& Vec);
    template<class T> inline const bool     IsIn(const std::vector<T>& Vec, 
                                                const T& val);
    template<class T> inline const std::vector<T>
                                            Seq(const T& from, const T& to);
    template<class T> inline const std::vector<T>  
                                            Uniq(const std::vector<T>& Vec);
    template<class T> inline std::map<T, std::vector<Scalar> >
                                            Runs(const std::vector<T>& Vec);

    //! \brief Pairwise maximum for class T.
    //!
    //! \param  a       scalar of class T
    //! \param  b       scalar of class T
    //! \return         maximum of a and b, class T.
    //!
    template<class T>
    inline const T
    Max(const T& a, const T& b)
    {
        return((a) > (b) ? (a) : (b));
    }

    //! \brief Pairwise minimum for class T.
    //!
    //! \param  a       scalar of class T
    //! \param  b       scalar of class T
    //! \return         minimum of a and b, class T.
    //!
    template<class T>
    inline const T
    Min(const T& a, const T& b)
    {
        return((a) < (b) ? (a) : (b));
    }

    //! \brief Mean of values in vector of class T, returned as a double.
    //!
    //! \param  Vec     vector of class T.
    //! \return         mean of Vec, double.
    //!
    template<class T>
    inline const double
    Mean(const std::vector<T>& Vec)
    {
        assert(Vec.size() > 0);
        return(Sum(Vec) / static_cast<double>(Vec.size()));
    }

    //! \brief Sum of values in vector of class T, returned as a double.
    //!
    //! \param  Vec     vector of class T.
    //! \return         sum of Vec, double.
    //!
    template<class T>
    inline const double
    Sum(const std::vector<T>& Vec)
    {
        assert(Vec.size() > 0);
        double ans = static_cast<double>(Vec[0]);
        for (long i = 1; i < Vec.size(); ++i) 
            { ans += static_cast<double>(Vec[i]); }
        return(ans);
    }

    //! \brief Sum of squared values in vector of class T, returned as a 
    //!        double.
    //!
    //! \param  Vec     vector of class T.
    //! \return         sum of squared values in Vec, double.
    //!
    template<class T>
    inline const double
    SumSquares(const std::vector<T>& Vec)
    {
        assert(Vec.size() > 0);
        double ans = (static_cast<double>(Vec[0]) * static_cast<double>(Vec[0]));
        for (long i = 1; i < Vec.size(); ++i) { 
            ans += (static_cast<double>(Vec[i]) * static_cast<double>(Vec[i]));
        }
        return(ans);
    }

    //! \brief Variance of values in vector of class T, returned as a double.
    //!
    //! \param  Vec     vector of class T.
    //! \param  sample  bool, calculates sample variance (SS/(n-1)) if true,
    //!                 and population variance (SS/n) if false.
    //! \return         variance of values in Vec, double.
    //!
    template<class T>
    inline const double
    Var(const std::vector<T>& Vec, bool sample = true)
    {
        long N = Vec.size();
        assert(N > 0);
        double sumsq = SumSquares(Vec);
        double sum = Sum(Vec);
        double ans = sumsq - ((sum * sum)/N);
        ans /= (N - (sample ? 1 : 0));
        return(ans);
    }

    //! \brief Maximum value in vector of class T.
    //!
    //! \param  Vec     vector of class T.
    //! \return         maximum in vector, class T.
    //!
    template<class T>
    inline const T
    Max(const std::vector<T>& Vec)
    {
        T ans;
        if (Vec.size() == 0) return(static_cast<T>(-9999999));
        if (Vec.size() > 1) { ans = Max(Vec[0], Vec[1]); } else { return(Vec[0]); }
        for (size_t i = 2; i < Vec.size(); ++i) { ans = Max(ans, Vec[i]); }
        return(ans);
    }

    //! \brief Minimum value in vector of class T.
    //!
    //! \param  Vec     vector of class T.
    //! \return         minimum in vector, class T.
    //!
    template<class T>
    inline const T
    Min(const std::vector<T>& Vec) 
    {
        T ans;
        if (Vec.size() == 0) return(static_cast<T>(9999999));
        if (Vec.size() > 1) { ans = Min(Vec[0], Vec[1]); } else { return(Vec[0]); }
        for (size_t i = 2; i < Vec.size(); ++i) { ans = Min(ans, Vec[i]); }
        return(ans);
    }

    //! \brief Pairwise absolute value for class T.
    //!
    //! \param  a       scalar of class T.
    //! \return         absolute value of a, class T.
    //!
    template<class T>
    inline const T
    Abs(const T a)
    {
        return((a) > ((T)0) ? (a) : -(a));
    }

    //! \brief Determine if value of class T is in vector of class T.
    //!
    //! \param  Vec     vector of class T.
    //! \param  val     value of class T.
    //! \return         bool, true of val is in Vec
    //!
    template<class T>
    inline const bool
    IsIn(const std::vector<T>& Vec, const T& val)
    {
        assert(Vec.size() > 0);
        for (size_t i = 0; i < Vec.size(); ++i) { if (val == Vec[i]) return(true); }
        return(false);
    }

    //! \brief Construct sequence of class T between two values.
    //!
    //! \param  from    start of sequence, class T.
    //! \param  to      end of sequence, class T.
    //! \return         vector of class T, sequence from .. to
    //!
    template<class T>
    inline const std::vector<T>
    Seq(const T& from, const T& to)
    {
        T by = (from > to) ? static_cast<T>(-1) : static_cast<T>(1);
        size_t N = ABS(to - from) + 1;
        std::vector<T> ans(N);
        T val = from;
        for (size_t i = 0; i < N; val += by, ++i) {
            ans[i] = val;
        }
        return(ans);
    }

    //! \brief Construct vector containing only unique values from vector 
    //!        of class T.
    //!
    //! \param  Vec     vector of class T.
    //! \return         vector of class T, unique values from Vec
    //!
    template<class T>
    inline const std::vector<T>
    Uniq(const std::vector<T>& Vec)
    {
        std::vector<T> ans; // vector to hold the unique values
        assert(Vec.size() > 0);
        ans.push_back(Vec[0]);
        for (size_t i = 1; i < Vec.size(); ++i) {
            if (! IsIn(ans, Vec[i])) { ans.push_back(Vec[i]); }
        }
        return(ans);
    }

    //! \brief Construct map of runs of unique values in a vector,
    //!        in order.
    //!
    //! A map of unique values in Vec, each mapping to a vector of
    //! Scalar (long) values containing the run lengths of that value,
    //! in order of appearance.  This is not a map of all runs in order,
    //! that has yet to be constructed.
    //!
    //! \param  Vec     vector of class T.
    //! \return         map<T, vector<Scalar> >, unique values from Vec
    //!                 mapping to a vector of run lengths.
    //!
    template<class T>
    inline std::map<T, std::vector<Scalar> >
    Runs(const std::vector<T>& Vec)
        // this RUNS() returns a map which maps each unique value to
        // a vector that contains its run lengths, in order
    {
        assert(Vec.size() > 0);
        std::map<T, std::vector<Scalar> > m;
        Scalar thisrun = 1;
        T prev = Vec[0];
        for (size_t i = 1; i < Vec.size(); ++i) {
            if (Vec[i] == prev) { ++thisrun; }
            else {
                m[prev].push_back(thisrun);
                thisrun = 1;
                prev = Vec[i];
            }
        }
        m[prev].push_back(thisrun);
        return(m); 

    }


}  // namespace VectorUtility


#endif // __VECTORUTILITY_H__

