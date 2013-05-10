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
    template<class T> inline void           Trim(std::vector<T>& Vec, 
                                                 const T& val = ((T)0));
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

    template<class T>
    inline const T
    Max(const T& a, const T& b)
    {
        return((a) > (b) ? (a) : (b));
    }

    template<class T>
    inline const double
    Mean(const std::vector<T>& Vec)
    {
        assert(Vec.size() > 0);
        return(static_cast<double>(Sum(Vec))/static_cast<double>(Vec.size()));
    }

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

    template<class T>
    inline const T
    Max(const std::vector<T>& Vec)
    {
        T ans;
        if (Vec.size() == 0) return(static_cast<T>(-9999999));
        if (Vec.size() > 1) { ans = Max(Vec[0], Vec[1]); } else { return(Vec[0]); }
        for (long i = 2; i < Vec.size(); ++i) { ans = Max(ans, Vec[i]); }
        return(ans);
    }

    template<class T>
    inline const T
    Min(const T& a, const T& b)
    {
        return((a) < (b) ? (a) : (b));
    }

    template<class T>
    inline const T
    Min(const std::vector<T>& Vec) 
    {
        T ans;
        if (Vec.size() == 0) return(static_cast<T>(9999999));
        if (Vec.size() > 1) { ans = Min(Vec[0], Vec[1]); } else { return(Vec[0]); }
        for (long i = 2; i < Vec.size(); ++i) { ans = Min(ans, Vec[i]); }
        return(ans);
    }

    template<class T>
    inline const T
    Abs(const T a)
    {
        return((a) > ((T)0) ? (a) : -(a));
    }

    template<class T>
    inline void 
    Trim(std::vector<T>& Vec, const T& val = ((T)0)) 
    {
        long trim = 0;
        assert(Vec.size() > 0);
        for (long i = Vec.size() - 1; i >= 0; --i) {
            if (Vec[i] == val) { ++trim; } else { break; }
        }
        Vec.resize(Vec.size() - trim);
    }

    template<class T>
    inline const bool
    IsIn(const std::vector<T>& Vec, const T& val)
    {
        assert(Vec.size() > 0);
        for (long i = 0; i < Vec.size(); ++i) { if (val == Vec[i]) return(true); }
        return(false);
    }

    template<class T>
    inline const std::vector<T>
    Seq(const T& from, const T& to)
    {
        T by = (from > to) ? static_cast<T>(-1) : static_cast<T>(1);
        long N = ABS(to - from) + 1;
        std::vector<T> ans(N);
        T val = from;
        for (long i = 0; i < N; val += by, ++i) {
            ans[i] = val;
        }
        return(ans);
    }

    template<class T>
    inline const std::vector<T>
    Uniq(const std::vector<T>& Vec)
    {
        std::vector<T> ans; // vector to hold the unique values
        assert(Vec.size() > 0);
        ans.push_back(Vec[0]);
        for (long i = 1; i < Vec.size(); ++i) {
            if (! IsIn(ans, Vec[i])) { ans.push_back(Vec[i]); }
        }
        return(ans);
    }

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
        for (long i = 1; i < Vec.size(); ++i) {
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

