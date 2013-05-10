#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <iostream>
#include <iomanip>
#include "VectorUtility.h"

template<class T_VALUE, class T_COUNT>
class Histogram {
    private:
        std::map<T_VALUE, T_COUNT> Hist;
        typedef typename std::map<T_VALUE, T_COUNT>::const_iterator HistCI;
        std::string value_name;
        std::string count_name;
        std::string freq_name;

    public:
        Histogram(const std::vector<T_VALUE>& Vec,
                  bool drop_zero = true,
                  T_VALUE min_val = T_VALUE(),
                  T_VALUE max_val = T_VALUE(),
                  bool use_min = false,
                  bool use_max = false) 
        {
            if (Vec.size() == 0) {
                std::cerr << "Histogram<T_VALUE,T_COUNT> : no values" 
                    << std::endl;
                exit(1);
            }
            names();
            fill(Vec, drop_zero, min_val, max_val, use_min, use_max);
        };

        void names(const std::string& nv = "", const std::string& nc = "",
                   const std::string& nf = "")
        { value_name = nv; count_name = nc; freq_name = nf; };

        long int size() { return(Hist.size()); }

        void fill(const std::vector<T_VALUE>& Vec,
                  bool drop_zero = true,
                  T_VALUE min_value = T_VALUE(),
                  T_VALUE max_value = T_VALUE(),
                  bool use_min = false,
                  bool use_max = false) 
        {
            if (! drop_zero) {
                // initialize all values in range to (T_COUNT)0
                // one way to extend this would be bringing it down to other
                // values so that the whole range is represented
                T_VALUE min, max;
                if (min_value != T_VALUE() || use_min) { min = min_value; }
                else { min = VectorUtility::Min(Vec); }
                if (max_value != T_VALUE() || use_max) { max = max_value; }
                else {max = VectorUtility::Max(Vec); }
                for (T_VALUE value = min; value <= max; ++value) {
                    Hist[value] = static_cast<T_COUNT>(0);
                }
            }
            for (long i = 0; i < Vec.size(); ++i) { Hist[Vec[i]]++; }
        };

        void print(std::ostream& os = std::cout,
                   const std::string& prefix = "") const
        { print_table(os, true, prefix); }

        void print_table(std::ostream& os = std::cout,
                         const bool header = true,
                         const std::string& prefix = "") const
        {
            if (header) {
                os << prefix;
                os << value_name;
                os << "\t" << count_name;
                os << "\t" << freq_name;
                os << std::endl;
            }
            T_COUNT count_sum = static_cast<T_COUNT>(0);
            for (HistCI p = Hist.begin(); p != Hist.end(); ++p) {
                count_sum += p->second; 
            }
            for (HistCI p = Hist.begin(); p != Hist.end(); ++p) {
                os << prefix;
                os << p->first;
                os << "\t" << p->second;
                os << "\t" << std::setprecision(4) 
                    << ((double)p->second)/((double)count_sum);
                os << std::endl;
            }
        };

        const std::vector<T_VALUE> values() const
        {
            std::vector<T_VALUE> ans;
            for (HistCI p = Hist.begin(); p != Hist.end(); ++p) {
                ans.push_back(p->first); 
            }
            return(ans);
        };

        const std::vector<T_COUNT> counts() const
        {
            std::vector<T_COUNT> ans;
            for (HistCI p = Hist.begin(); p != Hist.end(); ++p) {
                ans.push_back(p->second); 
            }
            return(ans);
        };

        friend std::ostream& operator<<(std::ostream& os, const Histogram& h) 
        { os << "Histogram:  "; h.print(os); os << std::endl; };
};

#endif // __HISTOGRAM_H__

