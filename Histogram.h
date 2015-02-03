#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <iostream>
#include <iomanip>
#include "VectorUtility.h"

/*! @class Histogram
    @brief Template to create a histogram of counts of unique values in a vector.
 */
template<class T_VALUE, class T_COUNT>
class Histogram {
    private:
        std::map<T_VALUE, T_COUNT> Hist;
        typedef typename std::map<T_VALUE, T_COUNT>::const_iterator HistCI;
        std::string value_name;
        std::string count_name;
        std::string freq_name;

    public:
        /*! Constructor

            @param Vec        vector of values, type used for template
            @param drop_zero  bool, whether to drop zero-valued values (true)
            @param min_val    minimum_value of histogram range
            @param max_val    maximum_value of histogram range
            @param use_min    bool, whether to use the minimum value (false)
            @param use_max    bool, whether to use the maximum value (false)
         */
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

        /*! Create names for value, counts, and frequencies.

            @param nv  string, name of value
            @param nc  string, name of count
            @param nf  string, name of frequency
         */
        void names(const std::string& nv = "", const std::string& nc = "",
                   const std::string& nf = "")
        { value_name = nv; count_name = nc; freq_name = nf; };

        //! Return histogram size
        size_t size() { return(Hist.size()); }

        /*! Fill the histogram, workhorse method

            If minimum and maximum values are not provided, use the min and
            max from Vec unless use_min or use_max are true, respectively, then
            use the parameter value.

            @param Vec        vector of values, type used for template
            @param drop_zero  bool, whether to drop zero-valued values (true)
            @param min_val    minimum_value of histogram range
            @param max_val    maximum_value of histogram range
            @param use_min    bool, whether to use the minimum value (false)
            @param use_max    bool, whether to use the maximum value (false)
         */
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
            for (size_t i = 0; i < Vec.size(); ++i) { Hist[Vec[i]]++; }
        };

        //! Simple print method
        void print(std::ostream& os = std::cout,
                   const std::string& prefix = "") const
        { print_table(os, true, prefix); }

        //! Print class contents as a table
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

        //! Return a vector of observed values
        const std::vector<T_VALUE> values() const
        {
            std::vector<T_VALUE> ans;
            for (HistCI p = Hist.begin(); p != Hist.end(); ++p) {
                ans.push_back(p->first); 
            }
            return(ans);
        };

        //! Return a vector of observed counts
        const std::vector<T_COUNT> counts() const
        {
            std::vector<T_COUNT> ans;
            for (HistCI p = Hist.begin(); p != Hist.end(); ++p) {
                ans.push_back(p->second); 
            }
            return(ans);
        };

        //! Suitable for inclusion in an ostream oparation
        friend std::ostream& operator<<(std::ostream& os, const Histogram& h) 
        { 
            os << "Histogram:  "; 
            h.print(os); 
            os << std::endl; 
            return(os);
        };
};

#endif // __HISTOGRAM_H__

