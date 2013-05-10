#ifndef SEQUENCERUNS_H
#define SEQUENCERUNS_H

#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cassert>
#include <typeinfo>
#include <iostream>
#include <sstream>
#include "VectorUtility.h"
#include "Histogram.h"

//template<class T_ITEM>
//class Runs : public class InternalRuns<T_ITEM, Scalar>;
//};
//
//template<class T_ITEM, class T_LENGTH>
//class InternalRuns {
//};

template<class T_ITEM, class T_COUNT>
class SequenceRuns {
    public:
        bool _debug_trace;

        class Run {
            public:
                long     run_index;
                T_COUNT  position;
                T_COUNT  length;
                T_ITEM   item;
                void print(std::ostream& os = std::cout) const
                { os << "(" << run_index << " @" << position << " " << item 
                    << ": " << length << " )"; };
                void print_table(std::ostream& os = std::cout) const
                { os << run_index << "\t" << position << "\t" << length 
                    << "\t" << item << std::endl; };
                friend std::ostream& operator<<(std::ostream& os, const Run& R)
                { R.print(os); };
        };

        typedef typename std::vector<T_ITEM>                     item_list_type;
        typedef typename std::vector<T_COUNT>                    run_list_type;
        typedef typename std::map<T_ITEM, std::vector<T_COUNT> > map_type;
        typedef typename map_type::const_iterator                map_type_CI;
        typedef typename std::map<T_ITEM, T_COUNT>               unique_item_type;
        typedef typename unique_item_type::const_iterator        unique_item_type_CI;
        typedef typename std::vector<Run>                        Runs_type;

        SequenceRuns(const std::vector<T_ITEM>& Vec)
            : _debug_trace(false), num_runs(0), total_items(0)
        {
            _trace("CONSTRUCTOR ( Vec )");
            if (Vec.size() == 0) {
                std::cerr << "SequenceRuns<>::CONSTRUCTOR : no values" << std::endl;
                exit(1);
            }
            names("item", "run_length");
            fill(Vec);
        };

    private:
        unique_item_type  unique_items;  // entry for each T_ITEM value
        long              num_runs;
        long              total_items;  // total number of items seen
        Runs_type         Runs;  // vector to which RunsList is converted
        map_type          Map;  // map of runs
        std::string       item_name;
        std::string       run_name;

        void _trace(const std::string& s) const {
            if (_debug_trace) {
                std::cerr << "SequenceRuns<" << typeid(T_ITEM).name() 
                    << "," << typeid(T_COUNT).name()
                    << ">::" << s << std::endl;
            }
        };

        void _fill(const std::vector<T_ITEM>& Vec, const long start_run_index);

    public:
        void names(const std::string& in = "", const std::string& rn = "")
        { _trace("names ( in, rn )"); item_name = in; run_name = rn; };

        void fill(const std::vector<T_ITEM>& Vec)
        {
            _trace("fill ( Vec )");
            Runs.clear();
            Runs.resize(Vec.size());
            _fill(Vec, 0);
            total_items = Vec.size();
        };

        void append(const std::vector<T_ITEM>& Vec)
        {
            _trace("append ( Vec )");
            if (Runs.size() < num_runs + Vec.size())
                { Runs.resize(num_runs + Vec.size()); }
            _fill(Vec, num_runs);
            total_items += Vec.size();
        };

        // Produce statistics on runs
        //
        void print_summary_stats(std::ostream& os = std::cout, 
                         const bool header = true,
                         const std::string& prefix = "") const
        {
            _trace("print_summary_stats ( )");
            if (header) {
                os << "SequenceRun:: Summary Statistics" << std::endl;
                os << "================================" << std::endl;
                os << prefix;
                os << "item_val";
                os << "\t" << "num_sites";
                os << "\t" << "freq";
                os << "\t" << "min_run";
                os << "\t" << "max_run";
                os << "\t" << "mean_run";
                os << "\t" << "var_run";
                os << std::endl;
            }
            for (map_type_CI p = Map.begin(); p != Map.end(); ++p) {
                if (p->second.size() == 0) {
                    os << prefix;
                    os << p->first << "\tNA\tNA\tNA\tNA\tNA\tNA" << std::endl;
                }
                // for item value p->first, produce statistics from
                // the vector of run lengths in p->second
                double sum = VectorUtility::Sum(p->second);
                T_ITEM min = VectorUtility::Min(p->second);
                T_ITEM max = VectorUtility::Max(p->second);
                double mean = VectorUtility::Mean(p->second);
                double var = VectorUtility::Var(p->second);
                os << prefix;
                os << p->first;
                os << "\t" << sum;
                os << "\t" << (sum / total_items);
                os << "\t" << min;
                os << "\t" << max;
                os << "\t" << mean;
                os << "\t" << var;
                os << std::endl;
            }
        };

        void print_histograms(std::ostream& os = std::cout, 
                             const bool header = true,
                             const std::string& prefix = "") const
        {
            _trace("print_histogram ( )");
            if (header) {
                os << "SequenceRun:: Runs Length Histogram" << std::endl;
                os << "===================================" << std::endl;
                os << prefix;
                os << "item_val";
                os << "\t" << "run_length";
                os << "\t" << "count";
                os << "\t" << "freq";
                os << std::endl;
            }
            for (map_type_CI p = Map.begin(); p != Map.end(); ++p) {
                if (p->second.size() == 0) {
                    os << prefix;
                    os << p->first << "\tNA\tNA\tNA" << std::endl;
                }
                Histogram<T_COUNT, T_COUNT> hist(p->second, false);
                std::ostringstream ost;
                ost << p->first << "\t";
                hist.print_table(os, false, ost.str());
            }
        };



        // Accessory data items //
        //
        void build_map();
        const map_type& get_map();
        const item_list_type get_item_list() const;
        const run_list_type get_run_list() const;

        // Print routines //
        //
        void print(std::ostream& os = std::cout) const;
        void print_table(std::ostream& os = std::cout,
                         const bool header = true) const;
        void print_runs(std::ostream& os = std::cout,
                        const int width = 3) const;
        void print_runs_table(std::ostream& os = std::cout, 
                              const bool header = true) const;
        void print_unique_items(std::ostream& os = std::cout,
                                const int width = 3) const;
        void print_unique_items_table(std::ostream& os = std::cout,
                                      const bool header = true) const;
        void print_map(std::ostream& os = std::cout) const;
        void print_map_table(std::ostream& os = std::cout,
                             const bool header = true) const;
        friend std::ostream& operator<<(std::ostream& os, const SequenceRuns& s)
        { 
            s._trace(" friend operator<< ( os, s )");
            s.print(os); os << std::endl; 
        };
};

template<class T_ITEM, class T_COUNT>
void
SequenceRuns<T_ITEM, T_COUNT>::_fill(const std::vector<T_ITEM>& Vec, 
                                     const long start_run_index)
{
    _trace("_fill ( Vec, start_run_index )");
    long run_index = start_run_index;
    T_ITEM prev = Vec[0];
    T_COUNT run_position = 0;
    T_COUNT thisrunlength = 1;
    for (long i = 1; i < Vec.size(); ++i) {
        if (Vec[i] == prev) { ++thisrunlength; }
        else {
            // note the run
            Runs[run_index].item = prev;
            Runs[run_index].length = thisrunlength;
            Runs[run_index].position = run_position;
            Runs[run_index].run_index = run_index;
            unique_items[prev]++;
            ++run_index;
            prev = Vec[i];
            run_position = i;
            thisrunlength = 1;
        }
    }
    Runs[run_index].item = prev;
    Runs[run_index].length = thisrunlength;
    Runs[run_index].position = run_position;
    Runs[run_index].run_index = run_index;
    unique_items[prev]++;
    num_runs = run_index + 1;
    build_map();
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::build_map()
{
    _trace("build_map ( )");
    Map.clear();
    for (long i = 0; i < num_runs; ++i) {
        Map[Runs[i].item].push_back( Runs[i].length );
    }
};

template<class T_ITEM, class T_COUNT>
const typename SequenceRuns<T_ITEM, T_COUNT>::map_type&
SequenceRuns<T_ITEM, T_COUNT>::get_map()
{
    _trace("get_map ( )");
    return(Map);
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_runs(std::ostream& os, const int width) const
{
    _trace("print_runs ( os, width )");
    os << "Runs: " << std::endl;
    int ww = 0;
    for (long i = 0; i < num_runs; ++i) {
        os << Runs[i];
        //os << prefix << i << midfix << "( " << Runs[i].item
        //    << " : " << Runs[i].length << " )" << infix;
        if (ww == (width - 1)) { os << std::endl; } else { os << "\t"; }
        ww = (ww + 1) % width;
    }
    if (ww != width) { os << std::endl; }
    //os << "__END" << std::endl;
    os << std::endl;
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_runs_table(std::ostream& os,
                                                const bool header) const
{
    _trace("print_runs_table ( os )");
    if (header) {
        os << "SequenceRuns:: Runs Table" << std::endl;
        os << "=========================" << std::endl;
        os << "run_index" << "\t" << "run_start" << "\t" << "run_length"
            << "\t" << "run_item" << std::endl;
    }
    for (long i = 0; i < num_runs; ++i) {
        Runs[i].print_table(os);
    }
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_unique_items(std::ostream& os,
                                                  const int width) const
{
    _trace("print_unique_items ( os, width )");
    os << "Unique items: ";
    int ww = 0;
    unique_item_type_CI p;
    for (p = unique_items.begin(); p != unique_items.end(); ++p) {
        os << p->first << ": " << p->second;
        if (ww == (width - 1)) { os << std::endl; } else { os << ", "; }
        ww = (ww + 1) % width;
    }
    if (ww != width) { os << std::endl; }
    //os << "__END" << std::endl;
    os << std::endl;
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_unique_items_table(std::ostream& os,
                                                        const bool header) const
{
    _trace("print_unique_items_table ( os )");
    if (header) {
        os << "SequenceRuns:: Unique Items Table" << std::endl;
        os << "=================================" << std::endl;
        os << "unique_item" << "\t" << "runs_count" << std::endl;
    }
    unique_item_type_CI p;
    for (p = unique_items.begin(); p != unique_items.end(); ++p) {
        os << p->first << "\t" << p->second << std::endl;
    }
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print(std::ostream& os) const
{
    _trace("print ( os )");
    print_runs(os, 3);
    print_unique_items(os, 3);
    os << "Map: " << std::endl;
    print_map(os); 
    //os << "__END" << std::endl;
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_table(std::ostream& os,
                                           const bool header) const
{
    _trace("print_table ( os )");
    print_runs_table(os, header);
    print_unique_items_table(os, header);
    print_map_table(os, header); 
};

template<class T_ITEM, class T_COUNT>
const typename SequenceRuns<T_ITEM, T_COUNT>::item_list_type
SequenceRuns<T_ITEM, T_COUNT>::get_item_list() const
{
    _trace("get_item_list ( )");
    item_list_type lst( num_runs );
    for (long i = 0; i < num_runs; ++i) { lst[i] = Runs[i].item; }
    return(lst);
};

template<class T_ITEM, class T_COUNT>
const typename SequenceRuns<T_ITEM, T_COUNT>::run_list_type
SequenceRuns<T_ITEM, T_COUNT>::get_run_list() const
{
    _trace("get_run_list ( )");
    run_list_type lst( num_runs );
    for (long i = 0; i < num_runs; ++i) { lst[i] = Runs[i].length; }
    return(lst);
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_map(std::ostream& os) const
{
    _trace("print_map ( os )");
    for (map_type_CI p = Map.begin(); p != Map.end(); ++p) {
        os << p->first << ": ";
        if (p->second.size() == 0) {
            os << "empty runs vector" << std::endl;
            continue;
        }
        os << p->second[0];
        for (long i = 1; i < p->second.size(); ++i) 
            { os << " " << p->second[i]; }
        os << std::endl;
    }
    //os << "__END" << std::endl;
    os << std::endl;
};

template<class T_ITEM, class T_COUNT>
void 
SequenceRuns<T_ITEM, T_COUNT>::print_map_table(std::ostream& os,
                                               const bool header) const
{
    _trace("print_map_table ( os )");
    if (header) {
        os << "SequenceRuns:: Runs Map Table" << std::endl;
        os << "=============================" << std::endl;
        os << "run_map_item" << "\t" << "run_map_index" << "\t"
            << "run_length" << std::endl;
    }
    for (map_type_CI p = Map.begin(); p != Map.end(); ++p) {
        if (p->second.size() == 0) {
            os << p->first << "\t" << "NA" << "\t" << "NA" << std::endl;
            continue;
        }
        for (long i = 0; i < p->second.size(); ++i) {
            os << p->first << "\t" << i << "\t" << p->second[i] 
                << std::endl;
        }
    }
};


#endif // SEQUENCERUNS_H

