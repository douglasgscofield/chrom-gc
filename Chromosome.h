#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#define RANDOM_SEED_FLAG false

#include "GC.h"
#include "VectorUtility.h"
#include "RandUniform.h"
#include "RandBinomial.h"
#include "RandGeometric.h"
#include "Histogram.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <deque>
#include <list>
#include <map>

class Chromosome {

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// The base-pair (bp) model and the sequence, interfaces and members
//
// // // // // // // // // // // // // // // // // // // // // // // //

    public:

        // a better choice is typedef sequence_type::size_type SeqSize;
        typedef long SeqSize; 

    private:

        SeqSize               _nbp; // nbp: number of bp to model

        // for set_heterozygosity()
        bool                  _random_seed;
        RandUniform           Unif;

    public:

        typedef short         bp;  // begin expanding our concept of bp here
        enum                  bp_state { HOMZ = 0, HETZ = 1 };
        bool                  is_homozygous(bp state) { return(state == HOMZ); }
        bool                  is_heterozygous(bp state) { return(state == HETZ); }
        typedef std::vector<bp>  sequence_type;

        sequence_type         X;

        void                  set_nbp(SeqSize n) { _nbp = n; };
        SeqSize               get_nbp() const { return(_nbp); };
        const sequence_type&  get_sequence() { return X; };
        SeqSize               size() const { check(); return(get_nbp()); };
        void                  fill(bp bpstate) { X.assign(get_nbp(), bpstate); };

        void                  init(SeqSize nnbp = -1) 
        {
            if (nnbp >= 0) { set_nbp(nnbp); X.resize(get_nbp()); }
            Unif.init(_random_seed);
            fill(HOMZ);
        };

        void                  set_heterozygosity(double het = 0.0) 
        {
            fill(HOMZ);
            for (SeqSize i = 0; i < X.size(); ++i)
                { if (Unif.draw() < het) X[i] = HETZ; }
        };

    private:

        void                  check() const 
        {
            if (get_nbp() != X.size()) {
                std::cerr << "Chromosome::check() : _nbp changed without init()"
                    << std::endl;
            }
            assert(get_nbp() == X.size());
        };

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// Constructors and destructors
//
// // // // // // // // // // // // // // // // // // // // // // // //

    public:

        Chromosome(const int sn = 0) 
            : _mu(0.0), _c(0.0), _did_mutate(false), _did_break(false),
              _random_seed(RANDOM_SEED_FLAG), _debug_trace(false)
        { _trace("CONSTRUCTOR ( sn )"); init(sn); };

        ~Chromosome() { /* empty */ };

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// Mutation-related interfaces and members
//
// // // // // // // // // // // // // // // // // // // // // // // //

    private:

        double            _mu;  // mutation rate per bp
        bool              _did_mutate;
        RandUniform       mutate_Uniform;

        // Keeping track of mutation events
        struct MutationEvent {
            long     event;
            double   event_threshold;
            double   event_draw;
            SeqSize  event_site;
            bp       val_orig;
            bp       val_new;

            static void   print_header(std::ostream& os = std::cout) 
            {
                os << "event\tevent_threshold\tevent_draw\tevent_site" 
                    << "\tval_orig\tval_new" << std::endl;
            };

            void          print(std::ostream& os = std::cout) const 
            {
                os << event << "\t" << event_threshold << "\t" 
                    << event_draw << "\t" << event_site << "\t" 
                    << val_orig <<"\t" << val_new << std::endl;
            };
        };

        std::deque<MutationEvent>    MutationLog;

        typedef std::deque<MutationEvent>::iterator         MutationDequeI;
        typedef std::deque<MutationEvent>::const_iterator   MutationDequeCI;

    public:

        void              mutate();

        double            get_mu() const { return(_mu); };
        void              set_mu(double m) { _mu = m; };
        bool              get_did_mutate() const { return(_did_mutate); };
        long              number_mutations() const { return(MutationLog.size()); };

        void              print_mutations(std::ostream& os = std::cout,
                                          bool header = true) const
        {
            if (header) { MutationEvent::print_header(os); }
            MutationDequeCI p;
            for (p = MutationLog.begin(); p != MutationLog.end(); ++p)
                { (*p).print(os); }
        };

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// Double-stranded break related interfaces and members
//
// // // // // // // // // // // // // // // // // // // // // // // //

    private:

        double             _c;   // gene conversion rate
        bool               _did_break;
        RandUniform        dsbreak_Uniform;

        // Keep track of double-stranded breaks.  We use the same
        // structure for recording two functionally different deques
        // of DSBs:
        //   - DSBreakLog, a log, like for mutations, that tracks
        //     all DSB events
        //   - DSBreakQueue, a queue that holds DSBs that must be
        //     repaired; this is filled and emptied as
        //     DSBs are created and repaired.  
        struct DSBreakEvent { 
            long     event;
            double   event_threshold;
            double   event_draw;
            SeqSize  event_site;

            static void print_header(std::ostream& os = std::cout) 
            {
                os << "event_threshold\tevent_draw\tevent_site" << std::endl;
            };

            void print(std::ostream& os = std::cout) const 
            {
                os << event_threshold << "\t" << event_draw << "\t"
                    << event_site << std::endl;
            };
        };

        std::deque<DSBreakEvent>   DSBreakLog;
        std::deque<DSBreakEvent>   DSBreakQueue;

        typedef std::deque<DSBreakEvent>::iterator         DSBreakDequeI;
        typedef std::deque<DSBreakEvent>::const_iterator   DSBreakDequeCI;

        enum { min_DSB_site = 1 };  // a named constant; we can't break beyond here

    public:

        void              dsbreak();

        double            get_c() const { return(_c); };
        void              set_c(double c) { _c = c; };
        bool              get_did_break() const { return(_did_break); };
        long              number_dsbreaks() const { return(DSBreakLog.size()); };

        void              print_dsbreaks(std::ostream& os = std::cout,
                                         bool header = true) const
        { 
            if (header) { DSBreakEvent::print_header(os); }
            DSBreakDequeCI p;
            for (p = DSBreakLog.begin(); p != DSBreakLog.end(); ++p)
                { (*p).print(os); }
        };

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// Repair-related interfaces and members
//
// // // // // // // // // // // // // // // // // // // // // // // //

    private:

        RandUniform       repair1_Uniform;
        RandGeometric     repair1_Geometric;
        //RandGeometric     repair1_Geometric((1.0 - 0.99717), -100);

    public:

        // repair0() implements a minimal repair that simply removes the
        // record of the break, effectively repairing the break with
        // tract length 0 and no errors.

        void              repair0();

        // repair1() implements a naive repair based on simply drawing a
        // direction from uniform and a tract length from a tract-length
        // distribution.  It will be extended to model a variety of
        // tract-length distributions.

        void              repair1();


// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// Other general interfaces and members
//
// // // // // // // // // // // // // // // // // // // // // // // //

    private:

        bool                  _debug_trace;

        void                  _trace(const std::string& s) const 
        { if (_debug_trace) { std::cerr << "Chromosome::" << s << std::endl; } };

    public:

        bool                  get_debug_trace() const { return(_debug_trace); };
        void                  set_debug_trace(bool dt) { _debug_trace = dt; };
        void                  print_stats(std::ostream& os = std::cout, 
                                          bool header = true) const;
        void                  print(std::ostream& os = std::cout, 
                                    SeqSize startbp = 0, 
                                    SeqSize endbp = (-1), 
                                    SeqSize width = 75, 
                                    const SeqSize markbp = (-1), 
                                    const bool header = true) const;
        void                  print_centered(std::ostream& os, 
                                    const SeqSize markbp, 
                                    const SeqSize stride = 30, 
                                    const SeqSize width = 75) const;

        friend std::ostream& operator<<(std::ostream& os, const Chromosome& c) 
        {
            c._trace(" friend operator<< ( os, c )");
            os << "Chromosome:  "; c.print(os); os << std::endl;
        };
};

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
//
// Chromosome:: methods defined outside the class declaration
//
// Also, some definition in separate files named Chromosome_<method>.cpp
//
// // // // // // // // // // // // // // // // // // // // // // // //

inline void 
Chromosome::print_stats(std::ostream& os, bool header) const
{
    _trace("print_stats ( os, header )");

    if (X.size() == 0) { os << "zero-length chromosome" << std::endl; return; }

    // site counts
    if (header) { 
        os << "Chromosome:: Summary Statistics" << std::endl; 
        os << "===============================" << std::endl;
    }
    Histogram<bp, Scalar> site_histogram ( X, false);
    site_histogram.names("bp_state", "num_sites", "freq_sites");
    site_histogram.print_table(os, header);
};


inline void
Chromosome::print(std::ostream& os, SeqSize startbp, SeqSize endbp, 
    SeqSize width, const SeqSize markbp, const bool header) const
{
    _trace("print ( os, startbp, endbp, width )");

    const bool note_markbp = true;  // to mark the line in which markbp is located
    const bool append_markbp = true;  // to append the markbp to its line
    const bool center_markbp = true; // to center the line on the markbp
    const std::string pad = " ";
    const SeqSize size = get_nbp();

    if (startbp < 0 || startbp > size - 1) { startbp = 0; }
    if (endbp < 0 || endbp > size - 1) { endbp = size - 1; }
    if (width < 0) { width = 70; }
    if (header) {
        os << "_nbp=" << get_nbp();
        os << "  _mu=" << get_mu();
        os << "  _c=" << get_c();
        os << std::endl;
    }
    for (SeqSize i = startbp; i <= endbp; i += width) {
        SeqSize endslice = VectorUtility::Min(i + width - 1, endbp);
        SeqSize left_pad = 0, right_pad = 0;
        if (center_markbp && markbp >= i && markbp <= endslice) {
            SeqSize left_width = markbp - i; 
            SeqSize right_width = endslice - markbp; 
            if (left_width < right_width)
                { left_pad = right_width - left_width; right_pad = 0; }
            else if (left_width > right_width) 
                { left_pad = 0; right_pad = left_width - right_width; }
        }
        if (note_markbp && markbp >= i && markbp <= endslice)
            { os << std::setw(7) << i << "* "; }
        else
            { os << std::setw(7) << i << "  "; }
        while (left_pad > 0) { os << pad; --left_pad; }
        for (SeqSize j = i; j <= endslice; ++j) {
            if (j == markbp) { os << "|"; }
            os << ((X[j] == true) ? "1" : "0");
            //if (j == markbp) { os << ">"; }
        }
        while (right_pad > 0) { os << pad; --right_pad; }
        if (append_markbp && markbp >= i && markbp <= endslice)
            { os << ":" << std::setw(7) << markbp; }
        os << std::endl;
    }
};


inline void
Chromosome::print_centered(std::ostream& os, const SeqSize markbp, 
                           const SeqSize stride, const SeqSize width) const {
    _trace("print_centered ( os, markbp, stride, width )");
    // stride notes the number of bases on either side of markbp
    // to print, while width has its meaning as in print()
    SeqSize startbp = markbp - stride;
    SeqSize endbp = markbp + stride;
    print(os, startbp, endbp, width, markbp, false);
};

// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // //

#endif // CHROMOSOME_H

