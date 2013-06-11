#include "GC.h"
#include "Chromosome.h"
#include "SequenceRuns.h"

typedef long Scalar;

int main () {
    Chromosome C(1000);
    C.set_mu(0.0000001);
    C.set_c(0.000001);
    //C._debug_trace = true;
    C.set_heterozygosity(0.4);
    //std::cout << C;
    C.print_stats();
    long num_events = 40;
    while (num_events > 0) {
        C.mutate();
        C.dsbreak();
        if (C.get_did_break()) { --num_events; }
        C.repair1();
    }
    std::cout << std::endl << "num mutations = " << C.number_mutations()
        << "  num dsbreaks = " << C.number_dsbreaks() << std::endl;
//    for (long i = 0; i < 1000000000; ++i) {
//        C.mutate();
//    }
    //C.print_stats();
    //C.print_mutations();
    //C.print_dsbreaks();


//    SequenceRuns<Chromosome::bp, Scalar> SR( C.get_sequence() );
//    SR.print_summary_stats();
//    SR.print_histograms();
//    SR.print_table();

}
