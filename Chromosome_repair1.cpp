#include "Chromosome.h"

/*
 * Chromosome::repair1
 *
 * This implements a naive policy for double-stranded break repair, which draws
 * a direction of repair from uniform and then draws a conversion tract length
 * from an empirically-determined distribution.  We (will be able to) handle
 * several different tract length distributions.  The one we currently
 * implement is due to Hilliker et al. 1994 Meiotic gene conversion tract
 * length distribution within the rosy locus of Drosophila melanogaster.
 * Genetics 137:1019-1026.
 *
 * There are several different policies that we need to consider, even within
 * this most naive case.  First, what do we do when the conversion tract is
 * going to extend past the end of the chromosome?  We can (a) draw tract
 * lengths until one does not; (b) kill the chromosome; (c) just convert out to
 * the end and leave it at that, effectively truncating the tract to the length
 * of chromosome available for it.
 *
 */

void 
Chromosome::repair1 ( )
{
    _trace("repair1 ( )");

    const bool dbgthis = true;
    bool       truncated = true;

    if (DSBreakQueue.size() == 0) { return; }
    if (dbgthis) {
        //  Remove or otherwise header always printed
        //std::cout << "Chromosome::repair1()" << std::endl;
        //std::cout << "action\titer_event\t";
        //DSBreakEvent::print_header(std::cout);
    }
    // DSBreakDequeI p;
    // DSBreakDequeCI cp;
    long count = 1;
    while (DSBreakQueue.size()) {
        if (count > 1) { 
            std::cerr << "Chromosome::repair1 : more than one DSB, not implemented"
                << std::endl;
            exit(1);
        }
        const DSBreakEvent& event = DSBreakQueue.back();
        if (dbgthis) {
            print_centered(std::cout, event.event_site);
        }
        if (event.event_site >= min_DSB_site) {
            // it is a valid DSB (placeholder)
            
            double dirdraw = repair1_Uniform.draw();
            int dir = (dirdraw < 0.5) ? (-1) : (1);
            SeqSize tract_length = repair1_Geometric.draw();
            if (tract_length > 0) {
                SeqSize tract_end = event.event_site + (tract_length * dir);
                // truncate the end of the tract to the end of the chromosome
                if (tract_end < 0) 
                    { tract_end = 0; truncated = true; }
                else if (tract_end >= get_nbp()) { 
                    tract_end = get_nbp() - 1; 
                    truncated = true; 
                }
                for (SeqSize i = event.event_site; i != tract_end; i += dir) {
                    X[i] = HOMZ;
                }
                X[tract_end] = HOMZ;
                if (dbgthis) {
                    std::cout << "Chromosome::repair1 : site = " 
                        << event.event_site 
                        << "  dir = " << dir 
                        << "  tract_length = " << tract_length
                        << "  tract_end = " << tract_end << std::endl;
                    print_centered(std::cout, event.event_site);
                    std::cout << std::endl;
                }
            }
            DSBreakQueue.pop_back();
        }
        ++count;
    }
};

