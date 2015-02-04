#include "Chromosome.h"

/*! Method implementing a minimal policy for double-stranded break repair.

  @sa repair1

  This implements a minimal policy for double-stranded break repair, which
  simply removes the break record and acts as if the break was repaired in
  situ with no conversion (tract length 0).  This might be better folded into
  repair1() with tract length 0, but for now it's separated.
 
 */

void 
Chromosome::repair0 ( )
{
    _trace("repair0 ( )");

    const int debug = 1;

    //! We have double-stranded breaks.  We need to go through the
    //! log and process them en masse... There might be more than one and 
    //! they might interfere, but that's a sophisticated approach.  I'll
    //! have to implement different mmr policies, I think.

    if (DSBreakQueue.size() == 0) { return; }
    if (debug >= 1) {
        std::cout << "Chromosome::repair0 - there are " << DSBreakQueue.size()
            << " breaks to repair" << std::endl;
        if (debug >= 2) {
            std::cout << "Chromosome::repair0()" << std::endl;
            std::cout << "action\titer_event\t";
            DSBreakEvent::print_header(std::cout);
        }
    }
    // DSBreakDequeI p;
    // DSBreakDequeCI cp;
    long count = 1;
    while (DSBreakQueue.size()) {
        if (count > 1) { 
            std::cerr << "Chromosome::repair0 : more than one DSB, not implemented"
                << std::endl;
            exit(1);
        }
        const DSBreakEvent& event = DSBreakQueue.back();
        if (debug >= 2) {
            print_centered(std::cout, event.event_site);
        }
        if (debug >= 1) {
            std::cout << "repair" << "\t" << count << "\t";
            event.print(std::cout);
        }
        if (event.event_site >= min_DSB_site) {
            // it is a valid DSB (placeholder), just remove the record
            DSBreakQueue.pop_back();
        }
        ++count;
    }
};

