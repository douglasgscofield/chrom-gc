#include "Chromosome.h"

void 
Chromosome::dsbreak()
{
    _trace("dsbreak ( )");

    // Chromosome::dsbreak() for the 0/1 homozygous/heterozygous chromosome model
    //
    // 1  Determine if an event occurred, via a draw from a uniform distribution
    //    that is less than (break rate per site * number of sites).  We might 
    //    eventually have to expand this to >1 events per tick if number of sites
    //    and/or break rate.
    // 2  Determine the site at which it occurred, via a draw from a uniform
    //    distribution multiplied by number of sites then rounded down.  Note
    //    that the number of break sites (since breaks only occur between bases)
    //    is the length of the chromosome - 1.  We might want to handle telomeres
    //    here.
    // 3  Add an entry to the DSBreakEventLog
    // 4  Save the decision to repair until Chromosome::mmr(), which fetches
    //    the DSBreakEventLog entry and does its thing.
    //

    SeqSize num_sites = get_nbp() - min_DSB_site; // number of potential break sites
    double break_event_threshold = (get_c() * num_sites);  // rate * num sites
    double event_draw;
    if ((event_draw = dsbreak_Uniform.draw()) < break_event_threshold) {
        SeqSize breaksite = 
            static_cast<SeqSize>((dsbreak_Uniform.draw() * num_sites)) 
            + min_DSB_site;
        // we only create the entry, Chromosome::mmr() fixes it
        DSBreakEvent event;
        event.event = DSBreakLog.size();
        event.event_threshold = break_event_threshold;
        event.event_draw = event_draw;
        event.event_site = breaksite;
        DSBreakLog.push_back(event);  // add to the global log
        DSBreakQueue.push_back(event);  // add to the (this-iteration) queue
        _did_break= true;
    } else { 
        _did_break = false; 
    }
};

