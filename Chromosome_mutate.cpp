#include "Chromosome.h"

void 
Chromosome::mutate()
{
    _trace("mutate ( )");

    // Chromosome::mutate() for the 0/1 homozygous/heterozygous chromosome model
    //
    // 1  Determine if an event occurred, via a draw from a uniform distribution
    //    that is less than (rate per site * number of sites).  We might 
    //    eventually have to expand this to >1 events per tick if number of sites
    //    gets big enough.
    // 2  Determine the site at which it occurred, via a draw from a uniform
    //    distribution multiplied by number of sites then rounded down.
    // 3a If the affected site was homozygous, make it heterozygous, and DONE.
    // 3b If the affected site was heterozygous, it has a 1/3 chance of turning
    //    homozygous ...
    // 4  So to save computing, check to see if our first uniform draw was less
    //    than 1/3 of (rate per site * number of sites).  If so, it's a
    //    homozygous-making event, make site homozygous, and DONE.
    // 5  Otherwise, leave site heterozygous, and DONE.
    //

    double mut_event_threshold = (get_mu() * get_nbp());  // site rate * num sites
    double event_draw;
    const double homz_fraction = (1.0/3.0);
    if ((event_draw = mutate_Uniform.draw()) < mut_event_threshold) {
        SeqSize mutsite = static_cast<SeqSize>(mutate_Uniform.draw() * get_nbp());
        bp site_old = X[mutsite];
        if (is_homozygous(X[mutsite])) { 
            X[mutsite] = HETZ; 
        } else {
            if (event_draw < (mut_event_threshold * homz_fraction)) {
                X[mutsite] = HOMZ;
            }
        }
        MutationEvent event;
        event.event = MutationLog.size();
        event.event_threshold = mut_event_threshold;
        event.event_draw = event_draw;
        event.event_site = mutsite;
        event.val_orig = site_old;
        event.val_new = X[mutsite];
        MutationLog.push_back(event);
        _did_mutate= true;
    } else { 
        _did_mutate = false; 
    }
};

