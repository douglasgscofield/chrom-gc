# chrom-gc

Simulate gene conversion on an arbitrary chromosome.  Quiet for a while and not
complete.  Eventually I would like to gain some insight into chromosome-wide
variation in GC content and other features as features affecting recombination
vary.

## Runs of homozygosity

The first target is simulating runs of homozygosity.  Mutation creates
heterozygosity, while gene conversion creates homozygous patches of varying
size.  What steady-state distributions of homozygosity runs can we expect?

* Sites are simulated as heterozygous (1) or homozygous (0).
* Mutational events occur at random.
* Double-stranded breaks (DSBs) also occur at random; direction and length of
  tract are draws as well.

TODO:

* 1000 and 10000bp chromosomes
* statistics at dynamic equilibrium
* statistics approaching equilibrium
* crossing over
* vary chromosome architecture
* Doxygen documentation on everything
