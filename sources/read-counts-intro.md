## Read counts from next-generation sequencing

But we don't actually measure the number of cells directly. Instead, we're measuring the
number of reads (or UMIs) which represent a random sampling of the population followed by
molecular biology handling and uneven sequencing per lane which decouples the relative
abundances for each timepoint. 

Below, you can simulate read counts for technical replicates of the growth curves above.
The simulation:

1. Randomly samples a defined fraction of the cell population (without replacement, i.e.
the [Hypergeometric distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution)).
Smaller samples from smaller populations are noisier.
2. Calculates the resulting proportional representation of every strain in every sample.
3. Multiplies that proportion by read depth.
4. Randomly samples sequencing read counts resulting from variations in library construction
and other stochasticity, according to the 
[Negative Binomial distribution](https://en.wikipedia.org/wiki/Negative_binomial_distribution), 
an established noise model for sequencing counts.
