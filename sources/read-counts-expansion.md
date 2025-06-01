### Accounting for sequencing subsampling per sample

Each sequencing sample $s$ could be over- or under-sampling the population relative to the first
timepoint by some factor $\phi_s$.

$$\log \frac{c_i(t)}{c_i(0)} = \log \phi_s\frac{n_i(t)}{n_i(0)} = \log \phi_s + \frac{w_i}{w_{wt}} \log \frac{n_{wt}(t)}{n_{wt}(0)}$$

Variables:
- $c_i(t)$: Read (or UMI) count of strain $i$ at time $t$
- $\phi_s$: The ratio of sampling depth at time $t$ to that at time $0$ for sample $s$

The factor $\phi_s$ is the ratio of _the ratio of read counts between samples_ 
and _the ratio of cell counts between samples_ for any strain (assuming each strain
is sampled without bias):

$$\log \phi_s = \log \frac{c_i(t)}{c_i(0)} - \log \frac{n_i(0)}{n_i(0)}$$

We can get rid of the nuisance parameter $\phi_s$ (which is difficult to measure becuase
we don't know the true number of cells for each strain and sample) using the following trick.

We have the equation for read counts for mutant $i$ (same as above):

$$
\log \frac{c_i(t)}{c_i(0)} = \log \phi_s + \frac{w_i}{w_{wt}} \log \frac{n_{wt}(t)}{n_{wt}(0)}
$$

And for the reference strain (relative fitness is 1):

$$
\log \frac{c_{wt}(t)}{c_{wt}(0)} = \log \phi_s + \log \frac{n_{wt}(t)}{n_{wt}(0)}
$$

We can make $\phi_s$ disappear by taking the difference:

$$
\log \frac{c_i(t)}{c_i(0)} - \log \frac{c_{wt}(t)}{c_{wt}(0)} = \frac{w_i}{w_{wt}} \log \frac{n_{wt}(t)}{n_{wt}(0)} - \log \frac{n_{wt}(t)}{n_{wt}(0)}
$$

This is equivalent to:

$$
\log \left( \frac{c_i(t)}{c_{wt}(t)}\frac{c_{wt}(0)}{c_i(0)} \right) = \left(\frac{w_i}{w_{wt}} - 1 \right) \log \frac{n_{wt}(t)}{n_{wt}(0)}
$$

So the ratio of _the count ratio of a strain to the reference strain at time t_ to 
_the count ratio of a strain to the reference strain at time 0_ is
dependent only on the relative fitness and the true fold-expansion of the reference strain.

Plotting the ratio of _the count ratio of a strain to the reference strain at time t_ to 
_the count ratio of a strain to the reference strain at time 0_ 
should give a straight line (on a log-log) plot, with intercept 0 and gradient equal to the relative fitness minus 1.
