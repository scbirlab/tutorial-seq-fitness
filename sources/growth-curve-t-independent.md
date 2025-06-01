### Removing time dependence

It can be difficult to get absolute fitness out of these curves, because 
when the pool approaches the carrying capacity, all the strains growth rates 
mutually affect each other. 

However, if we're only interested in the relative fitness of multiplexed strains relative
to a reference (e.g. wild-type) strain, then we can make this simplification:

$$
\frac{dn_{i}}{dt} / \frac{dn_{wt}}{dt} = \frac{dn_{i}}{dn_{wt}} = \frac{w_{i} n_i}{w_{wt} n_{wt}}
$$

The interdependency term cancels out, and time is removed, with the reference strain's
growth acting as the clock. Unlike the time-dependent Lotka-Volterra equations, this
has a closed-form integral:

$$
\log n_i(t) = \frac{w_i}{w_{wt}} \log \frac{n_{wt}(t)}{n_{wt}(0)} + \log{n_i(0)}
$$

So now the log of the number of cells of a mutant at any moment ($n_1(t)$) is dependent only
on its inoculum ($n_1(0)$), how much the reference strain has grown (i.e. fold-expansion, 
$\frac{n_{wt}(t)}{n_{wt}(0)}$), and the ratio of fitness between the mutant and the 
reference ($\frac{n_{wt}(t)}{n_{wt}(0)}$).

Plotting the number of cells of a mutant against the reference fold-expansion gives 
a straight line (on log scales), with intercept being the mutant inoculum and gradient
being its relative fitness with respect to the reference strain's fitness.

