## Multiplex growth curves

How do strains grow when competing with each other?

That's given by the [Lotka–Volterra competition model](https://en.wikipedia.org/wiki/Competitive_Lotka%E2%80%93Volterra_equations). 
For two strains:

$$
\frac{dn_{wt}}{dt} = w_{wt} n_{wt} \left( 1 - \frac{n_{wt} + n_{1}}{K} \right)
$$

$$
\frac{dn_1}{dt} = w_{1} n_1 \left( 1 - \frac{n_{wt} + n_{1}}{K} \right)
$$


- $n_i(t)$: abundance of species (or strain) $i$ at time $t$.
- $w_i$: intrinsic (exponential) growth rate of species $i$.
- $K$: carrying capacity.

We can generalize to many strains. For each one:

$$
\frac{dn_i}{dt} = w_{i} n_i \left( 1 - \frac{\Sigma_j n_j}{K} \right)
$$

It's not possible to algebraically integrate these equations, since they are
circularly dependent on each other. But we can numerically integrate, to simulate 
multiplexed growth curves.
