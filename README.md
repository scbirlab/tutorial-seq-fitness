---
title: Tutorial - Fitness estimation from pooled growth and NGS
emoji: 🧮
colorFrom: green
colorTo: blue
sdk: gradio
sdk_version: 5.22.0
app_file: app.py
pinned: false
license: mit
short_description: How do strains grow when competing with each other? And how can we infer their fitness from next-generation sequencing data?
tags:
    - biology
    - sequencing
---

# Tutorial – Fitness estimation from pooled growth

[![Open in Spaces](https://huggingface.co/datasets/huggingface/badges/resolve/main/open-in-hf-spaces-md-dark.svg)](https://huggingface.co/spaces/scbirlab/tutorial-seq-fitness)

This is the repository for the interactive tutorial, acccessed [here](https://huggingface.co/spaces/scbirlab/tutorial-seq-fitness). Non-interactive text from the tutorial is below.

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

### Using spike-in counts

But we don't actually know the true fold-expansion of the reference strain, since 
it's not directly observed. However, a non-growing fitness-zero control can help,
such as a heat-killed strain or a spike-in plasmid.

We start with the equation before,

$$
\log \left( \frac{c_i(t)}{c_{wt}(t)}\frac{c_{wt}(0)}{c_i(0)} \right) = \left(\frac{w_i}{w_{wt}} - 1 \right) \log \frac{n_{wt}(t)}{n_{wt}(0)}
$$

But for the fitness-zero control, $w_{spike} = 0$, so:

$$
\log \left( \frac{c_{spike}(t)}{c_{wt}(t)}\frac{c_{wt}(0)}{c_{spike}(0)} \right) = -\log \frac{n_{wt}(t)}{n_{wt}(0)}
$$

This means that, although we don't know how the reference strain grows directly, its
growth is given from the ratio of the spike counts to the reference counts, normalized
to the same ratio at time 0.

This leaves us with the overall equation:

$$
\log \left( \frac{c_i(t)}{c_{wt}(t)}\frac{c_{wt}(0)}{c_i(0)} \right) = \left(1 - \frac{w_i}{w_{wt}} \right) \log \left( \frac{c_{spike}(t)}{c_{wt}(t)}\frac{c_{wt}(0)}{c_{spike}(0)} \right)
$$

If we plot the left hand side against the right, we should get a straight line for 
each strain with intercept zero and gradient $1 - \frac{w_i}{w_{wt}}$.
