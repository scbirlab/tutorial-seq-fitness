But we don't actually know the true fold-expansion of the reference strain, since 
it's not directly observed. 

But we _do_ know some things about it. For example, when the pool is far from carrying 
capacity, the fold-expansion will be close to exponential:

$$
\log \frac{n_{wt}(t)}{n_{wt}(0)} = w_{wt} t
$$

And at the carrying capacity, the fold-expansion stops changing with time, arrested
at the carrying capacity minus the other cells in the pool:

$$
\log \frac{n_{wt}(t)}{n_{wt}(0)} = \log \frac{K - \Sigma_j n_j}{n_{wt}(0)}
$$

So at very early timepoints, long before carrying capacity is reached,

$$
\log \left( \frac{c_i(t)}{c_{wt}(t)}\frac{c_{wt}(0)}{c_i(0)} \right) = \left(\frac{w_i}{w_{wt}} - 1 \right) w_{wt} t
$$

or equivalently,

$$
\log \frac{c_i(t)}{c_{wt}(t)} = \log \frac{c_i(0)}{c_{wt}(0)} + (w_i - w_{wt}) t
$$

So at early timepoints, the log-ratio of a strain's counts to reference counts increases linearly 
over time with the fitness difference between the strain and the reference. The fitness difference 
is useful, but the ratio would be better. (Alternatively, we could use a known value of the reference 
fitness measured separately).

And after carrying capacity is reached,

$$
\log \frac{c_i(t)}{c_{wt}(t)} = \log \frac{c_i(0)}{c_{wt}(0)} + \left(\frac{w_i}{w_{wt}} - 1 \right) \log \frac{K - \Sigma_j n_j}{n_{wt}(0)}
$$

So in this regime, the log-ratio of a strain's counts to reference counts is fixed. 
Subtracting the log-ratio of a strain's counts to reference counts in the input leaves:

$$
\log \frac{c_i(t)}{c_{wt}(t)} - \log \frac{c_i(0)}{c_{wt}(0)} = \left(\frac{w_i}{w_{wt}} - 1 \right) \log \frac{K - \Sigma_j n_j}{n_{wt}(0)}
$$

But we still can't get the relative finess without dividing by the constant

$$
\log \frac{K - \Sigma_j n_j}{n_{wt}(0)}
$$

which we don't know directly. However, the final trick is to use a non-growing control, so that 

$$\frac{w_i}{w_{wt}} = 0$$

That means that for this control, 

$$
\log \frac{c_i(t)}{c_{wt}(t)} - \log \frac{c_i(0)}{c_{wt}(0)} = - \log \frac{K - \Sigma_j n_j}{n_{wt}(0)}
$$

We can then get the relative fitness ratio for the other strains directly.

