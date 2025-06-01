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
