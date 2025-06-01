
from functools import partial
import os

from carabiner.mpl import add_legend, grid, colorblind_palette
import gradio as gr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import multivariate_hypergeom, nbinom

# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
    color=("lightgrey", "dimgrey") + colorblind_palette()[1:],
) 

SEED: int = 42
MAX_TIME: float = 5.
SOURCES_DIR: str = "sources"


def inject_markdown(filename):
    with open(os.path.join(SOURCES_DIR, filename), 'r') as f:
        md = f.read()
    return gr.Markdown(
        md,
        latex_delimiters=[
            {"left": "$$", "right": "$$", "display": True},
            {"left": "$", "right": "$", "display": False},
        ],
    )


def lotka_volterra(t, y, w, K):
    remaining_capacity = np.sum(y) / K
    dy = w * y.flatten() * (1. - remaining_capacity)
    return dy


def grow(t, n0, w, K=None):
    """Deterministic population size at time t.

    n0 : initial cells
    w  : growth rate
    t  : time (arbitrary units)
    K  : carrying capacity  (None → pure exponential)

    """
    if K is None:
        return n0 * np.exp(w[None] * t[:,None])
    # logistic with shared K
    else:
        ode_solution = solve_ivp(
            lotka_volterra,
            t_span=sorted(set([0, max(t)])),
            t_eval=sorted(t),
            y0=n0,
            vectorized=True,
            args=(w, K),
        )
        # print(ode_solution)
        return ode_solution.y


def plotter_t(x, growths, scatter=False, **kwargs):
    fig, axes = grid(aspect_ratio=1.5)
    plotter_f = partial(axes.scatter, s=5.) if scatter else axes.plot
    for i, y in enumerate(growths):
        plotter_f(
            x.flatten(), y.flatten(), 
            label=f"Mutant {i-1}" if i > 1 else "_none",
        )
    axes.set(
        xlabel="Time",
        yscale="log",
        **kwargs,
    )
    add_legend(axes)
    return fig


def plotter_ref(x, growths, scatter=False, fitlines=None, text=None, **kwargs):
    fig, axes = grid(aspect_ratio=1.5 if text is None else 1.7)
    plotter_f = partial(axes.scatter, s=5.) if scatter else axes.plot
    for i, y in enumerate(growths):
        plotter_f(
            x.flatten(), y.flatten(), 
            label=f"Mutant {i-1}" if i > 1 else "_none",
        )
    if fitlines is not None:
        fit_x, fit_y = fitlines
        for i, b in enumerate(fit_y.flatten()):
            y = np.exp(np.log(fit_x) @ b[None])
            print(fit_x.shape, y.shape, b)
            axes.plot(
                fit_x.flatten(), y.flatten(),
                label="_none",
            )
    if text is not None:
        axes.text(
            1.05, .1, 
            text, 
            fontsize=10,
            transform=axes.transAxes,
        )
    axes.set(
        xscale="log",
        yscale="log",
        **kwargs,
    )
    add_legend(axes)
    return fig


def calculate_growth_curves(inoculum, inoculum_var, carrying_capacity, fitness, n_timepoints=100):
    inoculum_var = inoculum + inoculum_var * np.square(inoculum)
    p = inoculum / inoculum_var
    n = (inoculum ** 2.) / (inoculum_var - inoculum)
    w = [1., 0.] + list(fitness)
    n0 = nbinom.rvs(n, p, size=len(w), random_state=SEED)
    t = np.linspace(0., MAX_TIME, num=int(n_timepoints))
    growths = grow(t, n0, w, inoculum * carrying_capacity)
    ref_expansion = growths[0] / n0[0]
    return t, ref_expansion, growths

    
def growth_plotter(inoculum, inoculum_var, carrying_capacity, *fitness):
    t, ref_expansion, growths = calculate_growth_curves(inoculum, inoculum_var, carrying_capacity, fitness, n_timepoints=100)
    return [
        plotter_t(
            t, 
            growths, 
            ylabel="Number of cells per strain",
        ),
        plotter_ref(
            ref_expansion, 
            growths, 
            xlabel="Fold-expansion of wild-type",
            ylabel="Number of cells per strain",
        ),
    ]


def reads_sampler(population, sample_frac, seq_depth, reps, variance):
    samples = []
    for i, timepoint_pop in enumerate(np.split(population.astype(int), population.shape[-1], axis=-1)):
        sample_size = np.floor(timepoint_pop.sum() * sample_frac).astype(int)
        samples.append(
            multivariate_hypergeom.rvs(
                m=timepoint_pop.flatten(), 
                n=sample_size, 
                size=reps, 
                random_state=SEED + i,
            ).T
        )
    samples = np.stack(samples, axis=-2)
    read_means = np.floor(seq_depth * samples.shape[0] * samples / samples.sum(axis=0, keepdims=True))
    variance = read_means + variance * np.square(read_means)
    p = read_means / variance
    n = (read_means ** 2.) / (variance - read_means)
    return np.stack([
        nbinom.rvs(n[...,i], p[...,i], random_state=SEED + i) 
        for i in range(reps)
    ], axis=-1)


def fitness_fitter(read_counts, ref_expansion):
    
    read_count_expansion = read_counts / np.mean(read_counts[:,:1], axis=-1, keepdims=True)
    read_count_expansion_ref = read_count_expansion[:1]
    log_read_count_correction = np.log(read_count_expansion) - np.log(read_count_expansion_ref)

    ref_expansion = np.tile(
        np.log(ref_expansion)[:,None], 
        (1, log_read_count_correction.shape[-1]),
    ).reshape((-1, 1))
    betas = []
    for i, log_strain_counts_corrected in enumerate(log_read_count_correction):
        ols_fit = np.linalg.lstsq(a=ref_expansion, b=log_strain_counts_corrected.flatten())
        betas.append(ols_fit[0])

    return log_read_count_correction, np.asarray(betas)


def fitness_fitter_spike(log_read_count_corrected):
    log_spike_count_corrected = log_read_count_corrected[1,...].flatten()[...,None]
    betas = []
    for i, log_strain_counts_corrected in enumerate(log_read_count_corrected):
        ols_fit = np.linalg.lstsq(
            a=log_spike_count_corrected, 
            b=log_strain_counts_corrected.flatten(),
        )
        betas.append(ols_fit[0])

    return log_spike_count_corrected, np.asarray(betas)


def reads_plotter(
    sample_frac, seq_reps, seq_depth, read_var,
    inoculum, inoculum_var, carrying_capacity, *fitness
):
    t, ref_expansion, growths = calculate_growth_curves(inoculum, inoculum_var, carrying_capacity, fitness, n_timepoints=10)
    read_counts = reads_sampler(growths, sample_frac, seq_depth, seq_reps, read_var)
    log_read_count_correction, betas = fitness_fitter(read_counts, ref_expansion)
    plot_text = "\n".join(
        f"Mutant {i-1}: $w_{i-1}/w_{'{wt}'}={1. + b:.2f}$"
        for i, b in enumerate(betas.flatten()) if i > 1
    )
    log_spike_count_corrected, spike_betas = fitness_fitter_spike(log_read_count_correction)
    plot_text_spike = "\n".join(
        f"Mutant {i-1}: $w_{i-1}/w_{'{wt}'}={1. - b:.2f}$"
        for i, b in enumerate(spike_betas.flatten()) if i > 1
    )
    read_count_correction = np.exp(log_read_count_correction)
    return growth_plotter(inoculum, inoculum_var, carrying_capacity, *fitness) + [
        plotter_t(
            np.tile(t[:,None], (1, seq_reps)), 
            read_counts, 
            scatter=True, 
            ylabel="Read counts per strain",
        ),
        plotter_ref(
            np.tile(ref_expansion[:,None], (1, seq_reps)), 
            read_count_correction, 
            scatter=True, 
            fitlines=(ref_expansion[:,None], betas),
            text=plot_text,
            xlabel="Fold-expansion of wild-type",
            ylabel="$\\frac{c_1(t)}{c_{wt}(t)} / \\frac{c_1(0)}{c_{wt}(0)}$",
        ),
        plotter_ref(
            read_count_correction[1:2,...], 
            read_count_correction, 
            scatter=True, 
            fitlines=(read_count_correction[1:2,...].flatten()[...,None], spike_betas),
            text=plot_text_spike,
            xlabel="$\\frac{c_{spike}(t)}{c_{wt}(t)} / \\frac{c_{spike}(0)}{c_{wt}(0)}$",
            ylabel="$\\frac{c_1(t)}{c_{wt}(t)} / \\frac{c_1(0)}{c_{wt}(0)}$",
        ),
    ]

with gr.Blocks() as demo:
    inject_markdown("header.md")
    # Growth curves
    inject_markdown("growth-curve-intro.md")
    mut_fitness_defaults = [.5, 2., .2]
    with gr.Row():
        relative_fitness = [
            gr.Slider(0., 3., step=.1, value=w, label=f"Relative fitness, mutant {i + 1}")
            for i, w in enumerate(mut_fitness_defaults)
        ]
    with gr.Row():
        n_mutants = len(mut_fitness_defaults) 
        inoculum = gr.Slider(
            10, 1_000_000, 
            step=10, 
            value=1000, 
            label="Average inoculum per strain",
        )
        inoculum_var = gr.Slider(
            .001, 1., 
            step=.001, 
            value=.001, 
            label="Inoculum variance between strains",
        )
        carrying_capacity = gr.Slider(
            len(mut_fitness_defaults) + 1, 10_000, 
            step=1, value=10, 
            label="Total carrying capacity (x inoculum)",
        )
    plot_growth = gr.Button("Plot growth curves")
    growth_curves_t = gr.Plot(label="Growth vs time", format="png")

    inject_markdown("growth-curve-t-independent.md")
    growth_curves_ref = gr.Plot(label="Growth vs WT expansion", format="png")
    growth_curves = [growth_curves_t, growth_curves_ref]

    # Read counts
    inject_markdown("read-counts-intro.md")
    with gr.Row():
        sample_frac = gr.Slider(
            .001, 1., step=.001, 
            value=.1, 
            label="Fraction of population per sample",
        )
        seq_reps = gr.Slider(
            1, 10, 
            step=1, 
            value=3, 
            label="Technical replicates",
        )
        seq_depth = gr.Slider(
            10, 10_000, 
            step=10, 
            value=10_000, 
            label="Average reads per strain per sample",
        )
        read_var = gr.Slider(
            .001, 1., 
            step=.001, 
            value=.001, 
            label="Sequencing variance",
        )
    plot_reads = gr.Button("Plot read counts")
    read_curves_t = gr.Plot(label="Read counts vs time", format="png")

    inject_markdown("read-counts-expansion.md")
    read_curves_ref = gr.Plot(label="Read count diff vs WT expansion", format="png")

    inject_markdown("read-counts-spike.md")
    read_curves_t2 = gr.Plot(label="Read count diff vs spike count diff", format="png")

    read_curves = [
        read_curves_t, 
        read_curves_ref,
        read_curves_t2,
    ]

    # Events
    plot_growth.click(
        fn=growth_plotter,
        inputs=[inoculum, inoculum_var, carrying_capacity, *relative_fitness],
        outputs=growth_curves,
    )
    plot_reads.click(
        fn=reads_plotter,
        inputs=[sample_frac, seq_reps, seq_depth, read_var] + [inoculum, inoculum_var, carrying_capacity, *relative_fitness],
        outputs=growth_curves + read_curves,
    )

demo.launch(share=True)
