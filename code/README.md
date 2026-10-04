Ancillary code for "Cognitive Accessibility and the Fermi Problem" (Paper D of
Rungs, Voids and Effective Closure; Programme 1 of Effective Horizons and
Cognitive Spaces), working paper v1, October 2026.  Python 3, numpy, scipy,
matplotlib.  All scripts run from this directory.

marked_process.py  the marked birth-and-contact process (Section 3); run(fam, ...)
                   returns the graph statistics on V(t) with the active/passive
                   access gates (active=, kappa0=) of Section 4.
starttimes.py      cosmic start-time model and the Fermi-function figures (Section 2).
experiments.py     kernel robustness and the first factorial/phase grids.
sobol.py           Saltelli/Jansen Sobol indices (pooled, superseded) and the
                   local elasticities used by elast06.py.
sobol3.py          within-family Sobol indices with bootstrap intervals, grouped
                   effects, R_stoch and R_mean (Figure 3): python sobol3.py <family>
sobol3_fig.py      draws Figure 3 from the sobol3_*.json outputs.
elast06.py         local log-elasticities on V(t) (Figure 4).
phase4.py          persistence-bandwidth phase diagrams (Figures 5-6).
nscale.py, nscale2.py  finite-size sweep of the persistence crossover (Figure 7).
channels.py        active/passive channel table (Section 4).
filters.py         nested accessibility regions for one observer (Figure 8):
                   python filters.py hp  (high-persistence setting of the figure)

Results
-------
results/ holds the JSON outputs behind the figures and tables of Paper D, so
that the figures can be redrawn without rerunning the simulations.  To run a
script that reads an earlier output (phase4.py reads phase3.json;
sobol3_fig.py reads sobol3_*.json), copy the file from results/ into this
directory first.  Scripts write figures one directory up (../*_fig.pdf).

Requirements: Python 3.10+, numpy, scipy, matplotlib.  Runtimes range from
seconds (channels.py) to about five minutes per family (sobol3.py).
