# Effective Horizons and Cognitive Spaces

*What a mind can reach, in what time, and who can tell.*

**Programme 1: Rungs, Voids and Effective Closure**, four working papers and a synopsis, by Nial Goodman (version 1, October 2026).

The programme models the reach of a cognitive agent as a Turing ideal. The problems an agent can approach but not settle are the limits it can effectively generate from its own approximations, and iterating this *effective closure* produces a ladder of rungs (exactly the arithmetical hierarchy) with voids between them. The ladder exists only because convergence is uncertified: require a modulus and it collapses to one rung. The four papers develop the framework, its geometry, its dynamics in time and strategy, and an application to the Fermi problem.

## Papers

Start with the synopsis, then Paper B. A, C and D can be read in any order after B; A is the most technical and can be skipped on a first reading.

| | Paper | Pages |
|---|---|---|
| 0 | [**Synopsis**: where the idea came from, what each paper does, how they fit, how to read them](papers/0_Synopsis.pdf) | 3 |
| B | [**Rungs, Voids and Effective Closure**: cognitive spaces, uncertified convergence and promotion (the core framework)](papers/B_Rungs_Voids_Effective_Closure.pdf) | 22 |
| A | [**A Jump Metric on the Arithmetical Degrees**: rays, a shadow tree, dead ends, and the question of tree-likeness](papers/A_Jump_Metric_Arithmetical_Degrees.pdf) | 17 |
| C | [**Populations and Strategy in Cognitive Spaces**: lag, horizons, finite computation and the Gambler](papers/C_Populations_and_Strategy.pdf) | 12 |
| D | [**Cognitive Accessibility and the Fermi Problem**: a rate model, a Fermi function and a marked birth-and-contact process](papers/D_Cognitive_Accessibility_Fermi.pdf) | 16 |

arXiv identifiers will be added here when available.

## Selected results

- **Rung and Stall Theorems (B).** Iterated effective closure from the computable sets gives exactly the finite arithmetical levels. Their union is a fixed point, and true arithmetic lies outside it.
- **Indistinguishability (B).** No observer, of any computational power, can verify in the limit that a process has left the closure of a countable cognitive space. Next-rung content is verifiable from two rungs up.
- **Subsumption vs superiority (B).** An agent can come to believe correctly, in the limit, that a peer lies within its space. No agent can ever come to believe correctly that a peer exceeds it.
- **Hidden information is unbounded (A).** There are points whose first jump conceals arbitrarily much while their horofunction level is zero.
- **The Gambler (C).** An agent playing from a source it does not hold is unpredictable even when its cognition is fully subsumed. A random source can be audited in hindsight three rungs up but never predicted.
- **Five filters (D).** Between a physical population and an observed one lie existence, persistence, causal reach, channel access and cognitive access. For a marked Poisson population these give P(silence) = exp(−Λ_acc), so silence constrains an accessible intensity, not abundance.

## Code

[`code/`](code/) contains the simulation code for Paper D, with the JSON outputs behind its figures in [`code/results/`](code/results/). See [`code/README.md`](code/README.md) for which script produces which figure.

## Method

The ideas and questions are the author's. The formalisation, proofs, simulations and drafting were developed with Claude, an AI system by Anthropic, under the author's direction, and the drafts were read adversarially by a separate AI reviewer. Three errors caught in that process (a false lemma and two modelling artefacts) are kept in the papers as findings; the synopsis describes them. The proofs have not been independently refereed. Corrections, counterexamples and prior-art references are very welcome: please open an issue.

## Citation

> N. Goodman, *Rungs, Voids and Effective Closure* (Papers 0, A–D), working papers, version 1, October 2026. https://github.com/n-goodman/effective-horizons

## Licence

Papers: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).

Contact: nialgoodman@gmail.com
