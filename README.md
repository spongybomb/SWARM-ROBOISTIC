# SWARM-ROBOISTIC

Supporting code for the TAIG 2026 workshop paper (Paper 0012):

> Tianyi Zhang. *A Governance Reading of the Connectivity Floor: Bounding Human Oversight Latency in Decentralized Robotic Swarms.* TAIG Workshop 2026.

The paper argues that a floor on the algebraic connectivity λ₂ of a swarm's communication graph also bounds **oversight latency**: the number of communication rounds a human override needs to reach every agent. It uses the spectral diameter bound of Chung, Faber and Manteuffel (1994):

```
T_override = ecc(s) ≤ diam(G) ≤ floor( acosh(N−1) / acosh((λ_N + λ₂)/(λ_N − λ₂)) ) + 1
```

## Status

The paper is a position paper and reports no experimental results. This repository contains:

1. the analytic design curve used in the talk, and
2. a **preliminary static sanity check** on random graphs.

It does **not** contain the evaluation protocol proposed in §5 of the paper: moving robots, the barrier-function controller, flocking baselines, node loss and estimated λ₂. That work is ongoing.

## Files

| File | What it does |
|---|---|
| `plot_design_curve.py` | Analytic bound vs. connectivity floor λ_min for N = 200, λ_N ≤ 16 (at most 8 links per robot). Marks the floor needed for a 20-round requirement. |
| `sanity_check_200.py` | Static check on 20 random 200-robot graphs (see below). Prints a per-seed table and saves a figure. |
| `plot_sanity_for_slide.py` | Re-runs the sanity check and saves a 16:9 slide version of the figure. |
| `env_check.py` | Checks that the required packages import and that matplotlib can save a figure. |
| `results/` | Figures and console output from the runs described below. |

Code comments are in Chinese.

## Sanity check: method

For each of 20 random seeds:

1. Place N = 200 robots uniformly at random in a 20 m × 20 m area.
2. Connect two robots if they are within radio range R and each is among the other's 8 nearest in-range neighbours, so the maximum degree is ≤ 8 (unweighted graph).
3. Inject the order at the robot nearest the operator terminal at (0, 0).
4. Forward the order as a discrete flag, one hop per round, and count rounds until every robot has it.
5. Compare the measured rounds with the source eccentricity, the graph diameter and the spectral bound B(G) computed from the exact Laplacian eigenvalues.

The graph is frozen during propagation. Robots do not move, and no connectivity controller is applied.

## Sanity check: results

**R = 3.0 m** (default in the script). 16 of 20 graphs were connected. In all 16:

- measured rounds = eccentricity ≤ diameter ≤ bound held;
- measured rounds were 13–17, diameter 15–19 and the bound 36–57;
- λ₂ was 0.03–0.08, far below the floor of 0.5 used as the paper's example.

The bound holds on these graphs but is loose (roughly 2.5–3× the measured value). Of the 4 disconnected graphs, one let the order reach only 5 of 200 robots.

**R = 2.0 m** (the radio range in the paper's scenario). To reproduce, set `R = 2.0` in `sanity_check_200.py`. 0 of 20 graphs were connected, and on average the order reached about two thirds of the robots. Random placement alone does not guarantee that an override reaches everyone. That is the case for actively maintaining a connectivity floor.

Full per-seed output for R = 3.0 m: `results/sanity_check_output.txt`.

## How to run

Tested with Python 3.14 on Windows:

```bash
pip install -r requirements.txt
python env_check.py
python plot_design_curve.py
python sanity_check_200.py
python plot_sanity_for_slide.py
```

Figures are written to the current folder. They are deterministic for the fixed seeds 0–19.

## Note on authorship

The analysis scripts in this repository were written with AI assistance (Hermes Agent). The author ran them and reviewed the results. The paper itself and its argument are the author's own work.

## Contact

Tianyi Zhang, School of Data Science, Lingnan University. tianyizhang23@ln.hk

## Static experiments for the poster

`static_experiments.py` runs three frozen-topology experiments. Results go to `results/static_experiments.json`. `poster_figures.py` then draws them to `results/poster/`. The graph model is the same as in the sanity check: 200 robots in a 20 m × 20 m area, at most 8 mutual nearest-neighbour links, 20 seeds, and injection at the robot nearest (0, 0).

| Experiment | What varies | Result |
|---|---|---|
| E1 radio-range sweep | R = 1.5–5 m | 0/20 networks connected at 2 m (mean reach 67%) and 16/20 at 3 m. Above 3.5 m the 8-link cap saturates, so the graphs are identical. |
| Bound check | 45 distinct connected networks, R = 2.5–3.5 m | measured rounds = eccentricity ≤ diameter ≤ bound in 45/45 networks; the bound is 2.4–4.7× the measured value |
| E2 flag vs consensus | same networks and source; consensus with β = 1/9, threshold 0.95 | flag: median 16 rounds; consensus: median 4,468 rounds. Consensus rounds match ln(20)/(β·λ(L_g)) within 9%. The check λ(L_g) ≤ d_s/(N−1) held for all networks. |
| E3 node loss, 15% | none / random / regional (nearest to (15, 15) m), no re-linking | connected networks out of 16: 16 / 11 / 4 |

The consensus rounds are computed from the eigendecomposition of the grounded Laplacian. On two networks they were checked against direct iteration and matched exactly.

All of these are static snapshots: robots do not move and no connectivity controller is applied.

## Moving simulation: three disaster scenes, with and without a connectivity floor

`cbf_sim.py` is a reference implementation of the mechanism in §4 of the paper. `test_cbf_sim.py` checks three things:

- the analytic λ₂ gradient against finite differences;
- the per-robot QP against brute-force grid search;
- the spectral bound on 300 random weighted graphs.

`poster_sim_figures.py` draws the results.

**Setup.** 200 single-integrator robots start in a 6 m × 6 m entrance area and search for 60 s. The radio range is 2 m, and edge weights decay smoothly with distance. At 60 s the commander's evacuation order is injected at the robot nearest (0, 0). The order is then forwarded as a flag, one hop per control step, while the robots keep moving. There are three scenes, all geometric abstractions:

- a collapsed-building void, 20 × 20 m;
- a collapsed tunnel, 60 × 6.7 m;
- separated survivor pockets: three rooms in a 30 × 30 m area.

Each scene is run with no floor, λ_min = 0.05 and λ_min = 0.2. Each setting uses 20 seeds, and the start positions and targets are the same with and without the floor.

**Simplifications.**
- λ₂ and the Fiedler vector are computed exactly and centrally, so A3 is assumed.
- There is no 8-link cap. The bound is therefore evaluated with the weighted λ₂ and the λ_N of the support graph, which is valid because L_w ≼ L_support.
- The constraint is allocated as b_i ∝ ‖g_i‖².
- The floor is enforced with a 15% margin.
- A discrete-time safety step halves the step size whenever the next step would push λ₂ below the floor. This was needed because λ₂ became repeated: the minimum spectral gap was 5e-7, which is the regime that A2 warns about.
- There is no collision guarantee.

**Results.**

| Scene | Every robot got the order: no floor / 0.05 / 0.2 | Mean share of robots reached, no floor | Area searched: no floor / 0.05 / 0.2 |
|---|---|---|---|
| Building void | 3/20 / 20/20 / 20/20 | 90% | 96% / 56% / 38% |
| Collapsed tunnel | 2/20 / 20/20 / 20/20 | 55% | 61% / 18% / 18% |
| Survivor pockets | 0/20 / 20/20 / 20/20 | 34% | 69% / 57% / 47% |

- In the 120 runs with a floor, λ₂ never fell below the floor.
- The measured override rounds never exceeded the spectral bound at injection time.
- In 308 agent-steps (out of 7.2 million) the local QP was infeasible. These steps were logged and handled on a best-effort basis.

Run with `python cbf_sim.py 20 12` (20 seeds, 12 worker processes; about 6 minutes on 24 cores). Each run is appended to `results/cbf_sim_runs.jsonl`, so an interrupted job can be resumed.
