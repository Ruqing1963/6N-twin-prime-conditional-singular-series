# 6N Twin-Prime Conditional Singular Series (Part III)

A first-order **conditional singular series** for the twin-prime gap
distribution on the 6N ± 1 skeleton, resolving the open problem left by Part II,
together with its **second-order residual**.

**Background.** Part I (conditional density, factor ∏(q−1)/(q−3)) and Part II
(the ω-dependent shift of the gap distribution and the congruence-lockdown
mechanism) showed that the twin-gap distribution depends on ω₍>3₎(N), the number
of distinct prime factors > 3 of the centre N. Part II left the closed-form
conditional gap singular series as an open problem.

**This part.** We give its explicit first-order form. Conditioning the
Hardy–Littlewood correlation factor on the centre's factor set via the
congruence-lockdown rule (q | N forbids the gap residues d ≡ ±6⁻¹ mod q, i.e.
q | (6d ± 1)):

```
    S_{2,ω}(d)  ≈  C2(d) · L_ω(d)
    C2(d)   = ∏_{q>3} (1 − ν_q(d)/q)/(1 − 2/q)²        (unconditional H-L factor)
    L_ω(d)  = Pr[ q ∤ N for every q | (6d±1) | ω₍>3₎(N)=ω ]   (lockdown-allowance)
```

The admissibility of a gap is governed by the factorisation of 6d ± 1:
6·35∓1 = 209 = 11×19 (so 210 is locked by the common small primes 11, 19),
whereas 6·7∓1 = 41, 43 are both prime (so 42 is locked only by the rare primes
41, 43).

**Result (tested on the 23,988,173 twin centres of S₁₀).** The first-order
series reproduces the ω-trends of the gaps 30, 42, 60 (residuals ≤ 0.13) and the
**direction** of the collapse of 210. A systematic residual remains at high ω on
210: observed falls to 0.41 where the first-order series predicts 0.62
(residual ≈ −0.2 at ω = 5,6, stable across S₉ and S₁₀). We attribute it to the
**second-order correlation** between the factorisations of the two centres N and
N + d, which the first-order (per-prime independent) series omits, and pose the
two-centre conditional singular series as the next open problem.

> **Scope.** Experimental / computational number theory. The first-order series
> reduces to the unconditional Hardy–Littlewood series in the mean. No claim is
> made about the infinitude of twin primes or about any prime k-tuple conjecture.
> This advances the Part II open problem from fully open to first-order-solved
> with a quantified second-order remainder.

Part I: Zenodo doi:10.5281/zenodo.20470367 ·
Part II: Zenodo doi:10.5281/zenodo.20477664

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   ├── css_plot_S9.csv      6dN, d, omega, obs_r, pred_r, resid  (S9)
│   └── css_plot_S10.csv     same, S10 (the paper's headline numbers)
├── code/
│   ├── cond_singular_series.py     derive + evaluate S_{2,ω}(d); emits css_plot_S{K}.csv
│   └── make_cond_singular_fig.py   builds the 4-panel obs-vs-prediction figure
├── figures/                fig_paper3_cond_singular.{pdf,png}
└── paper/                  Chen_6N_Paper3.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Derive + evaluate the first-order conditional singular series.
#    Output: console table (obs/pred/resid) + css_plot_S{K}.csv
MAXK=9  python code/cond_singular_series.py    # S9  (~1-2 min)
MAXK=10 python code/cond_singular_series.py    # S10 (~15 min)

# 2. Build the showdown figure (reads ../data/css_plot_S10.csv)
cd code && python make_cond_singular_fig.py
```

### Definitions / conventions (same as Part II)

- Twin centre: N with 6N−1 and 6N+1 both prime. Gap ΔN = N_{i+1} − N_i between
  consecutive twin centres (centre-step units); physical distance is 6ΔN.
- The gap is attributed to ω₍>3₎ of the left centre.
- Only strata with ≥ 10⁴ gaps are reported (ω = 1…6 in S₁₀); ω ≥ 7 is
  noise-dominated and excluded.
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality, self-verified against sympy in Part I/II; the S₁₀ twin count
  23,988,173 matches Part I to the integer.

## License

MIT — see `LICENSE`.
