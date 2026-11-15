# Reproduction report — Beamforming with Joint Phase and Time Array

Paper: J. Mo, A. AlAmmouri, S. Dong, Y. Nam, W.-S. Choi, G. Xu, J. Zhang,
arXiv:2502.00139v1 (31 Jan. 2025).

Environment: MATLAB R2026b Update 1, Windows 11, 16 threads, 32 GB, no toolbox.
Gap numbers **G1–G8** are local to this repository.

---

## 1 What is covered

| Job | Script | Figures | Output |
|---|---|---|---|
| B-1 | `run_b1.m` | Fig. 7 from the unconstrained optimum of eq. (4) | `out/b1_fig7_{2,4,8,16}ue.png`, `out/b1_summary.txt` |
| B-1b | `run_b1b.m` | Fig. 7 for the three candidate design rules | console only |
| B-3 | `run_b3.m` | Fig. 8, delay law + symmetry + the eq. (5) read-out | `out/b3_fig8_delay.png`, `out/b3_fig8.csv`, `out/b3_summary.txt` |
| B-2/B-5 | `run_b1c.m` | Fig. 7 with the Fig. 8 delays, phase from the closed form | `out/b1c_fig7_eq5{q,}.png`, `out/b1c_summary.txt` |
| B-5 | `run_b1d.m` | TTD step-size sweep | `out/b1d_step_sweep.png`, `out/b1d_summary.txt` |
| B-5 | `run_b1e.m` | straight ramp versus sawtooth | `out/b1e_ramp_vs_saw.png` |
| B-4 | `run_b2.m` | Fig. 9 | `out/b2_fig9_throughput.png`, `out/b2_throughput.csv`, `out/b2_summary.txt` |

Figs. 4–6 of the paper are prototype measurements and are out of scope. Table II is a corollary
of Fig. 8 and is not plotted.

---

## 2 Model

**Array and weights.** 16 vertical × 16 horizontal, one TTD per column, so the hardware has 16
pairs $(\varphi_m,\tau_m)$, $m=1\ldots16$. The analog weight of subcarrier $k$ is

$$[\mathbf p_k]_m=\frac{1}{\sqrt{N_vN_h}}\,e^{j(\varphi_m+2\pi f_k\tau_m)}\quad(m\ \text{repeated for the }N_v\ \text{rows}),$$

and the gain toward angle θ is $|p_k^\dagger a(\theta,f_k)|^2$ with
$a_{m,n}=e^{-j2\pi f_k(n d\cos\theta)/c}$ per column and row.

**Eq. (4) decouples completely.** With $b_k$ the ideal steering vector of the user owning
subcarrier $k$,

$$\sum_k\|\mathbf p_k-b_k\|^2=\sum_m\underbrace{\sum_k|p_{k,m}-b_{k,m}|^2}_{\text{depends only on }\varphi_m,\tau_m},$$

so for a given $\tau_m$ the phase has a closed form,
$\varphi_m=\angle\sum_k\overline{b_{k,m}}e^{j2\pi f_k\tau_m}$, and only a one-dimensional search
over $\tau_m$ is left. The iterative algorithm of the paper's ref. [6] is therefore not needed —
but see gap **G7** for why this closed form is *not* what produces Fig. 7.

**Two delay rules are used in this repository.** (i) The unconstrained optimum of eq. (4),
obtained with one zero-padded FFT over the whole delay period $1/\Delta f=8.33\ \mu$s
(resolution 0.02 ns). (ii) The paper's own law from Fig. 8, which is what Fig. 7 must use.

---

## 3 Fig. 8: the delay profile

The figure shows four sawtooths, slope 2.5 ns per column, wrapping every $N_{UE}$ columns
(8/4/2/1 teeth for `N_UE` = 2/4/8/16), peak delay 2.5 / 7.5 / 17.5 / 35 ns.

Eq. (5) gives the slope $\sin(\Delta\theta/2)/W$: 1.25 ns for 2 and 4 users ($\pm30^\circ$,
$\Delta\theta=60^\circ$) and 2.05 ns for 8 and 16 users ($\pm55^\circ$). Neither is a multiple of
the 2.5 ns TTD step of Table I, so the hardware can only round **up**:

$$\tau_m=\underbrace{\left\lceil\frac{\sin(\Delta\theta/2)}{W}\right\rceil_{2.5\,\rm ns}}_{\text{slope of eq. (5), quantized}}\cdot\;\mathrm{mod}(m-1,\;N_{UE}),
\qquad \tau_{\max}=(N_{UE}-1)\cdot2.5\ \mathrm{ns}$$

| $N_{UE}$ | slope of eq. (5) | quantized | $\tau_{\max}$ with wrap | paper | difference |
|---|---|---|---|---|---|
| 2 | 1.25 | 2.50 | 2.50 | 2.50 | 0 |
| 4 | 1.25 | 2.50 | 7.50 | 7.50 | 0 |
| 8 | 2.05 | 2.50 | 17.50 | 17.50 | 0 |
| 16 | 2.05 | 2.50 | 37.50 | 35.00 | 2.50 (one TTD step) |

Two independent degrees of freedom must be hit at the same time — the wrap period $=N_{UE}$ and
the peak delay — and this law hits both. When $N_{UE}=M=16$ the modulo is the identity and the
sawtooth **is** eq. (5)'s straight ramp, so the paper's "with 16 users the type-1 beam
approaches the type-2 rainbow beam" is a corollary, and the two rules give bit-identical delay
vectors there (`run_b1e.m`).

`run_b3.m` additionally evaluates eq. (5) over the **whole** aperture (no wrap, no
quantization): 18.75 / 30.72 ns. That number is what an earlier version of this report used to
declare "eq. (5) cannot produce Fig. 8"; the claim is **withdrawn** — it compared one degree of
freedom out of two. The lesson is recorded because it is a general one: $\tau_{\max}$ alone does
not fix the convention, the wrap period must be matched at the same time.

---

## 4 Why the wrap is not cosmetic

Same slope, only the wrap removed:

| Worst in-band gain (dB) | 2 UEs | 4 UEs | 8 UEs | 16 UEs |
|---|---|---|---|---|
| straight ramp, eq. (5) as printed | −140.58 | −17.91 | −3.05 | −0.14 |
| sawtooth, as in Fig. 8 | **25.04** | **20.16** | **2.18** | −0.14 |

The wrap cuts the aperture into $16/N_{UE}$ teeth of 2.5 ns steps, and the jump between teeth is
where the multiple beams come from; a single ramp over the whole aperture forms one beam only.
Fig. 8's shape is a necessary condition for the multi-beam, not a plotting convention.

---

## 5 Fig. 7: beam patterns and the 3 dB claim

With the delays fixed by Fig. 8, the phase is determined in closed form, so there is no free
parameter left to fit:

$$\varphi_m=\angle\sum_k\overline{b_{k,m}}e^{j2\pi f_k\tau_m}$$

**(a) The 3 dB claim: two users hit, four and above do not.** Worst in-band gain of each user in
its own sub-band at its own angle (`out/b1c_summary.txt`):

| Rule | 2 UEs | 4 UEs | 8 UEs | 16 UEs |
|---|---|---|---|---|
| Fig. 8 delays (2.5 ns grid) | −140.58 | 0.65 | −24.03 | −12.34 |
| exact slope of eq. (5), continuous delay | **25.04** | 20.16 | 2.18 | −0.14 |
| per-subband ideal steering (upper bound, unrealisable) | 28.00 | 28.00 | 28.00 | 28.00 |

- Two users: 25.04 dB with 2.02 dB ripple — "at most 3 dB below 28 dB", exactly. This is the
  only scenario where B2/B5 is hit positively, and it is the scenario Fig. 2(a) uses to explain
  the type-1 codebook.
- Quantize the same slope to 2.5 ns and the same configuration collapses to −140.58 dB: the
  sub-band edge lands in a null. **Fig. 7(a)'s pattern can only come from a delay resolution
  finer than Table I allows** (gap **G5**).
- For 4/8/16 users, no rule reaches 25 dB; the best is 20.16 dB (4 users) and 8/16 users fall to
  2.18 / −0.14 dB. Even the ridge height is only 15–18 dB, so it is the multi-beam itself that
  is capped, not an unlucky user.
- Only the third row reaches 28 dB, and it needs one steering vector per sub-band; the hardware
  gives each element one frequency-flat phase shifter and one delay shared over the whole band,
  so $\mathbf p_k$ cannot change with $k$. Under a single $(\varphi,\tau)$ set the claim is
  unreachable (gap **G6**) — and that unreachable assumption is exactly what the system-level
  model of §6 adopts, as the paper does.

**(b) The pattern must be even, so Fig. 7 cannot be compared point by point.** The caption says
"JPTA **type-2** beam pattern". Reading the figure at high resolution, the 2-user panel has
about 25 dB at $-30^\circ$ and about 13 dB at $+30^\circ$ for $f-f_c=-150$ MHz — a 12 dB mirror
difference at one frequency. But for arbitrary weights

$$|AF(\theta,f)|=\Big|\sum_m p_m(f)\,e^{-j2\pi fdm\cos\theta/c}\Big|=|AF(-\theta,f)|,$$

verified numerically: $\max|G(\theta)-G(-\theta)|=0$ over 481 angles × 167 frequencies, and 0 for
random weights. Fig. 7 therefore contains something beyond eqs. (1)–(3) — a signed angle
convention, or weights applied per sub-band — so only the quantitative property (a) is compared
(gap **G4**). This reproduction's figures are symmetric in ±θ by construction.

**(c) The TTD step is not a smooth knob.** Sweeping the step as a continuous quantity
(`run_b1d.m`), the 2-user worst in-band gain oscillates between −10 and 27 dB on a 0.02–3 ns
grid: a small delay change rotates each element's phase across the band at a different rate, and
the beam alternates between pointing and falling into a null.

- the 1.25 ns of eq. (5) sits in a good window (25.04 dB);
- 2.5 ns — the only grid Table I allows — sits in a bad one (−140.58 dB);
- none of the four slopes (1.25, 1.25, 2.05, 2.05 ns) is a multiple of 2.5 ns, so no
  configuration can be quantized into a good window (4/8/16 users: 0.65 / −24.03 / −12.34 dB).

Verdict: **the delays of Fig. 8 and the patterns of Fig. 7(a) are mutually exclusive**, and the
2-user case already conflicts. This is not "a little less accurate": the type-1 codebook is
intrinsically sensitive to the delay value, which is also why the paper solves eq. (4)
iteratively instead of in closed form.

**(d) What the unconstrained optimum of eq. (4) actually gives.** Solving eq. (4) exactly
(`run_b1.m`) is not what produces Fig. 7:

| N_UE | $\tau_{\max}$ from eq. (4) | paper (Fig. 8) | worst in-band gain per user |
|---|---|---|---|
| 2 | 0.00 ns (all delays zero) | 2.50 | 4.48 / 1.17 dB |
| 4 | 2.50 | 7.50 | 14.17 / 24.88 / 24.04 / 17.01 |
| 8 | 7.50 | 17.50 | 11.0 / −1.45 / … / 11.16 |
| 16 | 17.50 | 35.00 | 9.79 / −2.21 / … / 7.05 |

For two users the "global optimum" is a plain phased array: all delays zero, one beam, worst
in-band 1.17 dB. The objective of eq. (4), minimising $\sum_k\|\mathbf p_k-b_k\|^2$, is happy to
give up the second beam because the two sub-bands are far apart and the squared error of the
dominant user dominates. This is the practical reason why the paper's Fig. 7 is evaluated **at
Fig. 8's delays**: the two figures describe one configuration, as the text says
("the corresponding configuration of the delay elements").

`run_b1b.m` compares the three design rules directly, and they separate exactly as predicted: at
8 users the "ls" rule returns $\tau_{\max}=8333.32$ ns — the period edge of G7 — and puts the
peaks at −47/−48/−47° instead of the requested −55/−39/−24°, whereas the sawtooth rule returns
17.50 ns, i.e. Fig. 8 exactly, and the rainbow rule 30.00 ns. Only the sawtooth reproduces the
configuration of Fig. 8.

---

## 6 Fig. 9: system level

**Link abstraction (Sec. IV-A).** Per-subcarrier SNR → EESM with η per MCS (Table 3 of the 2022
uplink paper) → BLER → throughput, with the MCS and the RB count chosen jointly:

$$\max_{m,\text{MCS}}\;m\cdot12\cdot\Delta f\cdot\mathrm{SE}\cdot(1-\mathrm{BLER})
\quad\text{s.t.}\quad \mathrm{BLER}\le10\%,\ m\ge4,$$

and the UE's **total** 23 dBm spread over the $m$ RBs it is granted, i.e.
$\mathrm{PSD}=P_{tx}-10\log_{10}(12m)$. That is the mechanism behind "fewer RBs can be faster",
and it is the strict form of the paper's $N_1$ slots × $N_2/k$ RBs versus $N_1/k$ slots × $N_2$
RBs.

**Two baselines.** PAA: all RBs, one slot out of $N_{UE}$ (`duty = 1/N_UE`), best of 16 beams.
JPTA: $1/N_{UE}$ of the RBs, every slot, in-band gain 28 dB dropping at most 3 dB at the
sub-band edge. Both therefore have the same average bandwidth $W/N_{UE}$; the difference is the
$N_{UE}$-fold energy concentration.

| Claim | Paper | Reproduced | Verdict |
|---|---|---|---|
| near distance: JPTA = PAA | equal | 0.98 / 1.00 / 0.98 / 1.00× | ✓ |
| 1500 m, 16 UEs | 830% (9.3×) | **8.05×** | ✓ 1.2× |
| 1500 m, gain vs users | 2/4/8/16× | 1.1 / 2.1 / 4.1 / 8.0× | trend ✓ |
| 1 Mbps coverage, 8 UEs | 410 → 790 m (1.93×) | 480 → 759 m (1.58×) | magnitude ✓ |
| coverage ratio | ≈$N_{UE}^{1/3}=2.0$ | 1.58 | 21% low |
| plateau throughput 30 m | ~150–200 / 85 / 50 / 28 | 180 / 90 / 45 / 23 Mbps | ✓ |
| cell edge 1500 m | ~0.013 Mbps | 0.0165 Mbps | ✓ |

The 8.05× depends on no fitted parameter, only on the three modelling decisions above (total
power, TDD duty 0.2, −3 dB parabola). It is the strongest evidence of this round.

**One structural difference.** In this model the JPTA curves merge at long distance; in the
paper they stay apart. Cause: concentrating the power in 4 RBs makes every JPTA user degenerate
to "4 RB × all slots" at the edge, independent of $N_{UE}$. The paper's curves staying apart
means its JPTA side does not concentrate power (RBs pinned at $\lfloor277/N_{UE}\rfloor$). Both
readings appear in the paper; the "select the RB count that maximises throughput" reading is
used here because only it produces the natural transition "equal near, JPTA wins far".

**A caveat that must travel with the 830%.** The JPTA side of this model assumes 28 dB in every
sub-band, which §5(a) shows to be unreachable for 4/8/16 users. So the system-level gain is
conditional on the beam design working, and the beam design does not work in this reproduction
under the stated hardware constraints. The honest statement is: with per-sub-band ideal beams
one gets 8.05×; with the beams the hardware can actually form (20.2 / 2.2 / −0.1 dB for
4/8/16 users) the JPTA advantage would be much smaller.

---

## 7 Gaps

| ID | Problem | Handling |
|---|---|---|
| **G1** | TDD uplink duty cycle not given while Fig. 9 is absolute Mbps | 0.2 (DDDSU); every absolute value scales with it, the ratios do not |
| **G2** | 10%-BLER SNR thresholds are in neither this paper nor ref. [7] | 3GPP CQI points interpolated in spectral efficiency, σ = 1.2 dB waterfall |
| **G3** | in-band shape of the JPTA gain not given ("at most 3 dB down") | parabola symmetric about the sub-band centre; affects the far-distance ratio only |
| **G4** | Fig. 7 not even in angle although a linear array must be | no point-by-point comparison; only the quantitative claim |
| **G5** | 3 dB claim needs 1.25 ns, Fig. 8 uses 2.5 ns | judged an internal inconsistency of the paper |
| **G6** | 28 dB per sub-band unreachable with one $(\varphi,\tau)$ pair | needs per-sub-band phase, impossible with one frequency-flat phase shifter per element |
| **G7** | eq. (4) has infinitely many optima ($\tau=1/\Delta f$ ≡ $\tau=0$) | delays taken from Fig. 8, eq. (4) used for the phase only |
| **G8** | PAA codebook not given | 16 beams uniform over ±55°, mismatch loss < 0.4 dB |

---

## 8 Pitfalls and bugs

| # | Symptom | Root cause |
|---|---|---|
| 1 | 2 users: both beams at −30°, in-band gain 4.5 dB | the delay was searched on the 2.5 ns grid while the required delays are sub-nanosecond; replaced by the FFT search over the whole period |
| 2 | every throughput zero | effective SNR stored as `Inf` and then offset by the path loss stays `Inf`; stable EESM (factor the `min` term out) |
| 3 | 16-user PAA gain −5.7 dB, oscillating with angle | `numel(p.Npaa)` returns 1 for a scalar, so the codebook had one beam |
| 4 | JPTA worse than PAA at short range | the sub-band parabola was written against the carrier instead of the sub-band centre, costing 20 dB in the upper sub-band |
| 5 | eq. (4) "optimum" at $\tau\approx8333$ ns, beam at 0° | the search covered $[0,1/\Delta f)$, whose endpoint is equivalent to $\tau=0$ (G7); the delays are taken from Fig. 8 instead |

---

## 9 Relation to the companion papers

| Dimension | this paper | 3-D beamforming (2024) | uplink coverage (2022) |
|---|---|---|---|
| array | 16V×16H, one TTD per column | 16×24, one TTD per element | 8V×16H, 16 sub-arrays |
| design variables | 16 pairs $(\varphi,\tau)$ | 384 pairs | 16 pairs, closed-form rainbow |
| solver | eq. (4) decouples; phase closed form, delay 1-D | LS / ∞-norm / greedy / Adam | closed form, no search |
| objective | 2-D type-1/type-2 beams + system level | 3-D sub-band gain | uplink coverage and throughput |
| carrier / bandwidth | 28 GHz / 400 MHz | 28 GHz / 95 MHz | 140 GHz / 2 GHz (and 28 GHz / 800 MHz) |
| BS gain | 28 dB | 25.84 dB | 24 dB |

The EESM parameter and the 29-row MCS table used here are Table 3 of the 2022 uplink paper; that
repository ships the same function as `jp_mcs_table.m`.
