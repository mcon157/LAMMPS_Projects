# Random surface orientation, single-impact model (SIM)

Following Jussila, Granberg & Nordlund, *Nucl. Mater. Energy* 17 (2018) 113: rotate a crystal by random Euler angles, cut a hemisphere whose flat face is the irradiated surface, and average the yield over many orientations. The result stands for a polycrystal with random grain orientations.

## Files

| File | What it does |
|---|---|
| `in_rand_W.step1`, `in_rand_WC.step1` | Build one surface `s`: crystal ball → rotate Rz(γ), Rx(β), Rz(α) → keep r ≤ R, z ≤ 0 → fixed / thermostat / NVE shells → minimise + 10 ps at T → restart. |
| `in_rand_W.step2`, `in_rand_WC.step2` | `nimp` single impacts on that surface (each re-reads the restart). One line per impact. |
| `run_rand.slurm` | Array job: task ID = surface ID; runs step 1 then step 2 at 50/100/200 eV. |
| `analyze_rand.py` | Per-surface yields, mean ± SEM over surfaces (Jussila's formula), the "clean" mean without fixed-shell hits, and yield on an inverse pole figure. |
| `rand_schematics.py` | The schematics in `figs/`. Rebuilds the geometry with the decks' own maths and random stream; atom counts match LAMMPS exactly. |
| `rand_surface_views.py` | Regular vs random surface, side by side: the unrotated crystal with the plane surface *s* cuts through it, plus top views of W(001) / basal WC and of surface *s* (`--s 7`). |
| `rand_orientation_gif.py` | Animated GIF, one frame per surface drawn in the unrotated lattice, Miller label top-left. Frame range, ms per frame, reference frames and resolution are set in the OPTIONS block at the top (or `--mat --first --last --ms`). |
| `rand_geom.py` | Shared geometry for the two scripts above: lattices, rotation, hemisphere cut, Miller labels (cubic and 4-index hexagonal). |
| `rand_angles_lammps.csv` | Euler angles and origin shifts for surfaces 1–500, pre-computed from the LAMMPS random stream so the figures show the surfaces the runs use, even without the LAMMPS Python module. |

All decks are for LAMMPS 29Oct2020. They have no `&` continuations and no if/else structure; `if` is used only for exit conditions, as in your CIM decks. Run them with `-in`.

## How the random orientation is drawn

- **Euler angles (ZXZ).** α and γ are uniform on 0–360°. β = arccos(2u − 1), so cos β is uniform; this gives a uniform distribution over all orientations, as in Jussila.
- **Random lattice origin.** The origin shift (o_x, o_y, o_z) is random in [0, 1), so the z = 0 cut samples every termination and step position. This matters for WC, where the top layer can be C or W.
- **Logged per surface** (`rand_<MAT>_surfaces.dat`): the Euler angles, the surface normal in the crystal frame n = Rᵀẑ = (sin β sin γ, sin β cos γ, cos β), the origin shift, the atom counts, z_s and (WC only) the top-layer C fraction.
- **Low-index references** (random azimuth only, Jussila's 30-surface sets) are selected by overriding the angles:
  - W: (001) `-var beta 0 -var gamma 0`; (011) `-var beta 45 -var gamma 0`; (111) `-var beta 54.7356 -var gamma 45`
  - WC: basal `-var beta 0 -var gamma 0`; prism (10-10) `-var beta 90 -var gamma 0`; a-plane (11-20) `-var beta 90 -var gamma 90`
- **Random-number note.** LAMMPS keeps one random stream per run, seeded by the first `random()` call; later seed arguments are ignored. This was checked in the 29Oct2020 source. The decks therefore put the surface ID (and energy) into that first seed. Re-running a given `s` rebuilds exactly the same surface.

## What counts as sputtered

The decks record three counts, so the results can be compared with both Jussila and our CIM decks:

- `cap`: crossed the capture plane z_s + 18 Å. Counted and deleted on the spot, as in the WC CIM deck.
- `end`: above z_s + d_cut when the impact ends. W: 5.5 Å, the EAM cutoff, which is Jussila's criterion. WC: 3.8 Å for W and 3.0 Å for C, the Tersoff R+D values used in your SIM deck.
- `lost`: left the box between polls.

Sputtered = cap + end + lost.

An impact runs until the He has thermalized **or left**, plus `t_settle` = 2 ps, with a cap of 10 ps. Unlike the old SIM deck, a reflected He no longer stops the count, because atoms it set moving can still leave.

`hit_fixed` = 1 marks impacts whose He reached the fixed shell. Jussila reports these as the lower error bar; `analyze_rand.py` gives the "clean" mean without them.

## Step 1: check W against Jussila first

Jussila's He energies are 200/300/500/1000 eV, so only 200 eV overlaps our set, and at 200 eV W sputtering is rare (about 1 W atom per 1000 impacts). A useful check therefore needs 500 eV as well.

- **Conditions:** 50 random surfaces, θ = 30°, 1000 impacts each, at 200 and 500 eV.
- **Jussila set-up:** `-var T 300 -var R 35 -var a <a at 300 K for your EAM4>`.
- **Our production set-up:** the same at 0.2 Tm, R = 50.
- **Also run 50 and 100 eV:** both should give zero W (below threshold), as a sanity check.

Targets, random-surface average at θ = 30°, read by eye from Jussila Fig. 6 (about ±20%):

| E (eV) | Y_W, random surfaces |
|---|---|
| 200 | ≈ 1.2 × 10⁻³ |
| 300 | ≈ 7 × 10⁻³ |
| 500 | ≈ 2 × 10⁻² |
| 1000 | ≈ 4 × 10⁻² |

Known set-up differences to state when comparing:

| | Jussila 2018 | This deck |
|---|---|---|
| W–W potential | Marinica EAM + Sand short-range modification | EAM4 (`w_eam4.fs`): check which version you have |
| He–W interaction | pure ZBL | Juslin–Wirth table |
| Electronic stopping | included (friction, KE ≥ 10 eV) | not included |
| Thermostat shell | Berendsen | Nosé–Hoover |
| Code | PARCAS | LAMMPS |

## Step 2: WC at 50/100/200 eV

**Pilot:** about 10 surfaces × 200 impacts at 100 eV. Measure the surface-to-surface spread σ_s and the time per impact. The standard error of the mean over N_s surfaces with N_i impacts each is

SEM² ≈ σ_s²/N_s + D·Ȳ/(N_i·N_s)

where D ≈ 1–2.5 is the burstiness of C sputtering measured from your CIM runs.

- **C yield ~0.2:** counting noise is small at 500–1000 impacts, so N_s sets the precision. With σ_s/Ȳ ≈ 0.3, 40 surfaces give about 5%.
- **W yield:** expect a handful of events per 1000 impacts. Report it pooled, with a Poisson upper limit.

**Reference sets:** run the basal, prism and a-plane references with the same deck. They double as low-index SIM points that use exactly the same sputtering criterion as CIM, which matters for the CIM-vs-SIM comparison.

**Covariates:** plot yield against the surface normal (IPF) and the top-layer C fraction. The question is whether orientation matters or termination does.

## Tested here (LAMMPS 23Jun2022 Python build, stand-in W potential, small R)

- **Builder:** runs; atom counts agree with an independent numpy rebuild (W 1556 / C 1554 at R = 25 Å).
- **Normal formula:** checked against LAMMPS `displace_atoms` rotations.
- **Impact loop:** runs. Capture, end and lost counts, He fate and the `hit_fixed` flag all work. A C atom captured at 300 eV, 75° was recorded with its kinetic energy.

The production potentials and the 29Oct2020 build were **not** run here; do a 20-impact pilot on the cluster first.
