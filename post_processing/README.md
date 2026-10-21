# post_processing: CIM analysis and renders

These scripts read the per-cycle dumps written by the CIM step 2 decks ([`WC_CIM`](../WC_CIM), [`WC_electronic_stopping`](../WC_electronic_stopping)). From them they produce:
- surface metrics, written to `erosion.csv`;
- composition and He profiles;
- He cluster statistics;
- OVITO stills and GIFs.

The cross-run comparison plots are in [`plotting/`](plotting).

| File | Purpose |
|---|---|
| `ovito_insertion_movies.py` / `.ipynb` | surface metrics per cycle, erosion / composition / He / order profiles, OVITO renders |
| `he_cluster_movies.py` / `.ipynb` | He cluster detection, tracking, depth profiles, He-only renders |
| `erosion_metrics.pdf` | definitions and equations of the surface metrics |
| `plotting/` | all energies or surfaces in one figure, from the `erosion.csv` files |

## ovito_insertion_movies.py

Input files, in the run folder:

| File | Holds | Used for |
|---|---|---|
| `dump.bulk-<i>.lammps` | W + C | bulk views |
| `dump.full-<i>.lammps` | W + C + He | He views and profiles |
| `dump.bulk-0.lammps` | the cycle-0 baseline from step 1 | reference (cycle 1 is used if it is missing) |

```bash
export QT_QPA_PLATFORM=offscreen                       # headless
python ovito_insertion_movies.py RUN_DIR               # everything
python ovito_insertion_movies.py RUN_DIR --fast        # metrics + plots, no renders
python ovito_insertion_movies.py RUN_DIR --surface prism --zoom 2.0
```

Pass the run folder. Without one, the script uses its own folder.

`--surface` takes `basal` (C-terminated), `basal_W` or `prism`. Everything else is set in the CONFIG block at the top, including:
- `HE_ENERGY_EV` (used in plot titles);
- `T_CYCLE_PS` (20);
- `IMPACT_FRAC` (2/3);
- `ANALYSIS_STRIDE`;
- which plots and renders to make (`JOBS`, `PLANES_*`).

Atom types follow the WC decks (W 1, C 2, He 3).

Surface metrics written to `erosion.csv` (equations in `erosion_metrics.pdf`):

| Metric | Definition |
|---|---|
| `z_surf` | Gibbs dividing surface of the W + C density, averaged over one chemical period |
| `z_cryst` | top of the ordered lattice, from the layering order parameter S = \|⟨exp(2πiz/d)⟩\| |
| growth | z_surf − z_cryst (≥ 0) |
| regression (front) | z_cryst(i) − z_cryst(0) (≤ 0) |
| regression (deficit) | −[D(i) − D(0)], D = integrated order deficit (sensitive to light damage) |

The layer period d is detected from the baseline: basal c/2 = 1.412 Å, prism a√3/6 = 0.846 Å.

Other `erosion.csv` columns, which feed `plotting/`:
- `sputtered_W`, `sputtered_C`: cumulative sputtering;
- `surface_W_at_pct`: near-surface W at.%;
- `he_retained`, `he_isolated`, `he_clusters_2_3`, `he_clusters_4_6`, `he_clusters_7plus`: He counts.

Outputs (in the run folder):
- `erosion.csv`, `erosion_combined.png`, `composition_profile.png`, `he_profile.png`, `stoichiometry_profile.png`, `order_profile.png`;
- `depth_slice.png`;
- stills and GIFs: `bulk_*`, `full_*`, `planes_*`.

Expect large GIFs for long runs. Raise `ANALYSIS_STRIDE`, or use `--fast` while checking.

## he_cluster_movies.py

How clusters are found and tracked:
- **Input.** He from `dump.full-<i>.lammps`. The host atoms are used only to locate the surface, so every He gets a depth.
- **Cluster rule.** Two He are in one cluster when closer than 2.5 Å (periodic in x and y; the same rule as above).
- **Tracking.** Clusters are followed between frames by shared atom IDs.
- **He type.** Detected automatically (3 in WC decks, 2 in W decks), or set with `--he-type`.

```bash
python he_cluster_movies.py --dir RUN_DIR --jobs 4             # impact area guessed from the path (prism/basal/w001)
python he_cluster_movies.py --dir RUN_DIR --area 3420 --stride 5 --no-gif
#   --color-by size|depth|track   --draw atoms|centroids   --hide-isolated   --gif-frames N
```

Outputs, in `RUN_DIR/he_clusters/`:
- CSVs: `he_clusters_per_cycle.csv`, `he_cluster_list.csv`, `he_cluster_tracks.csv`;
- figures: `he_cluster_evolution.png`, `he_depth_profile.png`, `he_clusters_final.png`;
- GIF: `he_clusters.gif`.

## Notebooks

The notebooks contain the same code. Their embedded outputs come from a WC basal (0001), C-terminated, 200 eV run. The full-cell GIFs are too large to embed, so only final-frame stills are shown.
