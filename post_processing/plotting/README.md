# plotting: cross-run erosion figures

These scripts compare the `erosion.csv` files written by [`../ovito_insertion_movies.py`](../ovito_insertion_movies.py) across He energies and surfaces. They need numpy, pandas and matplotlib.

| Script | Overlays | Main outputs |
|---|---|---|
| `plot_erosion_energies.py` / `.ipynb` | all energies, one figure per surface | `erosion_regression_<plane>.png`, `erosion_sputtering_<plane>.png` (or twin-axis `erosion_combined_<plane>.png`) |
| `plot_erosion_compare.py` / `.ipynb` | all surfaces at one energy | `surface_W_pct_<E>eV.png`, `he_retention_<E>eV.png`, `he_clusters_<E>eV.png`, `W001_surface_He.png`, `erosion_summary_equal_fluence.csv` |
| `make_nov2_figures.py` | runs both for a whole update folder | `figures/`, plus slide-sized `figures/slide2/` and `figures/slide3/` |

## Input names

The scripts expect one CSV per run:
- `erosion_WC_basal_C_<E>eV.csv`
- `erosion_WC_basal_W_<E>eV.csv`
- `erosion_WC_prism_<E>eV.csv`
- `erosion_W001_<E>eV.csv`

Plane keys are `bC`, `bW`, `pr` and `W001`; the file names are set in `PLANE_FILES`. `make_nov2_figures.py` builds these names itself:
1. Run it from a folder holding `CIM/<surface>/<E>eV/erosion.csv`.
2. It copies each file to `figures/csv/` under the expected name.
3. For W(001) it adds `surface_height_A`.
4. It then calls both plotting scripts with shared axis limits, so the panels line up on slides.

## Run

```bash
python plot_erosion_energies.py --dir CSVs --out figs --planes bC bW pr --mode separate
python plot_erosion_energies.py --dir CSVs --trend on            # add fits and print their equations
python plot_erosion_compare.py  --dir CSVs --out figs --energy 100 --compare-cycle 1000
python plotting/make_nov2_figures.py                             # from the update folder
```

Each script has a settings block at the top: colours, axis limits, smoothing, which columns to plot. The most used settings can also be changed from the command line; `-h` lists them. The notebooks show the figures made from the WC basal (C- and W-terminated) and prism runs and the W(001) runs.
