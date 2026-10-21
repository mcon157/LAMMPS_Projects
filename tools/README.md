# tools

`build_notebooks.py` turns the analysis scripts into the notebooks next to them. It then runs the notebooks, so that their figures, renders and GIFs are stored inside them.

```bash
pip install nbformat nbclient ipykernel pillow
python tools/build_notebooks.py                       # all notebooks
python tools/build_notebooks.py rand_schematics gifs  # selected ones
python tools/build_notebooks.py --no-exec             # code only, no outputs
python tools/build_notebooks.py --data paths.json     # data-driven notebooks on your own runs
```

`paths.json` maps a notebook to the settings used while it runs, for example:

```json
{"ovito_insertion_movies": {"NB_RUN": "/path/to/basal_C/200eV"},
 "plot_erosion_energies":  {"NB_CSV": "/path/to/figures/csv"}}
```

Only the outputs come from that data. The saved settings cell keeps its generic values, and build-machine paths are removed from the printed output.

How a notebook is built:
1. **Split.** The script is cut into cells: imports, the config block, sections (at banner comments) and the main call. `sys.exit(main())` becomes `main()`.
2. **Settings cell.** A cell at the top sets `ARGS`, which stands in for the command line. Its variables start with `NB_`, so the script's own globals cannot overwrite them.
3. **Execute.** Each notebook runs in a scratch copy of the repository.
   - Self-contained scripts run completely (`execute="all"`).
   - Scripts that need cluster data or OVITO keep their code unrun and only display example outputs (`execute="show"`).
4. **Embed.** Outputs are embedded as PNG; GIF bytes are stored under the same mimetype and still animate.
   - PNGs over 700 kB are scaled to 1200 px wide.
   - GIFs over the limit (1.5 MB by default) are first scaled down, and lose frames only if that is not enough.

To add a notebook, add an entry to `NOTEBOOKS` at the top of the script. Always edit the `.py`, then rebuild.
