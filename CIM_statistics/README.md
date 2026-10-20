# CIM seed replicas: setup and use

These scripts measure the scatter in your CIM yields. From each save point of a main CIM run, they start several short runs that differ only in their random seed. The scatter between seeds gives each CIM yield an error bar.

## Folder layout

You make the `potentials/` and `runs/` folders and copy in the files. The scripts create everything marked *(made)*.

```
CIM_replicas/
├── potentials/                  you add: WCH.tersoff WHe.table CHe.table HeHe.table w_eam4.fs
├── runs/
│   ├── WC_basal_C_50eV/         you add: wche0.restart he_cycle_500.restart he_cycle_1000.restart
│   │   ├── seed_b0000_s01/      (made) one folder per seed: cycles 1-100 from wche0.restart
│   │   ├── seed_b0000_s02/ ...
│   │   ├── seed_b0500_s01/ ...  (made) cycles 501-600 from he_cycle_500.restart
│   │   ├── replica_yields.csv   (made) one row per seed
│   │   ├── replica_summary.csv  (made) one row per save point: mean, SEM, 95 % CI, ...
│   │   └── replica_yields.png   (made) if matplotlib is available
│   ├── W001_100eV/              you add: w0.restart (+ he_cycle_500.restart, ... when they exist)
│   └── all_replica_summary.csv  (made) every run's summary in one table
├── in_CIM_rep.step2             WC input file
├── in_W_CIM_rep.step2           W input file
├── setup_replicas.sh            makes the seed folders
├── submit_all.sh                submits the waiting seeds (packed into as many jobs as you choose), then the summary job
├── run_replicas.slurm           the job script: runs the seeds of one job, one after another
├── summarize.slurm              runs summarize_replicas.py
├── summarize_replicas.py        writes the CSV files
└── logs/                        (made) SLURM output and the task lists (which seeds each job ran)
```

**Folder names:** a run folder must end in `_<E>eV`, because the He energy is read from the name. The material comes from the restart inside it: `wche0.restart` means WC, `w0.restart` means W.

**Save points:** every `he_cycle_<N>.restart` you add becomes a save point. The defaults are 0, 500 and 1000; a missing one is skipped with a message.

## Once: edit the job scripts

- In `run_replicas.slurm`, edit the four `EDIT` lines to match your main CIM jobs:
  - `module load lammps`
  - `LMP=`
  - `--ntasks`
  - `-t`
- In `summarize.slurm`, edit the `module load` line so you get a python3 with numpy. matplotlib is optional; without it you get the CSVs but no plot.
- If the files passed through Windows, fix the line endings once:
  ```
  dos2unix *.sh *.slurm *.step2 *.py
  chmod +x *.sh
  ```

## Run

```bash
cd CIM_replicas
./setup_replicas.sh          # 10 seeds x 100 cycles at save points 0/500/1000, every folder in runs/
./submit_all.sh -t           # test: submits ONLY the first waiting seed
```
When the test finishes, check that its `cycle_time.dat` reads 20.000 ps for every cycle. Note how long it took per cycle, then submit the rest:

```bash
./submit_all.sh -j 5         # all remaining seeds packed into 5 jobs, then the summary job
```

`setup_replicas.sh` ends by listing how many seeds are waiting in each run folder. It doesn't matter which run folders you named in that call.

### Packing seeds into jobs (`submit_all.sh`)

Your queue counts every array task as one submitted job, so 300 seeds as 300 tasks will hit a submission limit. Use `-j` to set how many jobs to use. Each job runs its seeds one after another, so set the wall time to cover them all.

| Command | What it submits |
|---|---|
| `./submit_all.sh -t` | the first waiting seed only, no summary job (a test) |
| `./submit_all.sh -j 1` | every waiting seed in one job, run one after another |
| `./submit_all.sh -j 5` | 5 jobs, seeds shared out evenly (the script prints how many each) |
| `./submit_all.sh -j 5 -r 2` | the same, but at most 2 of the 5 run at a time |
| `./submit_all.sh -j 5 -T 72:00:00` | the same, with a 72 h wall time per job |
| `./submit_all.sh -j 4 runs/WC_prism_100eV` | only the seeds of that run folder |
| `./submit_all.sh -n ...` | no summary job (one submission fewer) |
| `./submit_all.sh -d -j 5` | dry run: shows the split and the sbatch line, submits nothing |

Total submissions = the number of jobs, plus 1 for the summary unless you use `-n`. To find your limit: `sacctmgr show assoc user=$USER format=user,maxsubmit`, or `sacctmgr show qos format=name,maxsubmitpu`.

**Wall time:** a job that runs K seeds needs about K × (time per seed). Either set `#SBATCH -t` in `run_replicas.slurm` or pass `-T`. If a job runs out of time, the seed it was on stops part-way, and the seeds after it in that job never start. To recover:
- **The interrupted seed:** delete its folder, then rerun `./setup_replicas.sh` to recreate it.
- **The seeds that never started:** delete their `submitted` files.

Then run `./submit_all.sh` again. `logs/tasks_<date>.txt` shows which seeds each job had.

**Never submitted twice:** `submit_all.sh` writes a `submitted` file into every seed folder it submits, and skips any seed that has one or has a `log.lammps`. It always looks over all of `runs/` (or the folders you name), so you can run it again at any time.

### `setup_replicas.sh` options

| Option | Default | Meaning |
|---|---|---|
| `-n` | 10 | seeds per save point (1–99) |
| `-c` | 100 | cycles per seed |
| `-b` | `"0 500 1000"` | save points to use |
| `-d` | 50 | keep `dump.bulk/full` every N cycles |
| `-x` | – | extra LAMMPS arguments added to every seed, e.g. `-x "-var Tsub 739"`. Use only to match a main run that used non-default values. |
| `-f` | – | rebuild seed folders that already exist (overwrites their restart; only before they've run) |
| folders | all of `runs/` | e.g. `./setup_replicas.sh -n 20 runs/WC_basal_W_100eV` |

Seed folders that already exist are never touched. To add seeds, rerun with a larger `-n`, then run `./submit_all.sh` again; it submits only the new ones.

### Seeds and random numbers

The seed number passed to LAMMPS is `rep = save point + s`. For example, `seed_b0500_s03` is rep 503. The main run is rep 0.

LAMMPS uses a single random stream per run, seeded by the first `random()` call (I checked this in the 29Sep2021 source). The input files add `rep` to that first seed. Seeds started from the same save point therefore place their He differently from the first impact, and rerunning a seed folder repeats it exactly.

## Results

The summary job runs automatically after the array. You can also run it yourself at any time, for example while seeds are still running:

```bash
python3 summarize_replicas.py                       # all of runs/
python3 summarize_replicas.py runs/W001_100eV       # one run
python3 summarize_replicas.py --include-incomplete  # also use unfinished seeds
```

The statistics use the C yield for WC and the W yield for pure W. W sputtered from WC is only a few atoms per seed, so it is reported pooled with a 95 % upper limit. Only complete seeds count unless you pass `--include-incomplete`. The script prints a warning if any cycle is not 20 ps.

`replica_summary.csv` columns:

| Column | Meaning |
|---|---|
| `branch` | save point (0, 500, 1000) |
| `n_seeds`, `n_incomplete` | seeds used / seeds not yet finished |
| `window_centre_cycle`, `fluence_centre_m2` | where the window sits. Impact areas: basal 3567, prism 3420, W 3490 Å². |
| `Y_mean`, `Y_SD`, `Y_SEM`, `Y_CI95` | yield (atoms per He) and its spread; the CI uses Student t |
| `dispersion` | var/mean of the counts: 1 means Poisson, above 1 means bursty sputtering |
| `seeds_for_10pct_SEM` | how many seeds a 10 % SEM would need; add seeds where this exceeds `n_seeds` |
| `W_total`, `impacts_total`, `Y_W_pooled`, `Y_W_95UL` | pooled W yield and its 95 % upper limit |

`replica_yields.csv` has one row per seed: cycles done, W and C sputtered, yields, whether the seed is complete, and the cycle-length check.

## How many seeds

Starting estimates from your 10/04 main runs; `seeds_for_10pct_SEM` refines them after the first round.

| Surface | 50 / 100 / 200 eV |
|---|---|
| basal, C-terminated | 17 / 8 / 13 |
| basal, W-terminated | 51 / 25 / 18 (use `-n 20`, or `-c 200`) |
| prism | 13 / 6 / 11 |

**Pure W:** W(001) at 100 eV sputtered one atom in 424 cycles, so seeds there only give an upper limit. W seeds become useful at 200 eV. Until the main run passes cycle 500, only save point 0 is available.

## Tested here

Tested with LAMMPS 2022 on small slabs, using a stand-in W potential:
- **Setup:** folders are created correctly. Re-running only adds new seeds, and missing save points are skipped.
- **Seeds:** the job script ran the W seeds; the same rep gives the same He position and different reps give different ones.
- **Save points:** a seed started from `he_cycle_500.restart` continues at cycle 501.
- **Summary:** it produced the CSVs and the plot from those seeds plus dummy WC outputs.
- **Submitting:** with a stand-in sbatch, `setup_replicas.sh` listed the waiting seeds of two run folders, `-t` ran one seed, and `-j 2` ran the other 3 W seeds as 2 jobs. Seeds that had already been submitted or started were not picked up again.

The WC input was tested the same way on 5 Oct. The real potentials and your cluster build weren't run here, so time one seed first.
