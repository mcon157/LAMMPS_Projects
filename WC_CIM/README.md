# WC_CIM: WC continuous-impact model

He impacts WC one at a time, and the slab carries all damage, He and composition change forward from one cycle to the next. Each cycle lasts exactly 20 ps from one insertion to the next, so the flux is defined. He arrives at 30° from the normal, inside the central 2/3 × 2/3 of the cell.

| File | Purpose |
|---|---|
| `in_basal_092926.step1` | basal (0001) slab, C-terminated |
| `in_basal_W_092926.step1` | basal (0001) slab, W-terminated |
| `in_prism_092926.step1` | prism / m-plane (10-10) slab, C-terminated |
| `in_CIM_092926.step2` | the impact cycles (all surfaces) |
| `run.slurm` | self-chaining job for step 2 on a 3 h partition |

Each step 1 deck runs two stages and writes `wche0.restart` plus the cycle-0 baselines `dump.bulk-0` and `dump.full-0`:
- **Stage A.** NVT on the whole mobile slab.
- **Stage B.** The same NVE-top / NVT-bottom split as step 2.

All decks use T = 614.6 K (0.2 T<sub>m</sub>), a = 2.92949 Å and c = 2.82421 Å (from [`WC_therm`](../WC_therm)), and atom types W 1, C 2, He 3. They were written for LAMMPS 29Oct2020.

## Step 2 cycle

1. Read the restart.
2. Insert one He `standoff` Å above the current surface.
3. Follow the cascade until the He reflects, thermalizes, channels into the thermostat layer, or `t_impact` runs out.
4. Relax until the cycle has lasted `t_cycle`.
5. Record sputtering: atoms that pass the capture plane are counted and deleted.
6. Write the dumps and the restart.

| Output | Contents |
|---|---|
| `sputtering.dat` | every capture: cycle, step, time, n<sub>W</sub>, n<sub>C</sub>, kinetic energies |
| `sputter_summary.dat` | net W, C loss and He count per cycle |
| `cycle_time.dat` | cycle length check (must read 20 ps) |
| `dump.bulk-<i>.lammps`, `dump.full-<i>.lammps` | W + C, and W + C + He, every cycle |
| `he_cycle.restart`, `last_cycle.txt` | carry-over state and resume anchor |
| `he_cycle_<N>.restart` | permanent save points every `bak_every` (500) cycles, which seed [`CIM_statistics`](../CIM_statistics) |

| Knob (`-var`) | Default | Meaning |
|---|---|---|
| `e` | 100 | He energy (eV) |
| `istart`, `nimplants` | 1, 10000 | first and last cycle |
| `T`, `ang`, `t_cycle`, `t_impact` | 614.6 K, 30°, 20 ps, 10 ps | |
| `standoff`, `z_channel` | 8 Å, 14 Å | insertion height; depth counted as channelled |
| `he_dump_every`, `bak_every` | 50, 500 | He trajectory dump and save-point intervals |

## Run

One folder per surface and energy, e.g. `basal_C/100eV/`:

```bash
cp ../../potentials/{WCH.tersoff,WHe.table,CHe.table,HeHe.table} .
cp ../../WC_CIM/in_basal_092926.step1 ../../WC_CIM/run.slurm .
cp ../../WC_CIM/in_CIM_092926.step2 in.step2     # the deck jumps to in.step2; set 'variable e index 100'
lmp -in in_basal_092926.step1                    # once
sbatch run.slurm                                 # chains itself until NIMP cycles
```

`run.slurm` handles the chaining:
- **Each job.** It reads `last_cycle.txt`, queues its successor (`afterany`, so the chain survives the wall-time kill), and continues from the next cycle.
- **Stopping.** The chain stops when `NIMP` (1000) is reached or a needed restart is missing.
- **Before the first run.** Edit the `EDIT` lines (account, module, `lmp` binary).
- **Fresh start.** Run `rm -f last_cycle.txt he_cycle.restart sputtering.dat sputter_summary.dat cycle_time.dat`.

Analysis: [`post_processing/ovito_insertion_movies.py`](../post_processing) and `he_cluster_movies.py`.
