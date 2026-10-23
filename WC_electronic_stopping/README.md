# WC_electronic_stopping: CIM with composition-aware electronic stopping

This is the WC CIM step 2 of [`WC_CIM`](../WC_CIM) with `fix electron/stopping` added; the stopping blocks are tagged `### ES` in the deck. Every cycle the slab is split at z<sub>surf</sub> − `d_alt` into two regions:
- the altered surface layer;
- the deep slab.

The C/W ratio is counted in each region. `elstop.py` then writes one stopping table per region by Bragg's rule:

S<sub>p</sub>(E) = n<sub>W</sub> σ<sub>p@W</sub>(E) + n<sub>C</sub> σ<sub>p@C</sub>(E),  with n<sub>W</sub> = `nW_ref` and n<sub>C</sub> = (C/W) · `nW_ref`

The W sublattice is taken as intact, since C sputtering leaves C vacancies. Per cycle, the C/W ratios and the energy removed by stopping go to `elstop.dat`.

| File | Purpose |
|---|---|
| `in_CIM_ESdyn.step2` | the CIM step 2 deck with dynamic stopping tables |
| `elstop.py` / `.ipynb` | builds the per-element basis once, and the stopping tables from it (pure Python for `table`/`grid`) |
| `basis_rustbca.dat` | σ(E) for W, C, He in W and in C, 10–1000 eV, from the RustBCA electronic-stopping model |
| `elstop_WC.table` | an example fixed-composition table (stoichiometric WC) for a static `fix electron/stopping` |
| `build_lammps.sh` | builds a current stable LAMMPS with MANYBODY + EXTRA-FIX (+ OPENMP) as `lmp_elstop` |
| `run.slurm` | the self-chaining 3 h job (same logic as `WC_CIM/run.slurm`) |

## Build LAMMPS (once)

`fix electron/stopping` is not in the 29Oct2020 build used for the other decks.

```bash
# edit the module block first: CMake >= 3.20, a C++17 compiler, MPI
bash build_lammps.sh          # -> $HOME/software/lammps-install/bin/lmp_elstop, then checks the needed styles
```

## Run

The run folder needs these files:
- the potentials (`WCH.tersoff`, `WHe.table`, `CHe.table`, `HeHe.table`);
- `wche0.restart` from a `WC_CIM` step 1;
- `in_CIM_ESdyn.step2`, `elstop.py` and `basis_rustbca.dat`;
- `run.slurm`.

```bash
lmp -in in_CIM_ESdyn.step2 -var e 100                          # interactive
lmp -in in_CIM_ESdyn.step2 -var e 100 -var istart 501          # resume
lmp -in in_CIM_ESdyn.step2 -var e 100 -var ecut 5.0            # threshold test
sbatch run.slurm                                               # chained jobs (EDIT account and LMP path)
```

The deck calls `python3 elstop.py` through LAMMPS `shell` every cycle, so python3 must be available on the compute nodes. If `shell` is not allowed there, make a family of tables over C/W beforehand with `elstop.py grid`.

| Knob | Default | Meaning |
|---|---|---|
| `ecut` | 10 eV | no drag below this kinetic energy |
| `elstop_basis` | `basis_rustbca.dat` | the σ(E) basis |
| `d_alt` | 15 Å | altered-layer thickness below z<sub>surf</sub> |
| `nW_ref` | 0.04806 Å⁻³ | W density in WC (bcc W: 0.06306) |
| `pybin` | `python3` | python on the compute node |

## elstop.py

```bash
python3 elstop.py basis --backend rustbca -o basis_rustbca.dat                     # RustBCA model (no external code)
python3 elstop.py basis --backend srimfile --srim-file He@W "He in W.txt" ... -o basis_srim.dat   # SRIM SR output
python3 elstop.py table --basis basis_rustbca.dat --types W C He --target W:0.5,C:0.5 --rho 15.63 -o elstop_WC.table
```

Basis backends:
- **`srim`.** Runs SRIM's SR module through pysrim. Needs SRIM, and wine on Linux.
- **`srimfile`.** Reads SR tables you saved yourself.
- **`rustbca`.** Uses RustBCA's interpolated Lindhard–Scharff / Bethe model.

To use a SRIM basis in the deck, pass `-var elstop_basis basis_srim.dat`.
