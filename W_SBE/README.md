# W_SBE: W(001) surface binding energy

These decks benchmark the pure-W potential (`w_eam4.fs`, see [`potentials`](../potentials)) for sputtering thresholds.

| Deck | What it does | Appends to |
|---|---|---|
| `in.step1` | bulk bcc W, NPT at T: average a and thermal expansion vs a<sub>0</sub> = 3.152 Å at 300 K | `thermal_expansion_W_summary.dat` (+ `W_bulk_relaxed_<T>K.restart`) |
| `in.step2_pe` | W(001) slab at 0 K, minimized: SBE = −⟨PE of a surface atom⟩ | `sbe_PE_0K.dat` |
| `in.step2_kick` | W(001) slab at 0 K: kick the central surface atom upward with `Epush` and record whether it escapes | `kickout_0K_trials.dat` |

```bash
cp ../potentials/w_eam4.fs .
lmp -in in.step1          # set T in the deck first
lmp -in in.step2_pe
lmp -in in.step2_kick     # repeat for a range of Epush to find the kickout threshold
```

`T` and `Epush` are equal-style variables, so `-var` cannot override them. Edit the deck instead.

The slabs use a = 3.165 Å, as does `SIM_random/in_rand_W.step1`; 3.152 Å reproduces Jussila / Marinica at 300 K.
