# WC_SBE: WC surface binding energy and kickout threshold

`in.step2_kick` works on one WC slab (30 × 18 orthogonal cells) in three phases:

| Phase | What it does |
|---|---|
| 0 | build the slab, add vacuum, minimize |
| 1 (method A) | surface on NVE, interior thermostatted; average the PE/atom of the top layer → finite-T surface binding energy (`sbe_WC_<T>K.dat`) |
| 2 (method B) | re-equilibrate, find the top layer, kick one atom upward with `Epush`, follow it → escaped or not (`kickout_WC_T<T>_trials.dat`) |

```bash
cp ../potentials/{WCH.tersoff,WHe.table,CHe.table,HeHe.table} .
lmp -in in.step2_kick     # first pass (kickID = 0): writes top_atoms_before_kick_T<T>.dump only
```

Then pick a top atom from that dump; its type column says W or C. Set `kickID` and `mass_kick` in the deck and rerun, scanning `Epush` to find the threshold.

`T`, `Epush`, `kickID` and `mass_kick` are equal-style, so edit them in the deck rather than with `-var`. The deck's a, c (2.92329, 2.81805 Å) are the [`WC_therm`](../WC_therm) reference values at 300 K. For another T, use the `WC_therm` result at that T.
