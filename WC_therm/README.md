# WC_therm: WC lattice constants at temperature

`in.step1` runs NPT on bulk WC (12 × 12 × 12 orthogonal cells, x and y coupled). It averages a and c and computes the thermal expansion relative to a<sub>0</sub> = 2.92329 Å and c<sub>0</sub> = 2.81805 Å at 300 K.

Outputs:
- `avg_<T>.dat`: running averages;
- `thermal_expansion_summary.dat`: one line per temperature, appended; T, a, c, V<sub>cell</sub>, α<sub>a</sub>, α<sub>c</sub>, α<sub>V</sub>.

```bash
cp ../potentials/WCH.tersoff .
# set 'variable T equal ...' in in.step1 (equal-style, so -var cannot override it)
lmp -in in.step1
```

The slab decks use the result at 614.6 K (0.2 T<sub>m</sub>): a = 2.92949 Å, c = 2.82421 Å.
