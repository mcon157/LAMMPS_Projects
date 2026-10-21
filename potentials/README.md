# potentials

Copy the files a deck needs into its run folder; the decks read them by name.

| File | Style | Pairs | Source |
|---|---|---|---|
| `WCH.tersoff` | `tersoff/zbl` | W–W, W–C, C–C | Juslin et al., *J. Appl. Phys.* 98, 123520 (2005): analytical W–C–H bond-order potential with a ZBL core. He is mapped to `NULL`. The H entries are unused |
| `w_eam4.fs` | `eam/fs` | W–W (pure W decks) | "EAM4" of Marinica et al., *J. Phys.: Condens. Matter* 25, 395502 (2013) ([OpenKIM entry](https://openkim.org/id/EAM_Dynamo_MarinicaVentelonGilbert_2013EAM4__MO_046576227003_001)) |
| `WHe.table` | `table linear 10000`, keyword `WHe` | W–He | Juslin & Wirth, *J. Nucl. Mater.* 432, 61 (2013). The data are identical to LAMMPS `potentials/W_He_JW2013.table` |
| `CHe.table` | `table linear 10000`, keyword `CHe` | C–He | Group-internal table. See the region list below the table |
| `HeHe.table` | `table linear 10000`, keyword `HeHe` | He–He | Beck, *Mol. Phys.* 14, 311 (1968), with the short range of Morishita et al., *NIMB* 202, 76 (2003) |

`CHe.table` is built from three regions:
- ZBL for r < 0.8 Å;
- a quintic polynomial for 0.8–2.8 Å;
- a Lennard-Jones tail beyond that.

It has no independent validation, so it is the least certain of the five.

Typical `pair_coeff` blocks:

```
# WC decks (W 1, C 2, He 3)
pair_style  hybrid/overlay tersoff/zbl table linear 10000
pair_coeff  * * tersoff/zbl WCH.tersoff W C NULL
pair_coeff  1 3 table WHe.table WHe
pair_coeff  2 3 table CHe.table CHe
pair_coeff  3 3 table HeHe.table HeHe

# W decks (W 1, He 2)
pair_style  hybrid eam/fs table linear 10000
pair_coeff  * * eam/fs w_eam4.fs W NULL
pair_coeff  1 2 table WHe.table WHe
pair_coeff  2 2 table HeHe.table HeHe
```

Neither He–host table was fitted to WC: W–He was fitted for bcc W, and C–He comes from a separate group table. [`He_BE/model_schema.yaml`](../He_BE/model_schema.yaml) lists ranges, checks and sensitivities.
