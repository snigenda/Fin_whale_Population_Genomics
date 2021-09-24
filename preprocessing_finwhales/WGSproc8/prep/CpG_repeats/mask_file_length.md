# Different repeat masks

`DIRECTORY="/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/"`

| Maskfile                                                     | description                          | WGSproc8                          | Neutral Region | Canonical Coding Region | Length        |
| ------------------------------------------------------------ | ------------------------------------ | --------------------------------- | -------------- | ----------------------- | ------------- |
| `CpG_repeats_all.bed`                                        | `WM` output + `RM` output + UCSC CpG | Yes                               | Default        | Default                 | 1,247,900,490 |
| `GCF_000493695.1_BalAcu1.0_genomic_repeats.bed`              | `WM` output                          | included in `CpG_repeats_all.bed` | Default        | Default                 | 741,334,011   |
| `GCF_000493695.1_BalAcu1.0_rm.out.bed`                       | `RM` output                          | included in `CpG_repeats_all.bed` | Default        | Default                 | 969,461,823   |
| Overlap `GCF_000493695.1_BalAcu1.0_genomic_repeats.bed` and `GCF_000493695.1_BalAcu1.0_rm.out.bed` |                                      |                                   |                |                         | 487,884,195   |
| Merge `GCF_000493695.1_BalAcu1.0_genomic_repeats.bed` and `GCF_000493695.1_BalAcu1.0_rm.out.bed` |                                      |                                   |                |                         | 1,221,985,945 |

