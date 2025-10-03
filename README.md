# Germline SV Benchmarking Pipeline

A Snakemake workflow for running multiple structural‐variant (SV) callers on BAM/CRAM data and benchmarking their results against a single truth set with **[truvari](https://github.com/spiralgenetics/truvari)**.

---

##  Features

* Parallel execution on a compute cluster or locally.
* Modular caller rules—easy to add new SV callers.
* Automatic organization of outputs (VCFs, benchmark metrics, truvari stats).
* R-based downstream analysis.

---

##  Quick Start

```bash
snakemake --configfile config/pedigree.yaml --profile default
```

*Jobs are submitted to the cluster using the `default` profile.*

---

##  Project Layout

```
project/
├── config/
│   ├── versions.yaml        
│   └── config.yaml          # Example run configuration
├── workflow/
│   ├── Snakefile            # Entry point
│   ├── profiles/default/    # Cluster profile 
│   ├── rules/
│   │   ├── callers/         # One .smk per SV caller
│   │   ├── truvari.smk
│   │   └── analysis.smk
│   ├── r-scripts/           # R scripts for analysis
│   └── scripts/             # Helper scripts
└── results/                 # Created at runtime

```

---

##  Configuration

Edit `config/config.yaml` to set:
Replace every /PATH/TO/... placeholder with the actual locations on your system.

* `sample_path`: Read-write directory containing the input BAMs (`[ID].bam`).
Symlinking your BAMs into this directory is recommended. 
* `output`: Directory where all pipeline results will be written.
* `CALLERS`: List of SV callers to run (names must match the rules in workflow/rules/callers/).
* `truth_set`: Full path to the truth VCF used for benchmarking.


---

##  Results

The pipeline creates:

* **benchmarks/** – Snakemake `--benchmark-extended` outputs
* **callers/** – VCFs and intermediate files for each caller/sample
* **truvari/** – Per-run comparison outputs and statistics

> Some R scripts may run better locally; see `workflow/r-scripts/` for details.

---

##  Modifying the Pipeline

### Adding a New Caller

1. Create `workflow/rules/callers/{caller}.smk`

   * Ensure a rule outputs:
     `config["output"]/callers/{caller}/{sample}.vcf.gz`
2. Add `{caller}` to the `CALLERS` list in `config/config.yaml`.

### Adding or Removing Samples

* Place new `[ID].bam` files in `sample_path`.
* Add or remove sample names in `config/config.yaml`.

### Changing the Truth Set

* Update `truth_set` in the YAML.
* Adjust the `samples` list and optionally set a new output directory.

---

##  Notes

* Placeholder 


