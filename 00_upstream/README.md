# 00_upstream — alignment and variant calling (how the per-family WDL was run)

Everything downstream in this repository (`01_qc` → `13_phase_dnm`) starts from the outputs of PacBio's
**HiFi-human-WGS-WDL** run once per family. This module records exactly how that was done, so the callsets are
reproducible from unaligned HiFi reads: the workflow version pin, the container images it resolves to, our inputs format,
the Slurm/Singularity execution layer, the cohort-wide joint call, and the two operational lessons that matter for
anyone re-running it. Nothing here is a new method; it is the configuration and the operational wrapper around a
published, versioned workflow. The separate per-family clinical-interpretation workflow is not part of this repository.

## Workflow and versions

| item | value |
|---|---|
| workflow | [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) `workflows/family.wdl`, **v3.3.1**, commit `477ef39` |
| resource bundle | `hifi-wdl-resources-v3.1.0` (GRCh38 no-alt analysis set; TRGT catalog `adotto_strchive_20250827.hg38`; sawfish exclude / expected-CN beds; MethBat CpG islands; PharmCAT 2.15.4 positions) |
| executor | miniwdl with the `slurm_singularity` container backend (SingularityPro 4.1.2), one driver job per family |
| backend / preemption | `HPC`, `preemptible = false` |

Tools run inside the workflow are pinned by container digest in the checkout's `image_manifest.txt`; the digests this
cohort used (first 12 hex characters) were: pbmm2 `d27f27b94e7f`, DeepVariant `google/deepvariant:1.10.0`, GLnexus
`ce6fecf59ddd`, HiPhase `2c54932e4992` (1.6.0), sawfish `18ba096219fe` (2.2.1), TRGT `be0ed7c173d2` (5.0.0), mosdepth
`63f7a5d1a4a1`, paraphase `e7e1bd125019`, mitorsaw `bc1aaa4633d6`, pbstarphase `bad96d68a60d`, MethBat `0a5363af6a8b`,
pb-cpg-tools `afd5468a423f`, svpack `628e9851e425`, slivar `f71a27f756e2`, wgs_tertiary `410597030e0c`, pbtk
`67cd438ed9f3`, pb_wdl_base `4b889a1f21a6`. Human-readable tool versions are in the workflow's `docs/tools_containers.md`
at that tag. The family workflow performs per sample: pbmm2 alignment, merge of SMRT cells, mosdepth, paraphase,
mitorsaw, DeepVariant (make_examples → call_variants → postprocess) and, per family: GLnexus joint small-variant calling,
HiPhase joint phasing and haplotagging (small variants, SVs, TRs), sawfish discover + joint-call, TRGT genotyping,
methylation (pb-cpg-tools / MethBat), pbstarphase, and the tertiary annotation (slivar, svpack).

## Inputs

One `inputs.json` per family with five keys (`config/inputs.template.json`): `humanwgs_family.family` (family id; per
sample: id, one or more unaligned HiFi BAMs, `affected`, `sex`, and for children `father_id` / `mother_id`),
`ref_map_file`, `tertiary_map_file`, `backend`, `preemptible`. `wdl/make_inputs.py` writes them from a cohort manifest
(`manifest/build_manifest_*.py` build the manifest from the cohort's metadata and the discovered per-movie BAMs). The two
map files are `config/*.template.tsv` with `${RESOURCES_ROOT}` filled in.

## Running a family

```bash
cp config/upstream.env.example config/upstream.env      # site paths, account, partition, qos (gitignored)
wdl/make_inputs.py --manifest samples.tsv --outdir $INPUTS_DIR
wdl/run_family.sh FAM01                                  # submits driver + watchdog for one family
```
`run_family.sh` submits `driver_family.sb` (a 48 h Slurm job that runs `miniwdl run family.wdl -i <inputs> --dir
$RUNROOT/run_<FAMILY>`, writing `RESULT.status` atomically at the end) and `watchdog_family.sb` (detects driver death via
`squeue` and relaunches from the miniwdl call cache up to a restart budget). At most three families run concurrently
(`MAX_DRIVERS`); a restart of an existing run is exempt. The miniwdl configuration is `config/miniwdl.cfg.example`
(`~/.config/miniwdl.cfg`): 11 concurrent tasks, call cache on, every task submitted as its own Slurm job.

## Two operational lessons

1. **Per-task time limit is calibrated to coverage.** The default `--time=06:00:00` per task is ample at ~22–24× and
   fails at 46–48× (GIAB Revio release): the largest DeepVariant `make_examples` shard per sample exceeded it three
   times in a row (exit 143). Override per run rather than editing the global config:
   `export MINIWDL__SLURM__EXTRA_ARGS="--partition <p> --account=<a> --time=14:00:00 --nodes=1 --ntasks-per-node=1 --comment humanwgs"`
   before `run_family.sh`; `--export=ALL` carries it into the driver, and the call cache resumes every finished task.
2. **Write the terminal status first.** `driver_family.sb` writes `RESULT.status` before mail or any other trailing
   step, atomically (`mv`), because a status written last made completed runs look non-terminal and the watchdog
   restarted them.

## Cohort-wide joint call

`cohort/glnexus_full.sb` re-runs GLnexus over all per-sample gVCFs (105 here) with `cohort/glnexus_config.yml` to
produce the cohort small-variant callset ("freeze 1") used by `03_panel` onward and by `13_phase_dnm` for synthetic
trios; sizing rationale is in the script header. The per-family callsets from the WDL are what `13_phase_dnm` reviews.

## External trio example

`config/examples/HG002_trio.inputs.json` is the inputs file used for the GIAB Ashkenazi trio (public sample ids), run
through the identical recipe as a reproducibility anchor for `13_phase_dnm` (DESIGN P16).
