#!/usr/bin/env python3
import argparse
import shlex
import subprocess
import sys
from pathlib import Path

SECTION_TARGETS = {
    "download_data": [
        "data/raw/{donor}/fastqs.done",
    ],
    "download_data_and_qc": [
        "data/raw/{donor}/fastqs.done",
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
    ],
    "qc": [
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
    ],
    "ref": [
        "data/ref/whitelist.done",
        "data/ref/star_index.done",
    ],
    "all": [],
}


import yaml

def load_donors(configfile: str) -> list[str]:
    with open(configfile) as f:
        cfg = yaml.safe_load(f)
    return list(cfg["pbmc"].keys())

def q(x: str) -> str:
    return shlex.quote(x)

def run(cmd: list[str]) -> int:
    print("Running:\n  " + " \\\n  ".join(q(c) for c in cmd))
    return subprocess.call(cmd)

def main() -> int:
    p = argparse.ArgumentParser(
        prog="run_analysis.py",
        description="Docker wrapper for Snakemake workflow. You must choose what to run."
    )
    sub = p.add_subparsers(dest="section", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--image", default="scrnaseq-workflow", help="Docker image name/tag.")
        sp.add_argument("--snakefile", default="workflow/Snakefile", help="Path inside repo.")
        sp.add_argument("--configfile", default="config/config.yaml", help="Path inside repo.")
        sp.add_argument("--cpus", type=int, default=8, help="Docker --cpus value.")
        sp.add_argument("--cores", type=int, default=8, help="Snakemake --cores value.")
        sp.add_argument("-j", "--jobs", type=int, default=1, help="Snakemake -j (parallel jobs).")
        sp.add_argument("--set-threads", action="append", default=[],
                        help='Repeatable, e.g. --set-threads starsolo=8')
        sp.add_argument("--dry-run", action="store_true", help="Snakemake -n.")
        sp.add_argument("--rerun-incomplete", action="store_true", default=True)
        sp.add_argument("--rerun-triggers", default="mtime", help="Snakemake --rerun-triggers value.")
        sp.add_argument("--extra", nargs=argparse.REMAINDER,
                        help="Extra args passed to snakemake after '--'. Example: -- --keep-going")

    sp_dl = sub.add_parser("download_data", help="Download FASTQs for all donors (forces io.download_fastqs=true).")
    add_common(sp_dl)

    sp_dlqc = sub.add_parser("download_data_and_qc", help="Download FASTQs + FastQC + MultiQC (forces io.download_fastqs=true).")
    add_common(sp_dlqc)

    sp_qc = sub.add_parser("qc", help="Run FASTQ acquisition/presence + FastQC + MultiQC.")
    add_common(sp_qc)

    sp_ref = sub.add_parser("ref", help="Ensure whitelist + STAR index are present/built.")
    add_common(sp_ref)

    sp_align = sub.add_parser("align", help="Run STARsolo alignment.")
    sp_align.add_argument("--donor", action="append", help="Repeatable. Example: --donor donor1 --donor donor2")
    add_common(sp_align)

    sp_all = sub.add_parser("all", help="Run full workflow (explicit).")
    add_common(sp_all)
   

    args = p.parse_args()

    repo_root = Path.cwd()
    if not (repo_root / args.snakefile).exists():
        print(f"ERROR: Snakefile not found at {args.snakefile} (run from repo root).", file=sys.stderr)
        return 2

    targets: list[str] = []
    donors = load_donors(args.configfile)

    if args.section in ("download_data", "download_data_and_qc", "qc"):
        raw_targets = SECTION_TARGETS[args.section]
        for t in raw_targets:
            if "{donor}" in t:
                targets.extend([t.format(donor=d) for d in donors])
            else:
                targets.append(t)

    elif args.section == "ref":
        targets = SECTION_TARGETS["ref"]

    elif args.section == "align":
        if not args.donor:
            print("ERROR: align requires --donor (or --donor all).", file=sys.stderr)
            return 2

        if args.donor == ["all"]:
            donors = load_donors(args.configfile)
        else:
            donors = args.donor

        targets = [f"results/alignment/starsolo/{d}/starsolo.done" for d in donors]

    elif args.section == "all":
        targets = []


    # Build snakemake command (inside container)
    smk = [
        "snakemake",
        "-s", args.snakefile,
        "--configfile", args.configfile,
        "--cores", str(args.cores),
        "-j", str(args.jobs),
        "--reason",
        "-p",
    ]
    if args.dry_run:
        smk.append("-n")
    if args.rerun_incomplete:
        smk.append("--rerun-incomplete")
    if args.rerun_triggers:
        smk.extend(["--rerun-triggers", args.rerun_triggers])
    for st in args.set_threads:
        smk.extend(["--set-threads", st])
    if args.section in ("download_data", "download_data_and_qc"):
        smk.extend(["--config", "io.download_fastqs=true"])
    # Force STAR index build when running ref/all (or a dedicated command)
    if args.section == "ref":
        smk.extend([
            "--config", "ref.build_star_index=true",
        ])

    if args.section == "all":
        smk.extend([
            "--config", "ref.build_star_index=true",
            "--config", "io.download_fastqs=true",
        ])



    # Targets go after `--` to keep parsing clean
    smk.append("--")
    smk.extend(targets)

    # Pass-through extra args after `--`
    if args.extra:
        # args.extra already includes the leading '--' if user supplied it; strip it if present
        extra = args.extra
        if extra and extra[0] == "--":
            extra = extra[1:]
        smk.extend(extra)

    # Wrap in docker run
    docker = [
        "docker", "run", "--rm", "-it",
        "--cpus", str(args.cpus),
        "-v", f"{str(repo_root)}:/work",
        "-w", "/work",
        args.image,
    ]
    cmd = docker + smk
    return run(cmd)

if __name__ == "__main__":
    raise SystemExit(main())
