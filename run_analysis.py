#!/usr/bin/env python3
import argparse
import shlex
import subprocess
import sys
import os
from pathlib import Path
import yaml

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
    "align": [],
    "all": [],
    "all_no_download":[],
    "trim": [
        "data/trimmed/{donor}/trim.done",
    ],
    "trim_and_qc": [
        "data/trimmed/{donor}/trim.done",
        "results/qc/fastqc/trimmed/{donor}/fastqc.done",
        "results/qc/multiqc/trimmed/multiqc_report.html",
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
        
    ],

    "unlock": [],

}

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
        description="Docker wrapper for Snakemake workflow. You must choose what to run.",
    )

    p.add_argument("--list-donors", action="store_true")
    p.add_argument("--list-sections", action="store_true")

    sub = p.add_subparsers(dest="section")

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--image", default="scrnaseq-workflow")
        sp.add_argument("--snakefile", default="workflow/Snakefile")
        sp.add_argument("--configfile", default="config/config.yaml")
        sp.add_argument("--cpus", type=int, default=8)
        sp.add_argument("--cores", type=int, default=8)
        sp.add_argument("-j", "--jobs", type=int, default=1)
        sp.add_argument("--set-threads", action="append", default=[])
        sp.add_argument("--dry-run", action="store_true")
        sp.add_argument("--rerun-incomplete", dest="rerun_incomplete", action="store_true", default=True, help="Rerun incomplete jobs (default: on).")
        sp.add_argument("--no-rerun-incomplete", dest="rerun_incomplete", action="store_false",  help="Disable rerun of incomplete jobs.")
        sp.add_argument("--rerun-triggers", default="mtime")
        # IMPORTANT: --extra captures the rest; user should pass it LAST.
        sp.add_argument("--extra", nargs=argparse.REMAINDER)

    sp_dl = sub.add_parser("download_data")
    add_common(sp_dl)

    sp_dlqc = sub.add_parser("download_data_and_qc")
    add_common(sp_dlqc)

    sp_qc = sub.add_parser("qc")
    add_common(sp_qc)

    sp_ref = sub.add_parser("ref")
    add_common(sp_ref)

    sp_align = sub.add_parser("align")
    sp_align.add_argument("--donor", action="append")
    sp_align.add_argument("--trimmed", action="store_true")
    add_common(sp_align)

    sp_trim = sub.add_parser("trim")
    add_common(sp_trim)

    sp_trimqc = sub.add_parser("trim_and_qc")
    add_common(sp_trimqc)

    sp_all = sub.add_parser("all")
    sp_all.add_argument("--trimmed", action="store_true")
    add_common(sp_all)
    
    sp_all_no_download = sub.add_parser("all_no_download")
    sp_all_no_download.add_argument("--trimmed", action="store_true")
    add_common(sp_all_no_download)

    sp_unlock = sub.add_parser("unlock")
    add_common(sp_unlock)



    args = p.parse_args()

    if args.list_sections:
        print("Available sections:")
        for s in SECTION_TARGETS.keys():
            print(f"  - {s}")
        return 0

    if args.list_donors:
        for d in load_donors("config/config.yaml"):
            print(d)
        return 0
        


    if not args.section:
        p.error("You must choose a section to run.")

    repo_root = Path.cwd()
    if not (repo_root / args.snakefile).exists():
        print(f"ERROR: Snakefile not found at {args.snakefile} (run from repo root).", file=sys.stderr)
        return 2

    donors = load_donors(args.configfile)

    # Resolve targets
    targets: list[str] = []
    if args.section in ("download_data", "download_data_and_qc", "qc", "trim", "trim_and_qc"):
        for t in SECTION_TARGETS[args.section]:
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
        donors_for_align = load_donors(args.configfile) if args.donor == ["all"] else args.donor
        mode_dir = "trimmed" if getattr(args, "trimmed", False) else "raw"
        targets = [f"results/alignment/starsolo/{mode_dir}/{d}/starsolo.done" for d in donors_for_align]
    
    elif args.section in ("all", "all_no_download"):
        targets = []



    # Build snakemake command
    smk: list[str] = [
        "snakemake",
        "-s", args.snakefile,
        "--configfile", args.configfile,
        "--cores", str(args.cores),
        "-j", str(args.jobs),
        "--reason",
        "-p",
    ]

    # One consolidated --config
    config_overrides: list[str] = []
    
    if args.rerun_incomplete:
        smk.append("--rerun-incomplete")

    if getattr(args, "trimmed", False):
        config_overrides.append("trim_enabled=true")

    if args.section in ("trim", "trim_and_qc"):
        config_overrides.append("trim_enabled=true")

    if args.section in ("download_data", "download_data_and_qc", "all"):
        config_overrides.append("download_fastqs=true")

    if args.section in ("ref", "all"):
        config_overrides.append("build_star_index=true")
    
    if args.section == "all_no_download":
        config_overrides.append("download_fastqs=false")


    if config_overrides:
        smk.append("--config")
        smk.extend(config_overrides)

    if args.dry_run:
        smk.append("-n")

    if args.rerun_triggers:
        smk.extend(["--rerun-triggers", args.rerun_triggers])

    for st in args.set_threads:
        smk.extend(["--set-threads", st])

    # Extra snakemake args
    if args.extra:
        extra = args.extra
        if extra and extra[0] == "--":
            extra = extra[1:]
        smk.extend(extra)

    # Targets MUST be separated by '--' (prevents --rerun-triggers from eating them)
    if targets:
        smk.append("--")
        smk.extend(targets)

    if args.section == "unlock":
        # Build snakemake unlock command
        smk: list[str] = [
            "snakemake",
            "-s", args.snakefile,
            "--configfile", args.configfile,
            "--unlock",
        ]

        docker = [
            "docker", "run", "--rm", "-i",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "-e", "HOME=/tmp",
            "--init",
            "-e", "XDG_CACHE_HOME=/tmp/.cache",
            "-e", "XDG_CONFIG_HOME=/tmp/.config",
            "--cpus", str(args.cpus),
            "-v", f"{str(repo_root)}:/work",
            "-w", "/work",
            args.image,
        ]

        return run(docker + smk)


if __name__ == "__main__":
    raise SystemExit(main())
