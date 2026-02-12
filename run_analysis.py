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
    "all_no_download": [],
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

    # Internal target sets for upstream (do not expose as CLI sections)
    "upstream_raw": [
        "data/raw/{donor}/fastqs.done",
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
        "data/ref/whitelist.done",
        "data/ref/star_index.done",
        "results/alignment/starsolo/raw/{donor}/starsolo.done",
    ],
    "upstream_trimmed": [
        "data/raw/{donor}/fastqs.done",
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
        "data/trimmed/{donor}/trim.done",
        "results/qc/fastqc/trimmed/{donor}/fastqc.done",
        "results/qc/multiqc/trimmed/multiqc_report.html",
        "data/ref/whitelist.done",
        "data/ref/star_index.done",
        "results/alignment/starsolo/trimmed/{donor}/starsolo.done",
    ],
    "upstream_no_download_raw": [
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
        "data/ref/whitelist.done",
        "data/ref/star_index.done",
        "results/alignment/starsolo/raw/{donor}/starsolo.done",
    ],
    "upstream_no_download_trimmed": [
        "results/qc/fastqc/raw/{donor}/fastqc.done",
        "results/qc/multiqc/raw/multiqc_report.html",
        "data/trimmed/{donor}/trim.done",
        "results/qc/fastqc/trimmed/{donor}/fastqc.done",
        "results/qc/multiqc/trimmed/multiqc_report.html",
        "data/ref/whitelist.done",
        "data/ref/star_index.done",
        "results/alignment/starsolo/trimmed/{donor}/starsolo.done",
    ],

    "build_seurat_object_qc_untrimmed": [
        "results/downstream/seurat/untrimmed/{donor}/seurat_and_qc/seurat_qc.done",
    ],
    "build_seurat_object_qc_trimmed": [
        "data/trimmed/{donor}/trim.done",
        "results/downstream/seurat/trimmed/{donor}/seurat_and_qc/seurat_qc.done",
    ],


    "filter_and_normalize_seurat_untrimmed":[
        "results/downstream/seurat/untrimmed/{donor}/seurat_and_qc/seurat_qc.done",
        "results/downstream/seurat/untrimmed/{donor}/seurat_filt_normalized/seurat_filt_normalize.done",
    ],

    "filter_and_normalize_seurat_trimmed":[
        "results/downstream/seurat/trimmed/{donor}/seurat_and_qc/seurat_qc.done",
        "results/downstream/seurat/trimmed/{donor}/seurat_filt_normalized/seurat_filt_normalize.done",
    ],

    "cluster_annotate_seurat_untrimmed": [
        "results/downstream/seurat/untrimmed/{donor}/seurat_cluster_annot/seurat_cluster_annot.done",
    ],
    "cluster_annotate_seurat_trimmed": [
        "results/downstream/seurat/trimmed/{donor}/seurat_cluster_annot/seurat_cluster_annot.done",
    ],
    "deg_and_tost_untrimmed": [
        "results/downstream/deg_and_tost/untrimmed/deg_and_tost/deg_and_tost.done",
    ],
    "deg_and_tost_trimmed": [
        "results/downstream/deg_and_tost/trimmed/deg_and_tost/deg_and_tost.done",
    ],
    "networks_untrimmed": [
        "results/downstream/networks/untrimmed/networks.done",
    ],
    "networks_trimmed": [
        "results/downstream/networks/trimmed/networks.done",
    ],

    "downstream_untrimmed": [
        "results/downstream/seurat/untrimmed/{donor}/seurat_and_qc/seurat_qc.done",
        "results/downstream/seurat/untrimmed/{donor}/seurat_filt_normalized/seurat_filt_normalize.done",
        "results/downstream/seurat/untrimmed/{donor}/seurat_cluster_annot/seurat_cluster_annot.done",
        "results/downstream/deg_and_tost/untrimmed/deg_and_tost/deg_and_tost.done",
        "results/downstream/networks/untrimmed/networks.done",
    ],
    "downstream_trimmed": [
        "results/downstream/seurat/trimmed/{donor}/seurat_and_qc/seurat_qc.done",
        "results/downstream/seurat/trimmed/{donor}/seurat_filt_normalized/seurat_filt_normalize.done",
        "results/downstream/seurat/trimmed/{donor}/seurat_cluster_annot/seurat_cluster_annot.done",
        "results/downstream/deg_and_tost/trimmed/deg_and_tost/deg_and_tost.done",
        "results/downstream/networks/trimmed/networks.done",
    ],

    "unlock": [],
}

#HELPERS

def _docker_bind_mount_works(repo_root: Path, image: str) -> bool:
    """
    Returns True if Docker can see workflow/Snakefile when bind-mounting repo_root.
    """
    test_cmd = [
        "docker", "run", "--rm",
        "-v", f"{str(repo_root)}:/work",
        image,
        "sh", "-lc", "test -f /work/workflow/Snakefile"
    ]
    try:
        subprocess.run(
            test_cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=os.environ,
        )
        return True
    except Exception:
        return False


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
        #sp.add_argument("--image",
            #default="ghcr.io/inkasimo/scrnaseq-pbmc-workflow@sha256:80354b76e76405636c43e73902236e0399d26978a214227afbafa46fc0555bb8")
        sp.add_argument("--image",
            default="ghcr.io/inkasimo/scrnaseq-pbmc-workflow:v1.0.1")
        #sp.add_argument("--image", default="scrnaseq-workflow")
        sp.add_argument("--snakefile", default="workflow/Snakefile")
        sp.add_argument("--configfile", default="config/config.yaml")
        sp.add_argument("--cpus", type=int, default=8)
        sp.add_argument("--cores", type=int, default=8)
        sp.add_argument("-j", "--jobs", type=int, default=1)
        sp.add_argument("--set-threads", action="append", default=[])
        sp.add_argument("--dry-run", action="store_true")
        sp.add_argument("--rerun-incomplete", dest="rerun_incomplete", action="store_true", default=True,
                        help="Rerun incomplete jobs (default: on).")
        sp.add_argument("--no-rerun-incomplete", dest="rerun_incomplete", action="store_false",
                        help="Disable rerun of incomplete jobs.")
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

    # CLI sections: upstream + upstream_no_download (internal target sets are selected later)
    sp_upstream = sub.add_parser("upstream")
    sp_upstream.add_argument("--trimmed", action="store_true")
    add_common(sp_upstream)

    sp_upstream_no_dl = sub.add_parser("upstream_no_download")
    sp_upstream_no_dl.add_argument("--trimmed", action="store_true")
    add_common(sp_upstream_no_dl)

    sp_seurat = sub.add_parser("build_seurat_object_qc")
    sp_seurat.add_argument("--trimmed", action="store_true")
    sp_seurat.add_argument("--donor", action="append") 
    add_common(sp_seurat)

    sp_filt_norm_seurat = sub.add_parser("filter_and_normalize_seurat")
    sp_filt_norm_seurat.add_argument("--trimmed", action="store_true")
    sp_filt_norm_seurat.add_argument("--donor", action="append") 
    add_common(sp_filt_norm_seurat)

    sp_cluster_annot = sub.add_parser("cluster_annotate_seurat")
    sp_cluster_annot.add_argument("--trimmed", action="store_true")
    sp_cluster_annot.add_argument("--donor", action="append")
    add_common(sp_cluster_annot)

    sp_cluster_annot = sub.add_parser("deg_and_tost")
    sp_cluster_annot.add_argument("--trimmed", action="store_true")
    add_common(sp_cluster_annot)

    sp_networks = sub.add_parser("networks")
    sp_networks.add_argument("--trimmed", action="store_true")
    add_common(sp_networks)


    sp_downstream = sub.add_parser("downstream")
    sp_downstream.add_argument("--trimmed", action="store_true")
    add_common(sp_downstream)


    sp_unlock = sub.add_parser("unlock")
    add_common(sp_unlock)

    args = p.parse_args()

    if args.list_sections:
        print("Available sections:")
        for s in [
            "download_data",
            "download_data_and_qc",
            "qc",
            "ref",
            "trim",
            "trim_and_qc",
            "align",
            "upstream",
            "upstream_no_download",
            "all",
            "all_no_download",
            "build_seurat_object_qc",
            "filter_and_normalize_seurat",
            "cluster_annotate_seurat",
            "deg_and_tost",
            "networks",
            "downstream",
            "unlock",
        ]:
            print(f"  - {s}")
        return 0

    if args.list_donors:
        for d in load_donors("config/config.yaml"):
            print(d)
        return 0

    if not args.section:
        p.error("You must choose a section to run.")

    try:
        repo_root = Path.cwd()
    except FileNotFoundError:
        print(
            "ERROR: Current working directory no longer exists.\n"
            "You likely moved/renamed/deleted the folder you were in.\n"
            "cd into the repo root and re-run.",
            file=sys.stderr,
        )
        return 2

    if not (repo_root / args.snakefile).exists():
        print(f"ERROR: Snakefile not found at {args.snakefile} (run from repo root).", file=sys.stderr)
        return 2

    donors = load_donors(args.configfile)

    # Resolve targets
    targets: list[str] = []

    if args.section in ("download_data", "download_data_and_qc", "qc", "trim", "trim_and_qc"):
        for t in SECTION_TARGETS[args.section]:
            if "{donor}" in t:
                targets.extend(t.format(donor=d) for d in donors)
            else:
                targets.append(t)

    elif args.section == "ref":
        targets = SECTION_TARGETS["ref"]

    elif args.section == "align":
        if not args.donor:
            print("ERROR: align requires --donor (or --donor all).", file=sys.stderr)
            return 2
        donors_for_align = donors if args.donor == ["all"] else args.donor
        mode_dir = "trimmed" if getattr(args, "trimmed", False) else "raw"
        targets = [f"results/alignment/starsolo/{mode_dir}/{d}/starsolo.done" for d in donors_for_align]

    elif args.section in ("upstream", "upstream_no_download"):
        if args.section == "upstream":
            key = "upstream_trimmed" if getattr(args, "trimmed", False) else "upstream_raw"
        else:
            key = "upstream_no_download_trimmed" if getattr(args, "trimmed", False) else "upstream_no_download_raw"

        for t in SECTION_TARGETS[key]:
            if "{donor}" in t:
                targets.extend(t.format(donor=d) for d in donors)
            else:
                targets.append(t)


    elif args.section == "build_seurat_object_qc":
        key = "build_seurat_object_qc_trimmed" if getattr(args, "trimmed", False) else "build_seurat_object_qc_untrimmed"

        # donors selected
        selected = donors
        if getattr(args, "donor", None):
            selected = donors if args.donor == ["all"] else args.donor

        for t in SECTION_TARGETS[key]:
            if "{donor}" in t:
                targets.extend(t.format(donor=d) for d in selected)
            else:
                targets.append(t)
    
    elif args.section == "filter_and_normalize_seurat":
        key = "filter_and_normalize_seurat_trimmed" if getattr(args, "trimmed", False) else "filter_and_normalize_seurat_untrimmed"

        # donors selected
        selected = donors
        if getattr(args, "donor", None):
            selected = donors if args.donor == ["all"] else args.donor

        for t in SECTION_TARGETS[key]:
            if "{donor}" in t:
                targets.extend(t.format(donor=d) for d in selected)
            else:
                targets.append(t)

    elif args.section == "cluster_annotate_seurat":
        key = "cluster_annotate_seurat_trimmed" if getattr(args, "trimmed", False) else "cluster_annotate_seurat_untrimmed"

        selected = donors
        if getattr(args, "donor", None):
            selected = donors if args.donor == ["all"] else args.donor

        for t in SECTION_TARGETS[key]:
            if "{donor}" in t:
                targets.extend(t.format(donor=d) for d in selected)
            else:
                targets.append(t)

    elif args.section == "deg_and_tost":
        key = "deg_and_tost_trimmed" if getattr(args, "trimmed", False) else "deg_and_tost_untrimmed"
        targets = SECTION_TARGETS[key]

    
    elif args.section == "networks":
        key = "networks_trimmed" if getattr(args, "trimmed", False) else "networks_untrimmed"
        targets = SECTION_TARGETS[key]



    
    elif args.section == "downstream":
        key = "downstream_trimmed" if getattr(args, "trimmed", False) else "downstream_untrimmed"
        for t in SECTION_TARGETS[key]:
            if "{donor}" in t:
                targets.extend(t.format(donor=d) for d in donors)
            else:
                targets.append(t)



    elif args.section in ("all", "all_no_download", "unlock"):
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

    if args.section in ("download_data", "download_data_and_qc", "all", "upstream"):
        config_overrides.append("download_fastqs=true")

    if args.section in ("ref", "all", "upstream", "upstream_no_download"):
        config_overrides.append("build_star_index=true")

    if args.section in ("all_no_download", "upstream_no_download"):
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

    # Fail fast if Docker bind mount is empty/unavailable (prevents misleading Snakefile errors)
    if not _docker_bind_mount_works(repo_root, args.image):
        rr = str(repo_root)
        hint = ""
        if rr.startswith("/mnt/"):
            hint = "Hint: running from /mnt/* (Windows drive mount in WSL). Docker Desktop may block these mounts.\n"
        elif rr.startswith("/home/") or rr.startswith("/root/"):
            hint = "Hint: running from WSL Linux filesystem path. Some Docker Desktop setups cannot mount WSL paths.\n"

        print(
            (
                "ERROR: Docker bind mount failed.\n\n"
                f"Host path: {repo_root}\n"
                "Missing inside container: /work/workflow/Snakefile\n\n"
                f"{hint}"
                "Check Docker file-sharing settings or run on Linux/macOS.\n"
            ),
            file=sys.stderr,
        )
        sys.exit(2)



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

    if args.section == "unlock":
        smk = [
            "snakemake",
            "-s", args.snakefile,
            "--configfile", args.configfile,
            "--unlock",
        ]
        return run(docker + smk)

    # NORMAL PATH: run whatever smk you already built above
    return run(docker + smk)

if __name__ == "__main__":
    raise SystemExit(main())
