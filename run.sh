#!/usr/bin/env bash
set -euo pipefail

IMAGE="${IMAGE:-scrnaseq-workflow}"
WORKDIR="/work"

usage() {
  cat <<EOF
Usage:
  ./run.sh build
  ./run.sh dry [snakemake args...]
  ./run.sh run [snakemake args...]
  ./run.sh target <file> [snakemake args...]
  ./run.sh shell

Examples:
  ./run.sh build
  ./run.sh dry -n -p
  ./run.sh run -p -j 1 data/ref/whitelist.done
  ./run.sh target data/ref/star_index.done -p -j 1
EOF
}

docker_run() {
  docker run --rm -it \
    -v "$(pwd)":${WORKDIR} \
    -w ${WORKDIR} \
    "${IMAGE}" \
    "$@"
}

cmd="${1:-}"
shift || true

case "${cmd}" in
  build)
    docker build -t "${IMAGE}" .
    ;;
  dry)
    docker_run snakemake -n -p "$@"
    ;;
  run)
    docker_run snakemake "$@"
    ;;
  target)
    target="${1:-}"
    shift || true
    if [[ -z "${target}" ]]; then
      echo "ERROR: target file required"
      usage
      exit 1
    fi
    docker_run snakemake "${target}" "$@"
    ;;
  shell)
    docker_run bash
    ;;
  ""|-h|--help|help)
    usage
    ;;
  *)
    echo "ERROR: unknown command: ${cmd}"
    usage
    exit 1
    ;;
esac
