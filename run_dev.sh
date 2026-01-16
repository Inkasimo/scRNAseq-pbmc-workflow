cat > run_dev.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CORES="${1:-4}"

docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/work \
  scrnaseq-workflow \
  snakemake --cores "${CORES}" --rerun-triggers mtime
EOF

chmod +x run_dev.sh
