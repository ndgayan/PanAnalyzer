########################################################################################################################
# PanAnalyzer Server - A tool to visualize Anvio Pangenome Analysis results via a web interface
# Author: Gayan Nagoda
# Date: 2025-11-08
########################################################################################################################

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 --pan-db FILE --genomes-db FILE --conda-env NAME [--port PORT]

Required arguments:
  -p, --pan-db FILE        Path to the Anvio pan genome database (.db)
  -g, --genomes-db FILE    Path to the Anvio genomes storage database (.db)
  -e, --conda-env NAME     Conda environment that has Anvio installed

Optional arguments:
  -P, --port PORT          Port for anvi-display-pan (default: 8080)
  -h, --help               Show this help message and exit
EOF
}

resolve_path() {
    local target="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath "$target"
    else
        python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "$target"
    fi
}

PAN_DB=""
GENOMES_DB=""
CONDA_ENV=""
PORT="8080"

while [ $# -gt 0 ]; do
    case "$1" in
        -p|--pan-db)
            PAN_DB="$2"
            shift 2
            ;;
        -g|--genomes-db)
            GENOMES_DB="$2"
            shift 2
            ;;
        -e|--conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -P|--port)
            PORT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [ -z "$PAN_DB" ] || [ -z "$GENOMES_DB" ] || [ -z "$CONDA_ENV" ]; then
    echo "Error: --pan-db, --genomes-db, and --conda-env are required." >&2
    usage >&2
    exit 1
fi

if [ ! -f "$PAN_DB" ]; then
    echo "Error: pan genome database not found: $PAN_DB" >&2
    exit 1
fi

if [ ! -f "$GENOMES_DB" ]; then
    echo "Error: genomes database not found: $GENOMES_DB" >&2
    exit 1
fi

PAN_DB_PATH=$(resolve_path "$PAN_DB")
GENOMES_DB_PATH=$(resolve_path "$GENOMES_DB")

PAN_DB_DIR=$(dirname -- "$PAN_DB_PATH")
PAN_DB_FILE=$(basename -- "$PAN_DB_PATH")
GENOMES_DB_FILE=$(basename -- "$GENOMES_DB_PATH")

eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda activate "$CONDA_ENV"

echo ""
echo "=========================================="
echo "====== PAN-GENOME ANALYSIS RESULTS ======="
echo "=========================================="
echo ""
echo "Current directory: $(pwd)"
echo "Pan genome database: $PAN_DB_PATH"
echo "Genomes storage: $GENOMES_DB_PATH"
echo ""
echo "Checking for existing processes on port $PORT..."
PORT_PID=$(lsof -ti:"$PORT" 2>/dev/null || true)
if [ -n "$PORT_PID" ]; then
    echo "Killing existing process(es) on port $PORT: $PORT_PID"
    kill -9 $PORT_PID 2>/dev/null || true
    sleep 2
fi

echo ""
echo "Launching Anvio interactive interface..."
echo "Press Ctrl+C to stop the server"
echo ""
anvi-display-pan \
    -p "$PAN_DB_PATH" \
    -g "$GENOMES_DB_PATH" \
    --server-only \
    -P "$PORT"