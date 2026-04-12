#!/usr/bin/env bash
# run.sh — headless MATLAB entrypoint for Haber_PP_MEA
#
# Usage:
#   ./scripts/run.sh smoke                 # sanity-check utilities + toolboxes
#   ./scripts/run.sh preprocess            # build v2.0 caches for all datasets
#   ./scripts/run.sh preprocess doi        # DOI only
#   ./scripts/run.sh preprocess ket        # ketanserin only
#   ./scripts/run.sh preprocess all force  # rebuild even if cache exists
#   ./scripts/run.sh figures               # render every figure (DOI + ket)
#   ./scripts/run.sh figures doi           # DOI figures only
#   ./scripts/run.sh figures ket           # ketanserin figures only
#   ./scripts/run.sh figures doi no-conn   # skip connectivity figures
#   ./scripts/run.sh all                   # smoke + preprocess + figures
#   ./scripts/run.sh matlab "<expr>"       # run an arbitrary MATLAB expression
#
# Logs go to logs/<command>_<timestamp>.log and are tailed to stdout.

set -euo pipefail

# Locate repo root (this script lives in <repo>/scripts/).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Default MATLAB binary (override with MATLAB=/path/to/matlab).
: "${MATLAB:=/Applications/MATLAB_R2023b.app/bin/matlab}"
MATLAB_FLAGS="-nodisplay -nosplash -nodesktop"

# Common MATLAB path setup: add src (recursively) and scripts.
ADDPATH="addpath(genpath(fullfile('${REPO_ROOT}','src')));addpath(fullfile('${REPO_ROOT}','scripts'));"

# Log directory.
mkdir -p "$REPO_ROOT/logs"
stamp="$(date +%Y%m%d_%H%M%S)"

run_matlab() {
    local label="$1"
    local expr="$2"
    local log="$REPO_ROOT/logs/${label}_${stamp}.log"
    local diaryLog="$REPO_ROOT/logs/${label}_${stamp}_matlab_diary.log"
    echo "[run.sh] $label"
    echo "[run.sh] MATLAB   : $MATLAB"
    echo "[run.sh] tee log  : $log"
    echo "[run.sh] diary    : $diaryLog"
    echo ""
    # Wrap the user expression in diary + try/catch so warnings and any
    # uncaught errors are fully captured in the MATLAB-side diary.
    local wrapped
    wrapped="${ADDPATH}diary('${diaryLog}');diary on;try;${expr};catch ME;fprintf(2,'[run.sh] MATLAB error: %s\\n',ME.message);fprintf(2,'%s\\n',getReport(ME,'extended'));diary off;rethrow(ME);end;diary off;"
    (
        cd "$REPO_ROOT"
        "$MATLAB" $MATLAB_FLAGS -batch "$wrapped" 2>&1 | tee "$log"
    )
}

cmd="${1:-help}"
shift || true

case "$cmd" in
    smoke)
        run_matlab smoke "smoke_test"
        ;;

    preprocess)
        study="${1:-all}"
        force="${2:-}"
        if [[ "$force" == "force" ]]; then
            overwrite="true"
        else
            overwrite="false"
        fi
        run_matlab "preprocess_${study}" \
            "run_preprocess('study','${study}','overwrite',${overwrite})"
        ;;

    figures)
        study="${1:-all}"
        conn="${2:-}"
        if [[ "$conn" == "no-conn" ]]; then
            connFlag="false"
        else
            connFlag="true"
        fi
        run_matlab "figures_${study}" \
            "run_figures('study','${study}','connectivity',${connFlag})"
        ;;

    all)
        run_matlab smoke            "smoke_test"
        run_matlab preprocess_all   "run_preprocess('study','all','overwrite',false)"
        run_matlab figures_all      "run_figures('study','all','connectivity',true)"
        ;;

    matlab)
        expr="${1:-}"
        if [[ -z "$expr" ]]; then
            echo "[run.sh] error: ./scripts/run.sh matlab \"<expression>\""
            exit 2
        fi
        run_matlab adhoc "$expr"
        ;;

    help|-h|--help|"")
        cat <<USAGE
run.sh — headless MATLAB entrypoint for Haber_PP_MEA

Commands:
  smoke                          Run scripts/smoke_test.m
  preprocess [study] [force]     Build v2.0 caches (study: all|doi|ket)
  figures    [study] [no-conn]   Render figures (study: all|doi|ket)
  all                            smoke + preprocess + figures
  matlab     "<expr>"            Run an arbitrary MATLAB expression
  help                           Show this message

Environment variables:
  MATLAB       Path to matlab binary (default: $MATLAB)

Examples:
  ./scripts/run.sh smoke
  ./scripts/run.sh preprocess doi
  ./scripts/run.sh preprocess all force
  ./scripts/run.sh figures doi no-conn
  ./scripts/run.sh all
  ./scripts/run.sh matlab "disp(project_config().paths.cache)"

Logs are written to logs/<command>_<timestamp>.log.
USAGE
        ;;

    *)
        echo "[run.sh] unknown command: $cmd"
        echo "[run.sh] try: ./scripts/run.sh help"
        exit 2
        ;;
esac
