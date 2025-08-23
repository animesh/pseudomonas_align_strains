#!/usr/bin/env bash
set -euo pipefail

# Reproducible pipeline wrapper for the repo. Idempotent: skips steps if outputs exist.
# Writes a run.log in the repo root.

LOG=run.log
exec > >(tee -a "$LOG") 2>&1

timestamp(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
log(){ echo "[$(timestamp)] $*"; }

# Requirements check
REQUIRED=(nucmer delta-filter show-coords mummerplot gnuplot python3)
MISSING=()
for cmd in "${REQUIRED[@]}"; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    MISSING+=("$cmd")
  fi
done
if [ ${#MISSING[@]} -ne 0 ]; then
  log "Missing required commands: ${MISSING[*]}" >&2
  log "Please install the missing tools and re-run. Aborting."
  exit 2
fi

# Optional tools
MINIMAP2_OK=0
if command -v minimap2 >/dev/null 2>&1; then
  MINIMAP2_OK=1
  log "minimap2 found — optional sensitive mapping available"
else
  log "minimap2 not found — skipping optional minimap2 steps"
fi

# 0) Prepare FASTA files
ATCC=ATCC.fasta
ST=ST.fasta
ATCC_FIXED=ATCC.fixed.fasta
ST_FIXED=ST.fixed.fasta
BLA=blaVIM2.fasta

log "Step 0: FASTA preparation"
if [ ! -s "$ATCC_FIXED" ]; then
  if [ -s "$ATCC" ]; then
    if [ -x scripts/clean_fasta.sh ]; then
      log "Normalizing $ATCC -> $ATCC_FIXED"
      scripts/clean_fasta.sh "$ATCC" > "$ATCC_FIXED"
    else
      log "scripts/clean_fasta.sh missing or not executable; copying $ATCC -> $ATCC_FIXED"
      cp "$ATCC" "$ATCC_FIXED"
    fi
  else
    log "$ATCC not found; aborting" >&2
    exit 3
  fi
else
  log "$ATCC_FIXED exists — skipping FASTA cleaning for ATCC"
fi

if [ ! -s "$ST_FIXED" ]; then
  if [ -s "$ST" ]; then
    if [ -x scripts/clean_fasta.sh ]; then
      log "Normalizing $ST -> $ST_FIXED"
      scripts/clean_fasta.sh "$ST" > "$ST_FIXED"
    else
      log "scripts/clean_fasta.sh missing or not executable; copying $ST -> $ST_FIXED"
      cp "$ST" "$ST_FIXED"
    fi
  else
    log "$ST not found; aborting" >&2
    exit 3
  fi
else
  log "$ST_FIXED exists — skipping FASTA cleaning for ST"
fi

# 1) Pairwise alignment
log "Step 1: pairwise alignment (nucmer)"
if [ ! -s pa_comparison.delta ]; then
  log "Running nucmer --maxmatch -p pa_comparison $ATCC_FIXED $ST_FIXED"
  nucmer --maxmatch -p pa_comparison "$ATCC_FIXED" "$ST_FIXED"
else
  log "pa_comparison.delta already exists — skipping nucmer"
fi

# 2) Delta filtering and coords at different sensitivities
log "Step 2: delta-filter and show-coords (strict + relaxed)"
if [ -s pa_comparison.delta ]; then
  if [ ! -s pa_comparison.filtered.delta ]; then
    log "Creating strict filtered delta (i=90 l=100)"
    delta-filter -i 90 -l 100 -q pa_comparison.delta > pa_comparison.filtered.delta || true
  else
    log "pa_comparison.filtered.delta exists — skipping strict filter"
  fi
  if [ ! -s pa_comparison.coords ]; then
    if [ -s pa_comparison.filtered.delta ]; then
      show-coords -rclT pa_comparison.filtered.delta > pa_comparison.coords || true
      log "Wrote pa_comparison.coords"
    else
      log "Strict filtered delta missing; creating coords from pa_comparison.delta"
      show-coords -rclT pa_comparison.delta > pa_comparison.coords || true
      log "Wrote pa_comparison.coords (from raw delta)"
    fi
  else
    log "pa_comparison.coords exists — skipping show-coords (strict)"
  fi
  # relaxed
  if [ ! -s pa_comparison.relaxed.filtered.delta ]; then
    log "Creating relaxed filtered delta (i=90 l=20)"
    delta-filter -i 90 -l 20 -q pa_comparison.delta > pa_comparison.relaxed.filtered.delta || true
  else
    log "pa_comparison.relaxed.filtered.delta exists — skipping"
  fi
  if [ ! -s pa_comparison.relaxed.coords ]; then
    if [ -s pa_comparison.relaxed.filtered.delta ]; then
      show-coords -rclT pa_comparison.relaxed.filtered.delta > pa_comparison.relaxed.coords || true
      log "Wrote pa_comparison.relaxed.coords"
    fi
  else
    log "pa_comparison.relaxed.coords exists — skipping"
  fi
else
  log "pa_comparison.delta missing; cannot run delta-filter/show-coords" >&2
fi

# 3) Produce dotplots (mummerplot + sanitize gp + gnuplot)
log "Step 3: mummerplot and render PNGs"
GP_PREFIX=pa_comparison
if [ -s pa_comparison.relaxed.filtered.delta ]; then
  if [ ! -s ${GP_PREFIX}.relaxed.png ]; then
    log "Generating mummerplot files (relaxed)"
    mummerplot --fat --layout -p "$GP_PREFIX" pa_comparison.relaxed.filtered.delta || true
    if [ -s ${GP_PREFIX}.relaxed.gp ]; then
      sed '/set mouse/d' ${GP_PREFIX}.relaxed.gp > ${GP_PREFIX}.relaxed.sanitized.gp || true
      # ensure png terminal
      if ! grep -q "set terminal.*pngcairo" ${GP_PREFIX}.relaxed.sanitized.gp 2>/dev/null; then
        printf "%s\n" "set terminal pngcairo enhanced font 'Arial,10'" "set output '${GP_PREFIX}.relaxed.png'" | cat - ${GP_PREFIX}.relaxed.sanitized.gp > ${GP_PREFIX}.relaxed.sanitized.tmp && mv ${GP_PREFIX}.relaxed.sanitized.tmp ${GP_PREFIX}.relaxed.sanitized.gp
      fi
      gnuplot ${GP_PREFIX}.relaxed.sanitized.gp || true
      log "Rendered ${GP_PREFIX}.relaxed.png"
    else
      log "No sanitized gp generated for relaxed plot"
    fi
  else
    log "${GP_PREFIX}.relaxed.png exists — skipping plot render"
  fi
else
  # fallback to raw delta
  if [ -s pa_comparison.delta ] && [ ! -s ${GP_PREFIX}.png ]; then
    log "Generating mummerplot (raw delta)"
    mummerplot --fat --layout -p "$GP_PREFIX" pa_comparison.delta || true
    if [ -s ${GP_PREFIX}.gp ]; then
      sed '/set mouse/d' ${GP_PREFIX}.gp > ${GP_PREFIX}.sanitized.gp || true
      if ! grep -q "set terminal.*pngcairo" ${GP_PREFIX}.sanitized.gp 2>/dev/null; then
        printf "%s\n" "set terminal pngcairo enhanced font 'Arial,10'" "set output '${GP_PREFIX}.png'" | cat - ${GP_PREFIX}.sanitized.gp > ${GP_PREFIX}.sanitized.tmp && mv ${GP_PREFIX}.sanitized.tmp ${GP_PREFIX}.sanitized.gp
      fi
      gnuplot ${GP_PREFIX}.sanitized.gp || true
      log "Rendered ${GP_PREFIX}.png"
    fi
  else
    log "No delta for plotting found; skipping mummerplot" 
  fi
fi

# 4) Align blaVIM2 to ST and ATCC
log "Step 4: align blaVIM2 to ST and ATCC"
if [ -s "$BLA" ]; then
  if [ ! -s bla_ST.coords ]; then
    nucmer --prefix=bla_ST "$ST_FIXED" "$BLA" || true
    delta-filter -i 90 -l 20 -q bla_ST.delta > bla_ST.filtered.delta || true
    show-coords -rclT bla_ST.filtered.delta > bla_ST.coords || true
    log "Wrote bla_ST.coords"
  else
    log "bla_ST.coords exists — skipping"
  fi
  if [ ! -s bla_ATCC.coords ]; then
    nucmer --prefix=bla_ATCC "$ATCC_FIXED" "$BLA" || true
    delta-filter -i 90 -l 20 -q bla_ATCC.delta > bla_ATCC.filtered.delta || true
    show-coords -rclT bla_ATCC.filtered.delta > bla_ATCC.coords || true
    log "Wrote bla_ATCC.coords"
  else
    log "bla_ATCC.coords exists — skipping"
  fi
else
  log "$BLA missing — skipping bla mapping"
fi

# 5) Create rectangle overlays
log "Step 5: create rectangle overlays with scripts/make_rects.py"
if [ -x scripts/make_rects.py ] || [ -s scripts/make_rects.py ]; then
  # prefer relaxed.coords if present
  PA_COORDS=pa_comparison.relaxed.coords
  if [ ! -s "$PA_COORDS" ]; then
    PA_COORDS=pa_comparison.coords
  fi
  if [ -s "$PA_COORDS" ]; then
    python3 scripts/make_rects.py "$PA_COORDS" || true
    log "Ran make_rects.py on $PA_COORDS"
  else
    log "No pa_comparison coords file found for make_rects.py — skipping"
  fi
else
  log "scripts/make_rects.py missing — cannot create overlays"
fi

# 6) Optional: minimap2 per-contig mapping (if available)
if [ $MINIMAP2_OK -eq 1 ]; then
  log "Step 6: Optional minimap2 mapping of subset ST contigs (if subset exists)"
  if [ -s subset_st_contigs.fasta ]; then
    minimap2 -x asm5 -t2 "$ATCC_FIXED" subset_st_contigs.fasta > minimap_st_vs_atcc.paf || true
    log "minimap2 produced minimap_st_vs_atcc.paf (may be empty)"
  else
    log "subset_st_contigs.fasta not present — skipping minimap2 per-contig mapping"
  fi
fi

# 7) Collect artifacts
log "Step 7: collect artifacts into artifacts/"
ARTDIR=artifacts
mkdir -p "$ARTDIR"
for f in pa_comparison*.coords pa_comparison*.delta pa_comparison*.png pa_comparison*.gp pa_comparison.unqry bla_*.coords pa_comparison.rects.*; do
  if [ -e "$f" ]; then
    cp -a "$f" "$ARTDIR/" || true
  fi
done
cp -a scripts "$ARTDIR/" || true
if [ -s run.log ]; then
  cp -a run.log "$ARTDIR/"
fi

tar czf artifacts.tar.gz "$ARTDIR" || true
log "Packed artifacts to artifacts.tar.gz"

# Summary
log "Run complete. Summary of key files written (if present):"
ls -ltrh pa_comparison*.coords || true
ls -ltrh pa_comparison*.png || true
ls -ltrh bla_*.coords || true
log "run.log is available in the repo root (this script logs to run.log)"

exit 0
