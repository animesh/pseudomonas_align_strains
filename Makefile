PM := python3
SCRIPT := scripts/make_rects.py

.PHONY: rects rects-force

# Default (safe): write auto-named outputs to avoid overwriting any hand-edited pa_comparison.rects.gp
rects:
	$(PM) $(SCRIPT) --out-gp pa_comparison.rects.auto.gp --out-png pa_comparison.rects.auto.png --export-svg --export-pdf

# Force overwrite the canonical pa_comparison.rects.gp (use with care)
rects-force:
	$(PM) $(SCRIPT) --out-gp pa_comparison.rects.gp --out-png pa_comparison.rects.png --export-svg --export-pdf
