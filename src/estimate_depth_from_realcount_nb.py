#!/usr/bin/env python3
"""Run the project's NB depth estimator from the TKCN source directory."""

from pathlib import Path
import runpy
import sys


implementation = Path(__file__).resolve().parents[2] / "estimate_depth_from_realcount_nb.py"
if not implementation.is_file():
    raise SystemExit(f"ERROR: project implementation not found: {implementation}")

sys.argv[0] = str(implementation)
runpy.run_path(str(implementation), run_name="__main__")
