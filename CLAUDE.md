# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Build the package:**
```bash
R CMD build .
R CMD INSTALL .
```

**Check the package (includes basic tests and documentation validation):**
```bash
R CMD check .
```

**Run an interactive example in R:**
```r
library(enviPat)
data(isotopes)
data(chemforms)
results <- isopattern(isotopes, chemforms, threshold = 0.1, charge = FALSE)
```

**Recompile C code only (during development):**
```bash
R CMD INSTALL --preclean .
```

## Architecture

enviPat is an R package with a performance-critical C backend. The general flow is:

1. **R layer (`R/`)** handles user input validation, formula parsing, and result formatting.
2. **C layer (`src/`)** executes the computationally heavy isotope pattern algorithms, called from R via `.Call()` and `useDynLib(enviPat)`.

### Core workflow

- `check_chemform()` → parses and validates chemical formula strings into element counts
- `isopattern()` → takes validated formulas + isotope data, calls C code (`iso_pattern` or `iso_pattern_4`), returns stick patterns (m/z + intensity)
- `envelope()` → convolves stick patterns into Gaussian or Cauchy-Lorentz peak profiles at a given instrument resolution
- `isowrap()` → convenience wrapper combining `getR()`, `isopattern()`, `envelope()`, and optionally `vdetect()`

### C source layout (`src/`)

- `main.c` (~2700 lines) — core isotope pattern generation algorithms
- `combination.c` (~3900 lines) — isotope combination calculations (pruning tree)
- `parse.c/h` — chemical formula string parsing
- `profile.c/h` — peak envelope profile generation
- `peak.c/h` — peak processing utilities
- `element.c/h` — element data handling
- `isotope.c/h` — isotope operations
- `data.h` — isotope masses and abundances embedded as a header
- `preferences.h` — compile-time configuration macros
- `enviPat_init.c` — R registration of C routines

### Data

Four `.rda` files in `data/` are loaded via `data()`:
- `isotopes` — isotope masses and natural abundances (primary input to `isopattern`)
- `chemforms` — example chemical formulas for testing
- `adducts` — adduct definitions for mass spectrometry
- `resolution_list` — reference resolution data from various instruments

### No automated test suite

The repository has no `tests/` directory. Verification is done via `R CMD check` (which runs examples from `.Rd` files) and interactive R sessions.
