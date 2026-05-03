# AGENTS.md

Numerical Linear Algebra (数值代数) course implementations in C++17.

## Project Overview

Author: Quinn — a collection of numerical algorithms for homework assignments.

## Tech Stack

- **Language:** C++17 (`g++`)
- **Build:** GNU Make
- **Plotting:** Python (matplotlib)

## Build Commands

```bash
make                  # Build all targets
make run-hwN          # Build & run homework N (e.g. make run-hw2)
make run-inverseNorm  # Build & run inverse norm estimation test
make run-qr           # Build & run QR decomposition test
make run-timeConsumingTest START_N=10 END_N=500 STEP=10  # Parameterized benchmark
make clean            # Remove build artifacts
```

Build artifacts go to `build/`. Output files go to `output/`.

## Directory Layout

```
include/   — Headers (Matrix class, algorithm declarations)
src/       — Algorithm implementations (.cpp)
test/      — Homework programs & tests (homework02-08, test_*.cpp)
notes/     — Experiment notes (markdown)
output/    — Program output (gitignored)
build/     — Compiled binaries & objects (gitignored)
```

## Coding Conventions

- **Comments in English.**
- Mathematical formulas use LaTeX in comments/docstrings (e.g., `// Compute $\|A^{-1}\|_\infty$`).
- Follow existing naming style: `PascalCase` for classes/structs, `camelCase` for functions, `UPPER_CASE` for constants.
- Header guards use `#pragma once`.
- Keep algorithm implementations in `src/`, declarations in `include/`, homework drivers in `test/`.

## Git Workflow

- **main** — Stable, debugged code only.
- **feat/** — New features or homework branches. Merge back to main after validation.
- Always check `git status` before starting new work to avoid losing uncommitted changes.
- Remote is GitHub HTTPS (unreliable from mainland China); consider SSH key if persistent issues arise.

## Architecture Notes

- `Matrix` class (matrix.hpp) is the core data structure — all algorithms operate on it.
- `operations.hpp` provides basic matrix arithmetic (add, subtract, multiply, transpose, norms).
- Algorithms follow column-major or row-major ordering as needed; check existing code before adding new routines.
- Homework programs in `test/` are self-contained drivers that `#include` headers and link against compiled objects.

## Notes

- `plot_time_consuming.py` generates time-complexity line charts from benchmark output.
- `notes/` contains investigation writeups (e.g., Cholesky performance anomaly analysis).
