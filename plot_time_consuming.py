import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_csv(path: Path):
    n_values = []
    gauss = []
    pgauss = []
    cholesky = []
    modified_cholesky = []

    with path.open(newline="") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        for row in reader:
            n_values.append(float(row["n"].strip()))
            gauss.append(float(row["gaussSolve"].strip()))
            pgauss.append(float(row["PgaussSolve"].strip()))
            cholesky.append(float(row["choleskySolve"].strip()))
            modified_cholesky.append(float(row["modifiedCholeskySolve"].strip()))

    return n_values, gauss, pgauss, cholesky, modified_cholesky


def main() -> None:
    csv_path = Path("./output/time_consuming.csv")
    if not csv_path.exists():
        raise SystemExit(f"CSV not found: {csv_path}")

    n_values, gauss, pgauss, cholesky, modified_cholesky = read_csv(csv_path)

    plt.figure(figsize=(8, 5))
    plt.plot(n_values, gauss, label="gaussSolve", color="#1f77b4")
    plt.plot(n_values, pgauss, label="PgaussSolve", color="#ff7f0e")
    plt.plot(n_values, cholesky, label="choleskySolve", color="#2ca02c")
    plt.plot(n_values, modified_cholesky, label="modifiedCholeskySolve", color="#d62728")

    plt.xlabel("n")
    plt.ylabel("time (s)")
    plt.title("Time Consumption by Algorithm")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = Path("./output/time_consuming.png")
    plt.savefig(output_path, dpi=150)
    print(f"Saved plot to {output_path}")


if __name__ == "__main__":
    main()
