"""Replay saved 2024 entropy experiments for the October 2026 portfolio.

This analysis was added during portfolio preparation in October 2026. It reads
the archived trial data, checks the supplied aggregate statistics, and renders
them in dimensionless entropy units. It does not rerun particle dynamics or
fit a new model. Recorded timestamps are unchanged: the original simulation
records a post-step state at the initial zero timestamp. See PROVENANCE.md for
the source timing and model limitations.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
SCALE = 1.38e-23 * 1e22
TRIAL_COLUMNS = ["entropy", *[f"entropy_{i}" for i in range(1, 50)]]
AGGREGATE_ATOL = 1e-10
AGGREGATE_RTOL = 1e-10


@dataclass(frozen=True)
class Experiment:
    """Validated trial statistics, expressed as H = stored_entropy / SCALE."""

    filename: str
    time: np.ndarray
    mean: np.ndarray
    std: np.ndarray
    mean_discrepancy: float
    std_discrepancy: float

    def summary(self) -> dict:
        return {
            "source_file": f"data/{self.filename}",
            "trials": len(TRIAL_COLUMNS),
            "points": int(self.time.size),
            "recorded_time_range_seconds": [
                float(self.time[0]),
                float(self.time[-1]),
            ],
            "final_mean_dimensionless": float(self.mean[-1]),
            "final_population_std_dimensionless": float(self.std[-1]),
            "max_absolute_archive_discrepancy_stored_units": {
                "mean": self.mean_discrepancy,
                "population_std": self.std_discrepancy,
            },
        }


def load_experiment(path: Path) -> Experiment:
    """Validate all 50 trials and reproduce archived means and ddof=0 SDs."""
    frame = pd.read_excel(path, engine="openpyxl")
    found_trials = {
        column
        for column in frame.columns
        if isinstance(column, str) and re.fullmatch(r"entropy(?:_\d+)?", column)
    }
    expected_trials = set(TRIAL_COLUMNS)
    if found_trials != expected_trials:
        missing = sorted(expected_trials - found_trials)
        extra = sorted(found_trials - expected_trials)
        raise ValueError(
            f"{path.name}: expected exactly 50 trial columns; "
            f"missing={missing}, unexpected={extra}"
        )

    required = ["time", *TRIAL_COLUMNS, "entropy_ave", "entropy_std"]
    missing_columns = set(required) - set(frame.columns)
    if missing_columns:
        raise ValueError(f"{path.name}: missing columns {sorted(missing_columns)}")
    if len(frame) < 2:
        raise ValueError(f"{path.name}: at least two time points are required")

    try:
        numeric = frame.loc[:, required].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{path.name}: all required values must be numeric") from exc
    if not np.isfinite(numeric).all():
        raise ValueError(f"{path.name}: all required values must be finite")

    time = numeric[:, 0]
    if not (np.diff(time) > 0).all():
        raise ValueError(f"{path.name}: recorded time must be strictly increasing")

    trials = numeric[:, 1:51]
    mean = trials.mean(axis=1)
    std = trials.std(axis=1, ddof=0)
    archived_mean, archived_std = numeric[:, 51], numeric[:, 52]
    differences = {}
    for label, computed, archived in (
        ("mean", mean, archived_mean),
        ("population_std", std, archived_std),
    ):
        difference = float(np.max(np.abs(computed - archived)))
        differences[label] = difference
        if not np.allclose(
            computed, archived, atol=AGGREGATE_ATOL, rtol=AGGREGATE_RTOL
        ):
            raise ValueError(
                f"{path.name}: recomputed {label} disagrees with archive "
                f"(maximum absolute difference {difference:.6g} in stored units)"
            )

    return Experiment(
        filename=path.name,
        time=time,
        mean=mean / SCALE,
        std=std / SCALE,
        mean_discrepancy=differences["mean"],
        std_discrepancy=differences["population_std"],
    )


def render_figure(ideal_gas: Experiment, random_walk: Experiment, path: Path) -> None:
    """Draw independent axes and explicitly label the trial-to-trial spread."""
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "axes.edgecolor": "#C6CFD9",
        "axes.labelcolor": "#233247",
        "text.color": "#233247",
        "xtick.color": "#516074",
        "ytick.color": "#516074",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "grid.color": "#E7EDF3",
        "grid.linewidth": 0.7,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.6))
        panels = (
            (ideal_gas, "Particle transport through an aperture", "#126DAD", np.log(2), "ln(2)"),
            (random_walk, "Random walk in a reflecting cube", "#138678", np.log(1000), "ln(1000)"),
        )
        for ax, (experiment, title, color, bound, bound_label) in zip(axes, panels):
            ax.fill_between(
                experiment.time,
                experiment.mean - experiment.std,
                experiment.mean + experiment.std,
                color=color,
                alpha=0.22,
                linewidth=0,
                label="Mean ± 1 population SD",
                zorder=2,
            )
            ax.plot(
                experiment.time,
                experiment.mean,
                color=color,
                linewidth=1.8,
                label="Mean of 50 trials",
                zorder=3,
            )
            ax.axhline(
                bound,
                color="#7E6A50",
                linestyle=(0, (4, 3)),
                linewidth=1.2,
                label=f"Entropy upper bound: {bound_label}",
                zorder=1,
            )
            ax.set_title(title, loc="left", pad=14, fontweight="bold")
            ax.set_xlabel("Recorded time (s)", labelpad=9)
            ax.set_ylabel("Dimensionless entropy H", labelpad=9)
            ax.set_xlim(float(experiment.time[0]), float(experiment.time[-1]))
            low = min(0.0, float(np.min(experiment.mean - experiment.std)))
            high = max(float(bound), float(np.max(experiment.mean + experiment.std)))
            ax.set_ylim(low - 0.03 * (high - low), high + 0.10 * (high - low))
            ax.grid(axis="y")
            ax.set_axisbelow(True)
            ax.legend(loc="lower right", frameon=False, fontsize=8.5)
            ax.ticklabel_format(axis="x", style="plain", useOffset=False)

        fig.suptitle(
            "Entropy relaxation · archived 2024 experiments",
            x=0.065,
            y=0.965,
            ha="left",
            fontsize=17,
            fontweight="bold",
        )
        fig.text(
            0.065,
            0.045,
            "Recomputed from 50 saved trials per experiment. Shading shows ±1 SD, not a confidence interval.\n"
            "H = stored entropy / 0.138. Original recorded timestamps are preserved; no new simulation or fit.",
            fontsize=9,
            color="#516074",
            linespacing=1.5,
        )
        fig.subplots_adjust(left=0.065, right=0.98, top=0.83, bottom=0.23, wspace=0.24)
        fig.savefig(path, dpi=180, facecolor="white")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate and visualize the archived 2024 entropy trials; no simulation or fitting."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        metavar="DIR",
        help="write DIR/figures and DIR/results (default: this script's directory)",
    )
    args = parser.parse_args()

    try:
        ideal_gas = load_experiment(ROOT / "data" / "ideal_gas.xlsx")
        random_walk = load_experiment(ROOT / "data" / "random_walk.xlsx")
        output_root = args.output_dir.expanduser()
        figures = output_root / "figures"
        results = output_root / "results"
        figures.mkdir(parents=True, exist_ok=True)
        results.mkdir(parents=True, exist_ok=True)

        report = {
            "schema_version": 1,
            "analysis": "replay_of_archived_trials_without_new_simulation_or_fitting",
            "entropy_conversion": {
                "formula": "H = stored_entropy / (1.38e-23 * 1e22)",
                "stored_units_per_dimensionless_entropy": SCALE,
            },
            "standard_deviation_ddof": 0,
            "uncertainty_band": "mean_plus_or_minus_one_population_standard_deviation_not_a_confidence_interval",
            "timestamp_note": "Recorded timestamps are preserved. Original source records a post-step state at time zero.",
            "aggregate_validation_tolerance_stored_units": {
                "absolute": AGGREGATE_ATOL,
                "relative": AGGREGATE_RTOL,
            },
            "ideal_gas": ideal_gas.summary(),
            "random_walk": random_walk.summary(),
        }
        figure_path = figures / "entropy_summary.png"
        summary_path = results / "summary.json"
        render_figure(ideal_gas, random_walk, figure_path)
        summary_path.write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Error: {exc}\n")

    print("Validated both archives: 50 trials each; mean and population SD agree.")
    print(f"Figure: {figure_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
