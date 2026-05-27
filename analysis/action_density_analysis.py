#!/usr/bin/env python3
"""Action-density loading, summary, and plotting helpers."""

import os
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import matplotlib.pyplot as plt

from calculations import RunData, lattice_spacing_su2, lattice_spacing_su3


ACTION_DENSITY_FILENAMES = (
    "action_density_l.dat",
    "action_densityl.dat",
    "action_density.dat",
)


@dataclass
class ActionDensitySeries:
    run: RunData
    source_file: str
    sweeps: np.ndarray
    values: np.ndarray
    per_timeslice: bool
    note: str


def _candidate_action_file(run_dir: str) -> Optional[str]:
    out_dir = os.path.join(run_dir, "output")
    for name in ACTION_DENSITY_FILENAMES:
        candidate = os.path.join(out_dir, name)
        if os.path.isfile(candidate):
            return candidate
    return None


def _as_sorted_arrays(values_by_sweep: dict[int, list[float]]) -> tuple[np.ndarray, np.ndarray]:
    sweeps = np.array(sorted(values_by_sweep), dtype=int)
    values = np.array([np.mean(values_by_sweep[s]) for s in sweeps], dtype=float)
    return sweeps, values


def load_action_density(run: RunData) -> Optional[ActionDensitySeries]:
    """Load action density for one run.

    Supported formats:
      - two columns: sweep action_density (current output; already global)
      - three columns: sweep t action_density (future/local output; averaged
        over all t)
    """
    filepath = _candidate_action_file(run.run_dir)
    if filepath is None:
        return None

    global_values: dict[int, list[float]] = {}
    local_values: list[tuple[int, int, float]] = []

    with open(filepath) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                sweep = int(parts[0])
            except ValueError:
                continue

            if len(parts) >= 3:
                try:
                    t = int(parts[1])
                    value = float(parts[2])
                    local_values.append((sweep, t, value))
                    continue
                except ValueError:
                    pass

            try:
                value = float(parts[-1])
            except ValueError:
                continue
            global_values.setdefault(sweep, []).append(value)

    if local_values:
        values_by_sweep: dict[int, list[float]] = {}

        for sweep, t, value in local_values:
            values_by_sweep.setdefault(sweep, []).append(value)

        if not values_by_sweep:
            return None
        sweeps, values = _as_sorted_arrays(values_by_sweep)
        note = "per-timeslice input; averaged over all t"
        return ActionDensitySeries(
            run=run,
            source_file=filepath,
            sweeps=sweeps,
            values=values,
            per_timeslice=True,
            note=note,
        )

    if not global_values:
        return None

    sweeps, values = _as_sorted_arrays(global_values)
    note = "global two-column input"
    return ActionDensitySeries(
        run=run,
        source_file=filepath,
        sweeps=sweeps,
        values=values,
        per_timeslice=False,
        note=note,
    )


def collect_action_density(runs: List[RunData]) -> List[ActionDensitySeries]:
    series = []
    for run in runs:
        loaded = load_action_density(run)
        if loaded is not None:
            series.append(loaded)
    return series


def write_action_density_summary(series: List[ActionDensitySeries],
                                 output_dir: str,
                                 gauge_group: str) -> None:
    if not series:
        return

    path = os.path.join(output_dir, f"action_density_summary_{gauge_group}.txt")
    with open(path, "w") as f:
        f.write("# Action-density summary\n")
        f.write(
            "# beta boundary T L seed source n mean_all std_all "
            "start_conf n_post mean_post std_post note\n"
        )
        for item in sorted(series, key=lambda x: (x.run.beta, x.run.boundary, x.run.seed)):
            run = item.run
            post_mask = item.sweeps >= run.start_conf
            post_values = item.values[post_mask]
            mean_all = float(np.mean(item.values))
            std_all = float(np.std(item.values, ddof=1)) if len(item.values) > 1 else 0.0
            mean_post = float(np.mean(post_values)) if len(post_values) else float("nan")
            std_post = float(np.std(post_values, ddof=1)) if len(post_values) > 1 else 0.0
            f.write(
                f"{run.beta:.8g} {run.boundary} {run.T:d} {run.L:d} {run.seed:d} "
                f"{os.path.basename(item.source_file)} {len(item.values):d} "
                f"{mean_all:.12e} {std_all:.12e} {run.start_conf:d} "
                f"{len(post_values):d} {mean_post:.12e} {std_post:.12e} "
                f"{item.note}\n"
            )
    print(f"Saved: {os.path.basename(path)}")


def plot_action_density(series: List[ActionDensitySeries],
                        output_dir: str,
                        gauge_group: str) -> None:
    """Overlay action density versus Monte Carlo sweep for all runs."""
    if not series:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    series_sorted = sorted(series, key=lambda x: (x.run.beta, x.run.boundary, x.run.seed))

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = {
        "periodic": "steelblue",
        "open": "darkorange",
    }
    markers = {
        "periodic": "o",
        "open": "s",
    }

    for item in series_sorted:
        run = item.run
        a = lattice_spacing(run.beta)
        color = colors.get(run.boundary, "black")
        marker = markers.get(run.boundary, "o")
        label = f"$a={a:.4f}$ fm, {run.boundary}"

        ax.plot(item.sweeps, item.values, color=color, lw=0.8, alpha=0.75, label=label)
        step = max(1, len(item.sweeps) // 160)
        ax.scatter(item.sweeps[::step], item.values[::step], color=color, marker=marker,
                   s=8, alpha=0.65, linewidth=0)

        if run.start_conf > 0 and item.sweeps[0] <= run.start_conf <= item.sweeps[-1]:
            ax.axvline(run.start_conf, color=color, linestyle="--", linewidth=1.0, alpha=0.65)

    ax.set_xlabel("MC sweeps")
    ax.set_ylabel(r"$s = S/V$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)

    plt.tight_layout()
    path = os.path.join(output_dir, f"action_density_{gauge_group}.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {os.path.basename(path)}")


def analyze_action_density(runs: List[RunData], output_dir: str, gauge_group: str) -> None:
    series = collect_action_density(runs)
    if not series:
        print("No action-density files found.")
        return
    plot_action_density(series, output_dir, gauge_group)
    write_action_density_summary(series, output_dir, gauge_group)
