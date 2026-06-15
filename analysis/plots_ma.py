#!/usr/bin/env python3
"""Paper/thesis style plots for the selected MA figure set."""

from __future__ import annotations

import os
import sys
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401  Registers SciencePlots styles.
from scipy.optimize import curve_fit

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

from calculations import (  # noqa: E402
    AnalysisResult,
    RunData,
    CHI_T_FOURTH_ROOT_SU2,
    analyze_run,
    gaussian,
    gaussian_smooth,
    is_t160_l80,
    lattice_spacing_su2,
    load_run_data,
)
from result_paths import find_run_dirs  # noqa: E402
from tauint_fit import fit_power_law, power_law  # noqa: E402


GAUGE_GROUP = "su2"
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results_home")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "plots_MA")

PBC_COLOR = "#4477AA"
PBC_DARK = "#224466"
OBC_COLOR = "#EE7733"
PBC_COMP_COLOR = "#663399"
OBC_COMP_COLOR = "#AA3377"
FIT_PBC_COLOR = "#005AB5"
FIT_OBC_COLOR = "#D55E00"
GRAY = "0.45"


def _apply_style() -> None:
    plt.style.use(["science", "no-latex", "bright"])
    plt.rcParams.update({
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
    })


def _save(fig, filename: str) -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def _collect_data() -> tuple[list[RunData], list[RunData], list[AnalysisResult], list[AnalysisResult]]:
    runs_periodic: list[RunData] = []
    runs_open: list[RunData] = []
    results_periodic: list[AnalysisResult] = []
    results_open: list[AnalysisResult] = []

    for run_dir in find_run_dirs([RESULTS_DIR], GAUGE_GROUP):
        run_data = load_run_data(run_dir, GAUGE_GROUP)
        if run_data is None or run_data.beta > 4:
            continue
        result = analyze_run(run_data)
        if run_data.boundary == "periodic":
            runs_periodic.append(run_data)
            results_periodic.append(result)
        else:
            runs_open.append(run_data)
            results_open.append(result)

    return runs_periodic, runs_open, results_periodic, results_open


def _boundary_color(boundary: str) -> str:
    return PBC_COLOR if boundary == "periodic" else OBC_COLOR


def _boundary_label(boundary: str) -> str:
    return "pbc" if boundary == "periodic" else "obc"


def _t160_legend_label(key: tuple) -> str:
    boundary = key[3]
    if boundary == "periodic":
        return "pbc"
    if boundary == "open":
        return "obc"
    return boundary


def _comparison_style(label: str) -> tuple[str, str, str]:
    if label == "pbc":
        return "*", PBC_COMP_COLOR, r"pbc$_{\rm comp}$"
    if label == "obc":
        return "^", OBC_COMP_COLOR, r"obc$_{\rm comp}$"
    return "^", OBC_COMP_COLOR, rf"${label}_{{\rm comp}}$"


def _concat_groups(runs: Iterable[RunData]) -> list[tuple[tuple, list[RunData]]]:
    groups: dict[tuple, list[RunData]] = {}
    for run_data in runs:
        if os.path.basename(run_data.run_dir).startswith("qtarget_"):
            continue
        key = (
            run_data.beta,
            run_data.T,
            run_data.L,
            run_data.boundary,
            run_data.exclude_boundary_slices,
        )
        groups.setdefault(key, []).append(run_data)
    return [
        (key, sorted(group, key=lambda run: run.seed))
        for key, group in sorted(
            groups.items(),
            key=lambda item: (item[0][0], item[0][3], item[0][1], item[0][2]),
        )
    ]


def _is_t160_key(key: tuple) -> bool:
    return int(key[1]) == 160 and int(key[2]) == 80


def _panel_title(key: tuple, group: list[RunData]) -> str:
    beta, T, L, _, exclude = key
    title = rf"$a={lattice_spacing_su2(beta):.4f}\,\mathrm{{fm}}$, ${T}\times {L}^3$"
    if exclude > 0:
        title += rf", cut {exclude}"
    title += f"\n{len(group)} seed" + ("" if len(group) == 1 else "s")
    return title


def plot_susceptibility(results_periodic: list[AnalysisResult],
                        results_open: list[AnalysisResult]) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.5))

    for results, label, marker, color in [
        (results_periodic, "pbc", "o", PBC_COLOR),
        (results_open, "obc", "s", OBC_COLOR),
    ]:
        ordered = sorted(results, key=lambda r: lattice_spacing_su2(r.beta), reverse=True)
        a_values = np.array([lattice_spacing_su2(r.beta) for r in ordered])
        chi = np.array([r.chi_t_fourth_root_a / lattice_spacing_su2(r.beta) for r in ordered])
        err = np.array([r.chi_t_fourth_root_a_err / lattice_spacing_su2(r.beta) for r in ordered])
        highlight = np.array([is_t160_l80(r) for r in ordered])

        normal = ~highlight
        if np.any(normal):
            ax.errorbar(a_values[normal], chi[normal], yerr=err[normal],
                        fmt=marker, ms=4.2, capsize=2.0, color=color,
                        mfc=color, mec=color, elinewidth=0.8, label=label)
        if np.any(highlight):
            comp_marker, comp_color, comp_label = _comparison_style(label)
            ax.errorbar(a_values[highlight], chi[highlight], yerr=err[highlight],
                        fmt=comp_marker, ms=6.0, capsize=2.0,
                        color=comp_color, mfc=comp_color, mec=comp_color,
                        elinewidth=0.8, label=comp_label, zorder=4)

    ax.axhline(CHI_T_FOURTH_ROOT_SU2, color=GRAY, ls="--", lw=0.9,
               label=rf"$\chi_t^{{1/4}}={CHI_T_FOURTH_ROOT_SU2:.0f}\,\mathrm{{MeV}}$")
    ax.set_xlabel(r"$a$ [fm]")
    ax.set_ylabel(r"$\chi_t^{1/4}$ [MeV]")
    ax.set_xlim(0.008, 0.080)
    ax.set_ylim(-10, 275)
    ax.legend(frameon=False, loc="lower right", handlelength=1.6)
    ax.grid(True, color="0.88", lw=0.5)
    _save(fig, "susceptibility_su2.png")


def plot_tau_int(results_periodic: list[AnalysisResult],
                 results_open: list[AnalysisResult]) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.5))

    all_a = np.array([
        lattice_spacing_su2(r.beta)
        for results in (results_periodic, results_open)
        for r in results
    ])
    for results, label, marker, color, fit_color, line_style in [
        (results_periodic, "pbc", "o", PBC_COLOR, FIT_PBC_COLOR, "--"),
        (results_open, "obc", "s", OBC_COLOR, FIT_OBC_COLOR, "--"),
    ]:
        ordered = sorted(results, key=lambda r: lattice_spacing_su2(r.beta), reverse=True)
        a_values = np.array([lattice_spacing_su2(r.beta) for r in ordered])
        tau = np.array([r.tau_int for r in ordered])
        err = np.array([r.dtau_int for r in ordered])
        highlight = np.array([is_t160_l80(r) for r in ordered])

        normal = ~highlight
        if np.any(normal):
            ax.errorbar(a_values[normal], tau[normal], yerr=err[normal],
                        fmt=marker, ms=4.2, capsize=2.0, color=color,
                        mfc=color, mec=color, elinewidth=0.8, label=label)
        if np.any(highlight):
            comp_marker, comp_color, comp_label = _comparison_style(label)
            ax.errorbar(a_values[highlight], tau[highlight], yerr=err[highlight],
                        fmt=comp_marker, ms=6.0, capsize=2.0,
                        color=comp_color, mfc=comp_color, mec=comp_color,
                        elinewidth=0.8, label=comp_label, zorder=4)

        fit_mask = a_values >= 0.012
        unique_a_desc = np.array(sorted(set(a_values), reverse=True))
        fit_mask &= ~np.isin(a_values, unique_a_desc[:2])
        if np.sum(fit_mask) >= 3:
            fit = fit_power_law(a_values[fit_mask], tau[fit_mask], err[fit_mask],
                                n_boot=4000, seed=42)
            a_grid = np.linspace(all_a.min(), all_a.max(), 400)
            ax.plot(a_grid, power_law(a_grid, fit.c, fit.z), line_style,
                    color=fit_color, lw=1.1,
                    label=rf"{label} fit: $z={fit.z:.2f}$")

    ax.axhline(0.5, color=GRAY, ls=":", lw=0.9, label=r"$\tau_{\rm int}=0.5$")
    ax.set_xlabel(r"$a$ [fm]")
    ax.set_ylabel(r"$\tau_{\rm int}(Q^2)$")
    ax.set_xlim(0.008, 0.080)
    ax.set_ylim(-0.6, 35.5)
    ax.legend(frameon=False, loc="upper right", handlelength=2.0)
    ax.grid(True, color="0.88", lw=0.5)
    _save(fig, "tau_int_su2.png")


def _fit_histogram(ax, q_values: np.ndarray, bins: np.ndarray, color: str) -> None:
    counts, bin_edges, _ = ax.hist(
        q_values,
        bins=bins,
        density=False,
        color=color,
        edgecolor=color,
        alpha=0.42,
        linewidth=0.5,
    )
    try:
        sigma0 = max(np.std(q_values), 0.5)
        occupied_bins = np.sum(counts > 0)
        if len(q_values) >= 20 and occupied_bins >= 3 and sigma0 > 0.1:
            centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            popt, _ = curve_fit(
                gaussian,
                centers,
                counts,
                p0=[np.mean(q_values), sigma0, counts.max()],
                bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]),
            )
            x_fit = np.linspace(bins[0], bins[-1], 220)
            ax.plot(x_fit, gaussian(x_fit, *popt), color=color, lw=1.1)
    except Exception:
        pass


def plot_histograms(runs: list[RunData], boundary: str, filename: str) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if key[3] == boundary and not _is_t160_key(key)
    ]
    q_segments_by_group = [
        [run.Q_rescaled.astype(float) for run in group]
        for _, group in groups
    ]
    q_arrays = [
        np.concatenate(segments)
        for segments in q_segments_by_group
    ]
    q_unrounded_arrays = [
        np.concatenate([(run.alpha * run.Q_raw).astype(float) for run in group])
        for _, group in groups
    ]
    shared_min = int(np.floor(min(q.min() for q in q_arrays))) - 1
    shared_max = int(np.ceil(max(q.max() for q in q_arrays))) + 1

    n_cols = 3
    n_rows = (len(groups) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.6, 4.35), sharex=True)
    axes = np.array(axes).reshape(n_rows, n_cols)
    color = _boundary_color(boundary)

    for idx, ((key, group), q_values, q_unrounded) in enumerate(zip(groups, q_arrays, q_unrounded_arrays)):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]
        continuous_q = boundary == "open" or key[4] > 0
        if continuous_q:
            bins = np.linspace(shared_min - 0.5, shared_max + 0.5, 40)
        else:
            bins = np.arange(shared_min - 0.5, shared_max + 1.5, 1)
        _fit_histogram(ax, q_values, bins, color)
        if boundary == "periodic":
            ax.plot(q_unrounded, np.full_like(q_unrounded, 0.025), "|",
                    color=PBC_DARK, alpha=0.35, markersize=3.0,
                    transform=ax.get_xaxis_transform())
        ax.set_title(_panel_title(key, group), pad=3)
        ax.text(0.04, 0.92, chr(ord("A") + idx), transform=ax.transAxes,
                va="top", ha="left", fontweight="bold")
        ax.grid(False)

    for idx in range(len(groups), n_rows * n_cols):
        axes.flat[idx].set_visible(False)
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel(r"$Q$", fontsize=8)
    for ax in axes[:, 0]:
        if ax.get_visible():
            ax.set_ylabel(r"$N$", fontsize=8)
    for ax in axes.flat:
        if ax.get_visible():
            ax.set_xlim(shared_min - 0.5, shared_max + 0.5)
            ax.tick_params(direction="in")

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=color, ec=color, alpha=0.42),
        plt.Line2D([0], [0], color=color, lw=1.1),
    ]
    labels = [
        r"$Q_{\rm rounded}$" if boundary == "periodic" else r"$\alpha\hat{Q}$",
        "Gaussian fit",
    ]
    if boundary == "periodic":
        handles.append(plt.Line2D([0], [0], color=PBC_DARK, marker="|",
                                  linestyle="none", markersize=5, alpha=0.5))
        labels.append(r"$\alpha\hat{Q}$")
    fig.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        ncol=len(labels),
        frameon=False,
        bbox_to_anchor=(0.52, 1.02),
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, filename)


def plot_q_history_concatenated(runs: list[RunData], boundary: str, filename: str) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if key[3] == boundary and not _is_t160_key(key)
    ]
    q_segments_by_group = [
        [run.Q_rescaled.astype(float) for run in group]
        for _, group in groups
    ]
    q_arrays = [
        np.concatenate(segments)
        for segments in q_segments_by_group
    ]
    shared_min = int(np.floor(min(q.min() for q in q_arrays)))
    shared_max = int(np.ceil(max(q.max() for q in q_arrays)))

    n_cols = 2
    n_rows = (len(groups) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.6, 2.1 * n_rows), sharey=True)
    axes = np.array(axes).reshape(n_rows, n_cols)
    color = _boundary_color(boundary)

    for idx, ((key, group), q_segments) in enumerate(zip(groups, q_segments_by_group)):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]

        offset = 0
        for segment_idx, q_segment in enumerate(q_segments):
            x_segment = offset + np.arange(len(q_segment))
            ax.plot(x_segment, q_segment, color=color, lw=0.75, alpha=0.42)
            ax.scatter(x_segment, q_segment, s=5.5, color=color, alpha=0.82)
            offset += len(q_segment)
            if segment_idx < len(q_segments) - 1:
                ax.axvline(offset - 0.5, color=color, ls=":", lw=0.75, alpha=0.55)

        for q_int in range(shared_min, shared_max + 1):
            ax.axhline(q_int, color="0.86", ls="--", lw=0.4, zorder=0)

        ax.set_title(_panel_title(key, group), pad=3)
        ax.text(0.03, 0.92, chr(ord("A") + idx), transform=ax.transAxes,
                va="top", ha="left", fontweight="bold")
        ax.tick_params(direction="in")
        ax.grid(False)

    for idx in range(len(groups), n_rows * n_cols):
        axes.flat[idx].set_visible(False)
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel("MC sample", fontsize=8)
    for ax in axes[:, 0]:
        if ax.get_visible():
            ax.set_ylabel(r"$Q$", fontsize=8)
    for ax in axes.flat:
        if ax.get_visible():
            ax.set_ylim(shared_min - 0.5, shared_max + 0.5)

    fig.legend(
        handles=[plt.Line2D([0], [0], color=color, marker="o", lw=0.75,
                            markersize=3.5, alpha=0.82)],
        labels=[_boundary_label(boundary)],
        loc="upper center",
        frameon=False,
        bbox_to_anchor=(0.52, 1.01),
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, filename)


def plot_q_history_concatenated_open_gaussian(runs: list[RunData]) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if key[3] == "open" and not _is_t160_key(key)
    ]
    q_segments_by_group = [
        [run.Q_rescaled.astype(float) for run in group]
        for _, group in groups
    ]
    q_arrays = [np.concatenate(segments) for segments in q_segments_by_group]
    shared_min = int(np.floor(min(q.min() for q in q_arrays)))
    shared_max = int(np.ceil(max(q.max() for q in q_arrays)))

    sigma_values = np.array([4.0, 10.0, 20.0, 40.0])
    smooth_colors = ["#0072B2", "#D55E00", "#009E73", "#882255"]
    n_cols = 2
    n_rows = (len(groups) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.6, 2.1 * n_rows), sharey=True)
    axes = np.array(axes).reshape(n_rows, n_cols)

    for idx, ((key, group), q_segments) in enumerate(zip(groups, q_segments_by_group)):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]

        offset = 0
        for segment_idx, q_segment in enumerate(q_segments):
            x_segment = offset + np.arange(len(q_segment), dtype=float)
            ax.plot(x_segment, q_segment, color="0.45", lw=0.65, alpha=0.22)
            ax.scatter(x_segment, q_segment, s=4.0, color="0.35", alpha=0.14, zorder=2)

            local_x = np.arange(len(q_segment), dtype=float)
            for sigma, smooth_color in zip(sigma_values, smooth_colors):
                q_smooth = np.round(gaussian_smooth(local_x, q_segment, sigma))
                ax.plot(x_segment, q_smooth, color=smooth_color, lw=1.0, zorder=3)

            offset += len(q_segment)
            if segment_idx < len(q_segments) - 1:
                ax.axvline(offset - 0.5, color=OBC_COLOR, ls=":", lw=0.75, alpha=0.55)

        for q_int in range(shared_min, shared_max + 1):
            ax.axhline(q_int, color="0.86", ls="--", lw=0.4, zorder=0)

        ax.set_title(rf"$a={lattice_spacing_su2(key[0]):.4f}\,\mathrm{{fm}}$", pad=3)
        ax.text(0.03, 0.92, chr(ord("A") + idx), transform=ax.transAxes,
                va="top", ha="left", fontweight="bold")
        ax.tick_params(direction="in")
        ax.grid(False)

    for idx in range(len(groups), n_rows * n_cols):
        axes.flat[idx].set_visible(False)
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel("MC sample", fontsize=8)
    for ax in axes[:, 0]:
        if ax.get_visible():
            ax.set_ylabel(r"$Q$", fontsize=8)
    for ax in axes.flat:
        if ax.get_visible():
            ax.set_ylim(shared_min - 0.5, shared_max + 0.5)

    fig.legend(
        handles=[
            plt.Line2D([0], [0], color="0.45", marker="o", lw=0.65,
                       markersize=3.0, alpha=0.35),
            *[
                plt.Line2D([0], [0], color=color, lw=1.0)
                for color in smooth_colors
            ],
        ],
        labels=[
            "original",
            *[rf"$\sigma={sigma:g}$" for sigma in sigma_values],
        ],
        loc="upper center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.52, 1.01),
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, "Q_vs_mctime_concatenated_su2_open_gaussian_rounded_sigma40.png")


def plot_open_fine_spacing_overlay(runs: list[RunData]) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if key[3] == "open" and not _is_t160_key(key)
    ]

    selected = []
    for target_a in (0.0191, 0.0107):
        matches = [
            (key, group)
            for key, group in groups
            if abs(lattice_spacing_su2(key[0]) - target_a) < 5e-4
        ]
        if not matches:
            raise RuntimeError(f"missing open ensemble with a ~= {target_a:.4f} fm")
        selected.append(matches[0])

    series = [
        (
            key,
            group,
            [run.Q_rescaled.astype(float) for run in group],
            np.concatenate([run.Q_rescaled.astype(float) for run in group]),
        )
        for key, group in selected
    ]
    q_min = int(np.floor(min(values.min() for _, _, _, values in series)))
    q_max = int(np.ceil(max(values.max() for _, _, _, values in series)))
    max_len = max(len(values) for _, _, _, values in series)

    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    styles = [
        {"color": OBC_COLOR, "alpha_line": 0.28, "alpha_marker": 0.45,
         "markersize": 5.5, "label_suffix": "lighter"},
        {"color": "#C4551E", "alpha_line": 0.72, "alpha_marker": 0.90,
         "markersize": 7.0, "label_suffix": "fine"},
    ]

    for (key, group, q_segments, _), style in zip(series, styles):
        a_value = lattice_spacing_su2(key[0])
        label = rf"$a={a_value:.4f}\,\mathrm{{fm}}$"

        offset = 0
        for segment_idx, q_segment in enumerate(q_segments):
            x_segment = offset + np.arange(len(q_segment))
            ax.plot(x_segment, q_segment, color=style["color"], lw=0.9,
                    alpha=style["alpha_line"], zorder=2)
            ax.scatter(x_segment, q_segment, s=style["markersize"], color=style["color"],
                       alpha=style["alpha_marker"],
                       label=label if segment_idx == 0 else None, zorder=3)
            offset += len(q_segment)
            if segment_idx < len(q_segments) - 1:
                ax.axvline(offset - 0.5, color=style["color"], ls=":",
                           lw=0.8, alpha=0.35)

    for q_int in range(q_min, q_max + 1):
        ax.axhline(q_int, color="0.84", ls="--", lw=0.45, zorder=0)

    ax.set_xlabel("MC sample")
    ax.set_ylabel(r"$Q$")
    ax.set_xlim(-5, max_len + 5)
    ax.set_ylim(q_min - 0.5, q_max + 0.5)
    ax.set_title(r"obc fine-lattice comparison", pad=4)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(False)
    _save(fig, "Q_vs_mctime_obc_a0191_a0107_su2.png")


def plot_open_gaussian_smoothing(runs: list[RunData],
                                 target_a: float,
                                 a_tag: str,
                                 rounded: bool = False) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if key[3] == "open" and not _is_t160_key(key)
    ]
    matches = [
        (key, group)
        for key, group in groups
        if abs(lattice_spacing_su2(key[0]) - target_a) < 5e-4
    ]
    if not matches:
        raise RuntimeError(f"missing open ensemble with a ~= {target_a:.4f} fm")

    key, group = matches[0]
    q_segments = [run.Q_rescaled.astype(float) for run in group]
    if rounded:
        y_label = r"$Q$"
        title_suffix = "Gaussian smoothing, rounded"
        filename = f"Q_vs_mctime_obc_{a_tag}_gaussian_smoothing_rounded_su2.png"
        original_label = "original"
    else:
        y_label = r"$Q$"
        title_suffix = "Gaussian smoothing"
        filename = f"Q_vs_mctime_obc_{a_tag}_gaussian_smoothing_su2.png"
        original_label = "original"
    q_values = np.concatenate(q_segments)
    sigma_values = np.array([4.0, 10.0, 20.0, 40.0])
    smooth_colors = ["#0072B2", "#D55E00", "#009E73", "#882255"]

    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    offset = 0
    for idx, q_segment in enumerate(q_segments):
        x_segment = offset + np.arange(len(q_segment), dtype=float)
        ax.plot(
            x_segment, q_segment, color="0.45", lw=0.65, alpha=0.24,
            label=original_label if idx == 0 else None,
        )
        ax.scatter(x_segment, q_segment, s=5.0, color="0.35", alpha=0.16, zorder=2)

        local_x = np.arange(len(q_segment), dtype=float)
        for sigma, color in zip(sigma_values, smooth_colors):
            y_smooth = gaussian_smooth(local_x, q_segment, sigma)
            if rounded:
                y_smooth = np.round(y_smooth)
            ax.plot(
                x_segment, y_smooth, color=color, lw=1.25,
                label=rf"$\sigma={sigma:g}$" if idx == 0 else None,
                zorder=3,
            )

        offset += len(q_segment)
        if idx < len(q_segments) - 1:
            ax.axvline(offset - 0.5, color=OBC_COLOR, ls=":", lw=0.8, alpha=0.45)

    q_min = int(np.floor(q_values.min()))
    q_max = int(np.ceil(q_values.max()))
    for q_int in range(q_min, q_max + 1):
        ax.axhline(q_int, color="0.88", ls="--", lw=0.4, zorder=0)

    ax.set_xlabel("MC sample")
    ax.set_ylabel(y_label)
    ax.set_title(
        rf"$a={lattice_spacing_su2(key[0]):.4f}\,\mathrm{{fm}}$ obc {title_suffix}",
        pad=4,
    )
    ax.set_xlim(-5, len(q_values) + 5)
    ax.set_ylim(q_min - 0.6, q_max + 0.6)
    ax.legend(frameon=False, loc="upper right", ncol=2)
    ax.grid(False)
    _save(fig, filename)


def plot_open_gaussian_vs_periodic_rounded(runs_open: list[RunData],
                                           runs_periodic: list[RunData]) -> None:
    targets = (0.0299, 0.0191, 0.0107)
    row_sigmas = np.array([4.0, 20.0, 40.0])
    row_colors = ["#0072B2", "#009E73", "#882255"]

    open_groups = [
        (key, group)
        for key, group in _concat_groups(runs_open)
        if key[3] == "open" and not _is_t160_key(key)
    ]
    periodic_groups = [
        (key, group)
        for key, group in _concat_groups(runs_periodic)
        if key[3] == "periodic" and not _is_t160_key(key)
    ]

    rows = []
    for target_a in targets:
        open_match = [
            (key, group)
            for key, group in open_groups
            if abs(lattice_spacing_su2(key[0]) - target_a) < 5e-4
        ]
        periodic_match = [
            (key, group)
            for key, group in periodic_groups
            if abs(lattice_spacing_su2(key[0]) - target_a) < 5e-4
        ]
        if not open_match or not periodic_match:
            raise RuntimeError(f"missing obc or pbc ensemble with a ~= {target_a:.4f} fm")

        open_key, open_group = open_match[0]
        periodic_key, periodic_group = periodic_match[0]
        open_segments = [run.Q_rescaled.astype(float) for run in open_group]
        periodic_segments = [run.Q_rescaled.astype(float) for run in periodic_group]
        rows.append((target_a, open_key, open_segments, periodic_key, periodic_segments))

    all_values = [
        value
        for _, _, open_segments, _, periodic_segments in rows
        for segment in [*open_segments, *periodic_segments]
        for value in segment
    ]
    q_min = int(np.floor(min(all_values)))
    q_max = int(np.ceil(max(all_values)))

    fig, axes = plt.subplots(len(rows), 2, figsize=(6.9, 6.0), sharey=True)

    for row_idx, (target_a, _, open_segments, _, periodic_segments) in enumerate(rows):
        ax_obc = axes[row_idx, 0]
        ax_pbc = axes[row_idx, 1]
        sigma = row_sigmas[row_idx]
        smooth_color = row_colors[row_idx]

        offset = 0
        for segment_idx, q_segment in enumerate(open_segments):
            x_segment = offset + np.arange(len(q_segment), dtype=float)
            ax_obc.plot(x_segment, q_segment, color="0.45", lw=0.55, alpha=0.18)
            ax_obc.scatter(x_segment, q_segment, s=3.5, color="0.35", alpha=0.12, zorder=2)

            local_x = np.arange(len(q_segment), dtype=float)
            q_smooth = np.round(gaussian_smooth(local_x, q_segment, sigma))
            ax_obc.plot(x_segment, q_smooth, color=smooth_color, lw=1.1, zorder=3)

            offset += len(q_segment)
            if segment_idx < len(open_segments) - 1:
                ax_obc.axvline(offset - 0.5, color=OBC_COLOR, ls=":", lw=0.75, alpha=0.55)

        offset = 0
        for segment_idx, q_segment in enumerate(periodic_segments):
            x_segment = offset + np.arange(len(q_segment), dtype=float)
            ax_pbc.plot(x_segment, q_segment, color=PBC_COLOR, lw=0.95, alpha=0.72)

            offset += len(q_segment)
            if segment_idx < len(periodic_segments) - 1:
                ax_pbc.axvline(offset - 0.5, color=PBC_COLOR, ls=":", lw=0.75, alpha=0.55)

        for ax in (ax_obc, ax_pbc):
            for q_int in range(q_min, q_max + 1):
                ax.axhline(q_int, color="0.87", ls="--", lw=0.4, zorder=0)
            ax.set_ylim(q_min - 0.5, q_max + 0.5)
            ax.tick_params(direction="in")
            ax.grid(False)

        ax_obc.set_ylabel(r"$Q$", fontsize=8)
        ax_obc.text(0.03, 0.90, chr(ord("A") + 2 * row_idx),
                    transform=ax_obc.transAxes, va="top", ha="left", fontweight="bold")
        ax_pbc.text(0.03, 0.90, chr(ord("A") + 2 * row_idx + 1),
                    transform=ax_pbc.transAxes, va="top", ha="left", fontweight="bold")

        ax_obc.set_title(rf"$a={target_a:.4f}\,\mathrm{{fm}}$", pad=3)
        ax_pbc.set_title(rf"$a={target_a:.4f}\,\mathrm{{fm}}$", pad=3)

    axes[0, 0].text(0.5, 1.26, "obc Gaussian rounded", transform=axes[0, 0].transAxes,
                    ha="center", va="bottom")
    axes[0, 1].text(0.5, 1.26, "pbc rounded", transform=axes[0, 1].transAxes,
                    ha="center", va="bottom")
    for ax in axes[-1, :]:
        ax.set_xlabel("MC sample", fontsize=8)

    fig.legend(
        handles=[
            plt.Line2D([0], [0], color="0.45", marker="o", lw=0.55,
                       markersize=3.0, alpha=0.35),
            *[
                plt.Line2D([0], [0], color=color, lw=1.1)
                for color in row_colors
            ],
            plt.Line2D([0], [0], color=PBC_COLOR, marker="o", lw=0.75,
                       markersize=3.5, alpha=0.82),
        ],
        labels=[
            "obc original",
            *[rf"$\sigma={sigma:g}$" for sigma in row_sigmas],
            "pbc line",
        ],
        loc="upper center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.52, 1.02),
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "Q_vs_mctime_obc_gaussian_pbc_rounded_fine_su2.png")


def plot_t160_history(runs: list[RunData]) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if _is_t160_key(key)
    ]
    groups = sorted(groups, key=lambda item: 0 if item[0][3] == "periodic" else 1)
    q_segments_by_group = [
        [run.Q_rescaled.astype(float) for run in group]
        for _, group in groups
    ]
    q_arrays = [np.concatenate(segments) for segments in q_segments_by_group]
    shared_min = int(np.floor(min(q.min() for q in q_arrays)))
    shared_max = int(np.ceil(max(q.max() for q in q_arrays)))

    sigma = 40.0
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    for key, group, q_segments in zip([item[0] for item in groups],
                                      [item[1] for item in groups],
                                      q_segments_by_group):
        color = _boundary_color(key[3])
        label = _t160_legend_label(key)
        smooth_label = rf"{label}, rounded Gaussian $\sigma={sigma:g}$"
        raw_line_alpha = 0.26 if key[3] == "open" else 0.50
        raw_marker_alpha = 0.42 if key[3] == "open" else 0.85
        smooth_alpha = 0.72 if key[3] == "open" else 1.0
        offset = 0
        for segment_idx, q_segment in enumerate(q_segments):
            x_segment = offset + np.arange(len(q_segment), dtype=float)
            ax.plot(x_segment, q_segment, color=color, lw=0.9, alpha=raw_line_alpha)
            ax.scatter(x_segment, q_segment, s=9, color=color, alpha=raw_marker_alpha,
                       label=label if segment_idx == 0 else None)

            local_x = np.arange(len(q_segment), dtype=float)
            q_smooth = np.round(gaussian_smooth(local_x, q_segment, sigma))
            ax.plot(
                x_segment, q_smooth, color=color, lw=1.5, ls="--",
                alpha=smooth_alpha,
                label=smooth_label if segment_idx == 0 else None,
                zorder=4,
            )

            offset += len(q_segment)
            if segment_idx < len(q_segments) - 1:
                ax.axvline(offset - 0.5, color=color, ls=":", lw=0.8, alpha=0.55)

    for q_int in range(shared_min, shared_max + 1):
        ax.axhline(q_int, color="0.82", ls="--", lw=0.45, zorder=0)

    ax.set_xlabel(r"MC sample")
    ax.set_ylabel(r"$Q$")
    ax.set_title(r"$160\times 80^3$ boundary-condition comparison, cut 40", pad=4)
    ax.set_ylim(shared_min - 0.5, shared_max + 0.5)
    ax.legend(frameon=False, loc="upper right", handlelength=1.6)
    ax.grid(False)
    _save(fig, "Q_vs_mctime_T160_L80_su2.png")


def main() -> None:
    _apply_style()
    runs_periodic, runs_open, results_periodic, results_open = _collect_data()
    plot_tau_int(results_periodic, results_open)
    plot_susceptibility(results_periodic, results_open)
    plot_histograms(runs_open, "open", "histograms_concatenated_su2_open.png")
    plot_histograms(runs_periodic, "periodic", "histograms_concatenated_su2_periodic.png")
    plot_q_history_concatenated(runs_open, "open", "Q_vs_mctime_concatenated_su2_open.png")
    plot_q_history_concatenated(runs_periodic, "periodic", "Q_vs_mctime_concatenated_su2_periodic.png")
    plot_q_history_concatenated_open_gaussian(runs_open)
    plot_open_fine_spacing_overlay(runs_open)
    plot_open_gaussian_smoothing(runs_open, target_a=0.0191, a_tag="a0191")
    plot_open_gaussian_smoothing(runs_open, target_a=0.0191, a_tag="a0191", rounded=True)
    plot_open_gaussian_smoothing(runs_open, target_a=0.0107, a_tag="a0107")
    plot_open_gaussian_smoothing(runs_open, target_a=0.0107, a_tag="a0107", rounded=True)
    plot_open_gaussian_vs_periodic_rounded(runs_open, runs_periodic)
    plot_t160_history(runs_periodic + runs_open)


if __name__ == "__main__":
    main()
