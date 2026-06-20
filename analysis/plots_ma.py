#!/usr/bin/env python3
"""Paper/thesis style plots for the selected MA figure set."""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401  Registers SciencePlots styles.
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

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
MC_LABEL = r"$N_{\rm MC}$"
OVERVIEW_NMC_MAX = 800
THERMALIZATION_NMC_MAX = 500
WAVELENGTH_EDGE_SIGMA = 3.0
WAVELENGTH_MIN_CYCLES = 2.0
WAVELENGTH_METHODS = ["peak_spacing", "fft"]
FINE_TARGETS = (0.0299, 0.0191, 0.0107)
FINE_DISPLAY_SIGMAS = np.array([10.0, 20.0, 40.0])
FINE_SIGMA_COLORS = ["#0072B2", "#009E73", "#882255"]


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
    return "PBC" if boundary == "periodic" else "OBC"


def _t160_legend_label(key: tuple) -> str:
    boundary = key[3]
    if boundary == "periodic":
        return "PBC"
    if boundary == "open":
        return "OBC"
    return boundary


def _comparison_style(label: str) -> tuple[str, str, str]:
    if label == "PBC":
        return "*", PBC_COMP_COLOR, r"PBC$_{\rm comp}$"
    if label == "OBC":
        return "^", OBC_COMP_COLOR, r"OBC$_{\rm comp}$"
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
    title = rf"$a={lattice_spacing_su2(beta):.4f}\,\mathrm{{fm}}$"
    title += rf", ${T}\times {L}^3$"
    if exclude > 0:
        title += rf", cut {exclude}"
    return title


def _thermalization_panel_title(key: tuple) -> str:
    beta, T, L, boundary, exclude = key
    title = rf"{_boundary_label(boundary)}, $a={lattice_spacing_su2(beta):.4f}\,\mathrm{{fm}}$"
    title += rf", ${T}\times {L}^3$"
    if exclude > 0:
        title += rf", cut {exclude}"
    return title


def _thermalization_legend_label(key: tuple) -> str:
    beta, _, _, boundary, _ = key
    return rf"{_boundary_label(boundary)}, $a={lattice_spacing_su2(beta):.4f}$"


def plot_susceptibility(results_periodic: list[AnalysisResult],
                        results_open: list[AnalysisResult]) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.5))

    for results, label, marker, color in [
        (results_periodic, "PBC", "o", PBC_COLOR),
        (results_open, "OBC", "s", OBC_COLOR),
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
        (results_periodic, "PBC", "o", PBC_COLOR, FIT_PBC_COLOR, "--"),
        (results_open, "OBC", "s", OBC_COLOR, FIT_OBC_COLOR, "--"),
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
            ax.set_ylabel(r"$N_{\rm MC}$", fontsize=8)
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
        ax.set_xlim(0, OVERVIEW_NMC_MAX)

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
            ax.set_xlabel(MC_LABEL, fontsize=8)
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


def plot_thermalization(runs: list[RunData]) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups(runs)
        if any(len(run.plaq_data) > 0 and len(run.mc_time) > 0 for run in group)
    ]
    if not groups:
        return

    palette = [
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
        "#E69F00", "#882255", "#44AA99", "#AA4499", "#332288",
        "#88CCEE", "#DDCC77", "#117733", "#999933", "#661100",
    ]
    y_values = []
    for _, group in groups:
        for run in group:
            if len(run.plaq_data) == 0 or len(run.mc_time) == 0:
                continue
            mask = run.mc_time <= THERMALIZATION_NMC_MAX
            if np.any(mask):
                y_values.extend(run.plaq_data[mask])
    if not y_values:
        return

    y_min = float(np.min(y_values))
    y_max = float(np.max(y_values))
    y_pad = max(0.01 * (y_max - y_min), 0.002)

    fig, ax = plt.subplots(figsize=(6.9, 3.7))
    style_keys = sorted({
        (_boundary_label(key[3]), round(lattice_spacing_su2(key[0]), 4))
        for key, _ in groups
    }, key=lambda item: (item[1], item[0]))
    color_by_style = {
        style_key: palette[idx % len(palette)]
        for idx, style_key in enumerate(style_keys)
    }
    legend_by_style = {}
    all_starts = set()

    for key, group in groups:
        style_key = (_boundary_label(key[3]), round(lattice_spacing_su2(key[0]), 4))
        color = color_by_style[style_key]
        visible_runs = [
            run for run in group
            if len(run.plaq_data) > 0 and len(run.mc_time) > 0
        ]

        for run_idx, run in enumerate(visible_runs):
            mask = run.mc_time <= THERMALIZATION_NMC_MAX
            if not np.any(mask):
                continue
            alpha = 0.82 if len(visible_runs) == 1 else 0.55 + 0.18 * (run_idx % 2)
            ax.plot(
                run.mc_time[mask],
                run.plaq_data[mask],
                color=color,
                lw=0.9,
                alpha=alpha,
            )
            if run.start_conf <= THERMALIZATION_NMC_MAX:
                all_starts.add(run.start_conf)

        legend_by_style.setdefault(
            style_key,
            (
                plt.Line2D([0], [0], color=color, lw=1.0),
                _thermalization_legend_label(key),
            ),
        )

    for start in sorted(all_starts):
        ax.axvline(start, color="0.25", ls="--", lw=0.75, alpha=0.75)

    ax.set_xlabel(MC_LABEL)
    ax.set_ylabel(r"$\langle P_{\mu\nu}\rangle$")
    ax.set_xlim(0, THERMALIZATION_NMC_MAX)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.grid(True, color="0.88", lw=0.45)
    ax.tick_params(direction="in")

    ax.legend(
        handles=[
            *[item[0] for item in legend_by_style.values()],
            plt.Line2D([0], [0], color="0.25", ls="--", lw=0.75),
        ],
        labels=[*[item[1] for item in legend_by_style.values()], "analysis start"],
        loc="upper right",
        ncol=1,
        fontsize=5.0,
        frameon=False,
        handlelength=1.5,
        borderaxespad=0.4,
    )
    fig.tight_layout()
    _save(fig, "thermalization_su2.png")


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

    ax.set_xlabel(MC_LABEL)
    ax.set_ylabel(r"$Q$")
    ax.set_xlim(-5, max_len + 5)
    ax.set_ylim(q_min - 0.5, q_max + 0.5)
    ax.set_title(r"OBC fine lattices", pad=4)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(False)
    _save(fig, "Q_vs_mctime_obc_a0191_a0107_su2.png")


def plot_open_gaussian_smoothing(runs: list[RunData],
                                 target_a: float,
                                 a_tag: str) -> None:
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
    filename = f"Q_vs_mctime_obc_{a_tag}_gaussian_smoothing_rounded_su2.png"
    q_values = np.concatenate(q_segments)
    sigma_values = np.array([4.0, 10.0, 20.0, 40.0])
    smooth_colors = ["#0072B2", "#D55E00", "#009E73", "#882255"]

    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    offset = 0
    for idx, q_segment in enumerate(q_segments):
        x_segment = offset + np.arange(len(q_segment), dtype=float)
        ax.plot(
            x_segment, q_segment, color="0.45", lw=0.65, alpha=0.24,
            label="original" if idx == 0 else None,
        )
        ax.scatter(x_segment, q_segment, s=5.0, color="0.35", alpha=0.16, zorder=2)

        local_x = np.arange(len(q_segment), dtype=float)
        for sigma, color in zip(sigma_values, smooth_colors):
            y_smooth = np.round(gaussian_smooth(local_x, q_segment, sigma))
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

    ax.set_xlabel(MC_LABEL)
    ax.set_ylabel(r"$Q$")
    ax.set_title(
        rf"OBC, $a={lattice_spacing_su2(key[0]):.4f}\,\mathrm{{fm}}$, rounded Gaussian",
        pad=4,
    )
    if len(q_values) >= 0.75 * OVERVIEW_NMC_MAX:
        ax.set_xlim(0, OVERVIEW_NMC_MAX)
    else:
        ax.set_xlim(-5, len(q_values) + 5)
    ax.set_ylim(q_min - 0.6, q_max + 0.6)
    ax.legend(frameon=False, loc="upper right", ncol=2)
    ax.grid(False)
    _save(fig, filename)


def _detrend_for_period(y: np.ndarray) -> np.ndarray:
    x = np.arange(len(y), dtype=float)
    if len(y) >= 3:
        coeff = np.polyfit(x, y, deg=1)
        y = y - np.polyval(coeff, x)
    y = y - np.mean(y)
    return y


@dataclass(frozen=True)
class WavelengthEstimate:
    wavelength: float
    error: float
    n_features: int
    n_usable: int
    status: str
    amplitude: float = np.nan
    amplitude_error: float = np.nan
    edge_sensitive: bool = False


def _usable_for_wavelength(y_smooth: np.ndarray, sigma: float) -> tuple[np.ndarray, int] | None:
    edge = int(np.ceil(WAVELENGTH_EDGE_SIGMA * sigma))
    n_usable = len(y_smooth) - 2 * edge
    if n_usable < 12:
        return None
    return y_smooth[edge:len(y_smooth) - edge], n_usable


def _missing_wavelength_status(n_samples: int, sigma: float) -> str:
    edge = int(np.ceil(WAVELENGTH_EDGE_SIGMA * sigma))
    return "edge_too_short" if n_samples - 2 * edge < 12 else "no_estimate"


def _quality_checked_estimate(
    wavelength: float,
    error: float,
    n_features: int,
    n_usable: int,
    sigma: float,
    amplitude: float = np.nan,
    amplitude_error: float = np.nan,
    edge_sensitive: bool = False,
) -> WavelengthEstimate:
    if wavelength <= WAVELENGTH_EDGE_SIGMA * sigma:
        status = "lambda_too_close_to_sigma"
    elif n_usable < WAVELENGTH_MIN_CYCLES * wavelength:
        status = "too_few_cycles"
    else:
        status = "ok"
    return WavelengthEstimate(
        float(wavelength),
        float(error),
        int(n_features),
        int(n_usable),
        status,
        float(amplitude),
        float(amplitude_error),
        bool(edge_sensitive),
    )


def _opposite_extrema_amplitudes(extrema: list[tuple[int, int, float]]) -> np.ndarray:
    amplitudes = [
        abs(extrema[idx + 1][2] - extrema[idx][2])
        for idx in range(len(extrema) - 1)
        if extrema[idx + 1][1] != extrema[idx][1]
    ]
    return np.asarray(amplitudes, dtype=float)


def _estimate_peak_spacing(y_smooth: np.ndarray, sigma: float) -> WavelengthEstimate | None:
    edge = int(np.ceil(WAVELENGTH_EDGE_SIGMA * sigma))
    n_usable = len(y_smooth) - 2 * edge
    if n_usable < 12:
        return None
    local_x = np.arange(len(y_smooth), dtype=float)
    y_broad = gaussian_smooth(local_x, y_smooth, max(float(sigma), 1.0))
    y = _detrend_for_period(y_broad)
    amplitude = float(np.ptp(y))
    if len(y) < 6 or amplitude < 1e-8:
        return None

    prominence = max(0.18 * amplitude, 0.05)
    min_distance = max(2, int(round(2.0 * sigma)))
    peak_idx, _ = find_peaks(y, distance=min_distance, prominence=prominence)
    trough_idx, _ = find_peaks(-y, distance=min_distance, prominence=prominence)
    extrema = sorted(
        [(int(idx), 1, float(y[idx])) for idx in peak_idx]
        + [(int(idx), -1, float(y[idx])) for idx in trough_idx]
    )
    edge_sensitive = any(idx < edge or idx >= len(y_smooth) - edge for idx, _, _ in extrema)
    amplitudes = _opposite_extrema_amplitudes(extrema)
    mean_amplitude = float(np.mean(amplitudes)) if len(amplitudes) else np.nan
    amplitude_error = float(np.std(amplitudes, ddof=1)) if len(amplitudes) > 1 else 0.0

    spacings = []
    if len(peak_idx) >= 2:
        spacings.extend(np.diff(peak_idx))
    if len(trough_idx) >= 2:
        spacings.extend(np.diff(trough_idx))
    if spacings:
        spacings = np.asarray(spacings, dtype=float)
        return _quality_checked_estimate(
            float(np.median(spacings)),
            float(np.std(spacings, ddof=1) if len(spacings) > 1 else 0.0),
            len(spacings),
            n_usable,
            sigma,
            mean_amplitude,
            amplitude_error,
            edge_sensitive,
        )

    half_spacings = [
        extrema[idx + 1][0] - extrema[idx][0]
        for idx in range(len(extrema) - 1)
        if extrema[idx + 1][1] != extrema[idx][1]
    ]
    if not half_spacings:
        return None

    spacings = 2.0 * np.asarray(half_spacings, dtype=float)
    return _quality_checked_estimate(
        float(np.median(spacings)),
        float(np.std(spacings, ddof=1) if len(spacings) > 1 else 0.0),
        len(spacings),
        n_usable,
        sigma,
        mean_amplitude,
        amplitude_error,
        edge_sensitive,
    )


def _estimate_autocorr_period(y_smooth: np.ndarray, sigma: float) -> WavelengthEstimate | None:
    usable = _usable_for_wavelength(y_smooth, sigma)
    if usable is None:
        return None
    y_smooth, n_usable = usable
    y = _detrend_for_period(y_smooth)
    if len(y) < 8 or np.std(y) < 1e-8:
        return None

    acf = np.correlate(y, y, mode="full")[len(y) - 1:]
    acf = acf / np.arange(len(y), 0, -1)
    if acf[0] <= 0:
        return None
    acf = acf / acf[0]

    min_lag = max(2, int(round(2.0 * sigma)))
    max_lag = len(y) // 2
    if min_lag >= max_lag:
        return None
    peaks, props = find_peaks(acf[min_lag:max_lag + 1], prominence=0.03, height=0.0)
    if len(peaks) == 0:
        return None

    lags = peaks + min_lag
    best = int(np.argmax(props["prominences"]))
    lag = float(lags[best])
    return _quality_checked_estimate(lag, 0.0, len(peaks), n_usable, sigma)


def _raw_fft_period(y_smooth: np.ndarray, sigma: float) -> tuple[float, int] | None:
    y = _detrend_for_period(y_smooth)
    n = len(y)
    if n < 8 or np.std(y) < 1e-8:
        return None

    window = np.hanning(n)
    spectrum = np.fft.rfft(y * window)
    freqs = np.fft.rfftfreq(n, d=1.0)
    power = np.abs(spectrum) ** 2
    with np.errstate(divide="ignore"):
        periods = 1.0 / freqs

    min_period = max(2.0, 2.0 * sigma)
    mask = (freqs > 0) & (periods >= min_period) & (periods <= n)
    if not np.any(mask):
        return None
    selected = np.where(mask)[0]
    best = selected[int(np.argmax(power[mask]))]
    return float(periods[best]), int(np.sum(mask))


def _fft_spectrum_for_segment(q_segment: np.ndarray, sigma: float) -> tuple[np.ndarray, np.ndarray, float] | None:
    local_x = np.arange(len(q_segment), dtype=float)
    y_smooth = gaussian_smooth(local_x, q_segment.astype(float), sigma)
    usable = _usable_for_wavelength(y_smooth, sigma)
    if usable is None:
        return None

    y_smooth, _ = usable
    y = _detrend_for_period(y_smooth)
    n = len(y)
    if n < 8 or np.std(y) < 1e-8:
        return None

    window = np.hanning(n)
    spectrum = np.fft.rfft(y * window)
    freqs = np.fft.rfftfreq(n, d=1.0)
    power = np.abs(spectrum) ** 2
    with np.errstate(divide="ignore"):
        periods = 1.0 / freqs

    min_period = max(2.0, 2.0 * sigma)
    mask = (freqs > 0) & (periods >= min_period) & (periods <= n)
    if not np.any(mask):
        return None

    periods = periods[mask]
    power = power[mask]
    if np.max(power) <= 0:
        return None
    dominant_period = float(periods[int(np.argmax(power))])
    power = power / np.max(power)
    order = np.argsort(periods)
    return periods[order], power[order], dominant_period


def _estimate_fft_period(y_smooth: np.ndarray, sigma: float) -> WavelengthEstimate | None:
    usable = _usable_for_wavelength(y_smooth, sigma)
    if usable is None:
        return None
    y_smooth, n_usable = usable
    raw = _raw_fft_period(y_smooth, sigma)
    if raw is None:
        return None
    period, n_modes = raw
    return _quality_checked_estimate(period, 0.0, n_modes, n_usable, sigma)


def _sinusoid_with_trend(x: np.ndarray, offset: float, slope: float,
                         amplitude: float, period: float, phase: float) -> np.ndarray:
    x_mid = x - np.mean(x)
    return offset + slope * x_mid + amplitude * np.sin(2.0 * np.pi * x / period + phase)


def _estimate_sinusoid_period(y_smooth: np.ndarray, sigma: float) -> WavelengthEstimate | None:
    usable = _usable_for_wavelength(y_smooth, sigma)
    if usable is None:
        return None
    y_smooth, n_usable = usable
    x = np.arange(len(y_smooth), dtype=float)
    n = len(y_smooth)
    if n < 12 or np.std(y_smooth) < 1e-8:
        return None

    fft_guess = _raw_fft_period(y_smooth, sigma)
    if fft_guess is None:
        return None

    min_period = max(2.0, 2.0 * sigma)
    max_period = float(n)
    period0 = min(max(fft_guess[0], min_period), max_period)
    amplitude0 = max(float(np.std(y_smooth)) * np.sqrt(2.0), 1e-4)
    p0 = [float(np.mean(y_smooth)), 0.0, amplitude0, period0, 0.0]
    lower = [-np.inf, -np.inf, 0.0, min_period, -2.0 * np.pi]
    upper = [np.inf, np.inf, max(float(np.ptp(y_smooth)) * 2.0, amplitude0), max_period, 2.0 * np.pi]

    try:
        popt, pcov = curve_fit(
            _sinusoid_with_trend,
            x,
            y_smooth,
            p0=p0,
            bounds=(lower, upper),
            maxfev=20000,
        )
    except (RuntimeError, ValueError):
        return None

    period = float(popt[3])
    period_err = float(np.sqrt(pcov[3, 3])) if np.isfinite(pcov[3, 3]) else 0.0
    return _quality_checked_estimate(period, period_err, 1, n_usable, sigma)


def _wavelength_estimates_for_segment(q_segment: np.ndarray, sigma: float) -> dict[str, WavelengthEstimate | None]:
    local_x = np.arange(len(q_segment), dtype=float)
    y_smooth = gaussian_smooth(local_x, q_segment.astype(float), sigma)
    return {
        "peak_spacing": _estimate_peak_spacing(y_smooth, sigma),
        "fft": _estimate_fft_period(y_smooth, sigma),
    }


def _sigma_scan_values(q_segments: list[np.ndarray]) -> np.ndarray:
    min_len = min(len(segment) for segment in q_segments)
    if min_len <= 60:
        return np.array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0])
    if min_len <= 120:
        return np.array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 24.0])
    return np.array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0])


def _method_label(method: str) -> str:
    labels = {
        "peak_spacing": r"$\lambda_{\rm peak}$",
        "fft": r"$\lambda_{\rm FFT}$",
    }
    return labels.get(method, method)


def _symmetric_relative_change(left: float, right: float) -> float:
    denominator = abs(left) + abs(right)
    if denominator <= 1e-12:
        return 0.0
    return float(2.0 * abs(right - left) / denominator)


def _scan_sigma_for_wavelengths(runs_open: list[RunData],
                                runs_periodic: list[RunData]) -> None:
    groups = [
        (key, group)
        for key, group in _concat_groups([*runs_periodic, *runs_open])
        if not _is_t160_key(key)
    ]
    methods = WAVELENGTH_METHODS
    colors = {
        "peak_spacing": "#1f77b4",
        "fft": "#2ca02c",
    }

    detailed_rows = []
    median_rows = []
    suggestion_rows = []

    n_cols = 3
    n_rows = (len(groups) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.2, 1.85 * n_rows), sharex=False, sharey=True)
    axes = np.array(axes).reshape(n_rows, n_cols)

    for panel_idx, (key, group) in enumerate(groups):
        ax = axes.flat[panel_idx]
        boundary = _boundary_label(key[3])
        target_a = lattice_spacing_su2(key[0])
        q_segments = [run.Q_rescaled.astype(float) for run in group]
        sigma_values = _sigma_scan_values(q_segments)
        medians_by_method: dict[str, list[tuple[float, float]]] = {method: [] for method in methods}

        for sigma in sigma_values:
            values_by_method: dict[str, list[float]] = {method: [] for method in methods}
            for segment_idx, q_segment in enumerate(q_segments):
                estimates = _wavelength_estimates_for_segment(q_segment, sigma)
                for method in methods:
                    estimate = estimates[method]
                    status = _missing_wavelength_status(len(q_segment), sigma)
                    wavelength = ""
                    amplitude = ""
                    amplitude_error = ""
                    edge_sensitive = ""
                    if estimate is not None:
                        wavelength_value = estimate.wavelength
                        wavelength = f"{wavelength_value:.6g}"
                        status = estimate.status
                        if np.isfinite(estimate.amplitude):
                            amplitude = f"{estimate.amplitude:.6g}"
                            amplitude_error = f"{estimate.amplitude_error:.6g}"
                        edge_sensitive = "yes" if estimate.edge_sensitive else "no"
                        if estimate.status == "ok":
                            values_by_method[method].append(wavelength_value)
                    detailed_rows.append((
                        boundary,
                        target_a,
                        sigma,
                        segment_idx,
                        len(q_segment),
                        method,
                        wavelength,
                        amplitude,
                        amplitude_error,
                        edge_sensitive,
                        status,
                    ))

            for method in methods:
                values = np.array(values_by_method[method], dtype=float)
                if len(values) == 0:
                    continue
                median = float(np.median(values))
                spread = float(np.max(values) - np.min(values)) if len(values) > 1 else 0.0
                medians_by_method[method].append((sigma, median))
                median_rows.append((boundary, target_a, sigma, method, len(values), median, spread))

        for method in methods:
            points = np.array(medians_by_method[method], dtype=float)
            if len(points) == 0:
                continue
            ax.plot(points[:, 0], points[:, 1], marker="o", ms=3.2, lw=1.0,
                    color=colors[method], label=_method_label(method))

        candidates = []
        for sigma in sigma_values:
            local_lambdas = []
            local_scores = []
            for method in methods:
                points = np.array(medians_by_method[method], dtype=float)
                if len(points) < 2 or sigma not in points[:, 0]:
                    continue
                idx = int(np.where(points[:, 0] == sigma)[0][0])
                lam = points[idx, 1]
                if lam < 3.0 * sigma:
                    continue
                deltas = []
                if 0 < idx < len(points) - 1:
                    deltas.append(_symmetric_relative_change(points[idx - 1, 1], lam))
                    deltas.append(_symmetric_relative_change(lam, points[idx + 1, 1]))
                elif idx == 0:
                    deltas.append(_symmetric_relative_change(lam, points[idx + 1, 1]))
                else:
                    deltas.append(_symmetric_relative_change(points[idx - 1, 1], lam))
                local_lambdas.append(lam)
                local_scores.append(float(np.mean(deltas)))
            if len(local_scores) >= 2:
                candidates.append((float(np.mean(local_scores)), float(np.std(local_lambdas)), sigma, len(local_scores)))

        if candidates:
            stable = [candidate for candidate in candidates if candidate[0] <= 0.05]
            if stable:
                score, method_spread, sigma_star, n_methods = min(
                    stable,
                    key=lambda item: (item[2], item[0], item[1]),
                )
                status = "smallest_stable_candidate"
                suggestion_rows.append((boundary, target_a, sigma_star, score, method_spread, n_methods, status))
                ax.axvline(sigma_star, color="0.25", ls="--", lw=0.8, zorder=0)
            else:
                score, method_spread, sigma_star, n_methods = min(
                    candidates,
                    key=lambda item: (item[0], item[1], abs(item[2] - np.median(sigma_values))),
                )
                suggestion_rows.append((boundary, target_a, "", score, method_spread, n_methods, "unstable_sigma"))
        else:
            suggestion_rows.append((boundary, target_a, "", "", "", 0, "no_stable_candidate"))

        ax.set_title(rf"$a={target_a:.4f}\,\mathrm{{fm}}$, {boundary}", pad=3)
        ax.set_xlabel(r"$\sigma$")
        ax.grid(True, color="0.88", lw=0.5)
        ax.tick_params(direction="in")

    for idx in range(len(groups), n_rows * n_cols):
        axes.flat[idx].set_visible(False)
    for ax in axes[:, 0]:
        if ax.get_visible():
            ax.set_ylabel(r"$\lambda$ [$N_{\rm MC}$]")

    handles = [
        plt.Line2D([0], [0], color=colors[method], marker="o", lw=1.0,
                   markersize=3.2, label=_method_label(method))
        for method in methods
    ]
    handles.append(plt.Line2D([0], [0], color="0.25", ls="--", lw=0.8,
                              label=r"suggested $\sigma$"))
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.52, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, "wavelength_sigma_scan_su2.png")

    detailed_path = os.path.join(OUTPUT_DIR, "wavelength_sigma_scan_detailed_su2.csv")
    with open(detailed_path, "w", encoding="utf-8") as handle:
        handle.write(
            "boundary,a_fm,sigma,segment,n_samples,method,wavelength_mc,"
            "amplitude_q,amplitude_std_q,edge_sensitive,status\n"
        )
        for (
            boundary, target_a, sigma, segment_idx, n_samples, method,
            wavelength, amplitude, amplitude_error, edge_sensitive, status,
        ) in detailed_rows:
            handle.write(
                f"{boundary},{target_a:.4f},{sigma:g},{segment_idx},{n_samples},"
                f"{method},{wavelength},{amplitude},{amplitude_error},{edge_sensitive},{status}\n"
            )
    print(f"Saved: {detailed_path}")

    median_path = os.path.join(OUTPUT_DIR, "wavelength_sigma_scan_medians_su2.csv")
    with open(median_path, "w", encoding="utf-8") as handle:
        handle.write("boundary,a_fm,sigma,method,n_segments,median_wavelength_mc,segment_range_mc\n")
        for boundary, target_a, sigma, method, n_values, median, spread in median_rows:
            handle.write(
                f"{boundary},{target_a:.4f},{sigma:g},{method},{n_values},"
                f"{median:.6g},{spread:.6g}\n"
            )
    print(f"Saved: {median_path}")

    suggestion_path = os.path.join(OUTPUT_DIR, "wavelength_sigma_suggestions_su2.csv")
    with open(suggestion_path, "w", encoding="utf-8") as handle:
        handle.write("boundary,a_fm,suggested_sigma,stability_score,method_spread_mc,n_methods,status\n")
        for boundary, target_a, sigma_star, score, method_spread, n_methods, status in suggestion_rows:
            handle.write(
                f"{boundary},{target_a:.4f},{sigma_star},{score},{method_spread},{n_methods},{status}\n"
            )
    print(f"Saved: {suggestion_path}")


def _write_pbc_wavelength_estimates(rows: list[tuple[str, float, float, int, int, str, WavelengthEstimate | None]]) -> None:
    valid_targets = {
        target_a
        for _, target_a, _, _, _, _, estimate in rows
        if estimate is not None and estimate.status == "ok"
    }
    rows = [row for row in rows if row[1] in valid_targets]

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = os.path.join(OUTPUT_DIR, "pbc_smoothed_wavelength_estimates_su2.csv")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            "boundary,a_fm,sigma,segment,n_samples,method,wavelength_mc,error_mc,n_features,"
            "amplitude_q,amplitude_std_q,edge_sensitive,status\n"
        )
        for boundary, target_a, sigma, segment_idx, n_samples, method, estimate in rows:
            if estimate is None:
                status = _missing_wavelength_status(n_samples, sigma)
                handle.write(
                    f"{boundary},{target_a:.4f},{sigma:g},{segment_idx},{n_samples},"
                    f"{method},,,,,,,{status}\n"
                )
                continue
            amplitude = f"{estimate.amplitude:.6g}" if np.isfinite(estimate.amplitude) else ""
            amplitude_error = f"{estimate.amplitude_error:.6g}" if np.isfinite(estimate.amplitude) else ""
            edge_sensitive = "yes" if estimate.edge_sensitive else "no"
            handle.write(
                f"{boundary},{target_a:.4f},{sigma:g},{segment_idx},{n_samples},{method},"
                f"{estimate.wavelength:.6g},{estimate.error:.6g},{estimate.n_features},"
                f"{amplitude},{amplitude_error},{edge_sensitive},{estimate.status}\n"
            )
    print(f"Saved: {path}")


def _valid_wavelength_rows(
    rows: list[tuple[str, float, float, int, int, str, WavelengthEstimate | None]],
) -> list[tuple[str, float, float, int, int, str, WavelengthEstimate | None]]:
    valid_targets = {
        target_a
        for _, target_a, _, _, _, _, estimate in rows
        if estimate is not None
    }
    return [row for row in rows if row[1] in valid_targets]


def _plot_pbc_wavelength_options(
    rows: list[tuple[str, float, float, int, int, str, WavelengthEstimate | None]],
) -> None:
    methods = WAVELENGTH_METHODS
    method_labels = ["peak", "FFT"]
    rows = _valid_wavelength_rows(rows)
    targets = sorted({row[1] for row in rows}, reverse=True)
    boundaries = ["pbc", "obc"]
    colors = {"pbc": PBC_COLOR, "obc": OBC_COLOR}
    markers = {"pbc": "o", "obc": "x"}
    offsets = {"pbc": -0.08, "obc": 0.08}

    fig, axes = plt.subplots(1, len(targets), figsize=(5.4, 3.1), sharey=True)
    if len(targets) == 1:
        axes = [axes]

    for ax, target_a in zip(axes, targets):
        for boundary in boundaries:
            color = colors[boundary]
            marker = markers[boundary]
            offset = offsets[boundary]
            for method_idx, method in enumerate(methods):
                estimates_for_method = [
                    estimate
                    for row_boundary, row_a, _, _, _, row_method, estimate in rows
                    if (
                        row_boundary == boundary
                        and np.isclose(row_a, target_a)
                        and row_method == method
                        and estimate is not None
                    )
                ]
                if len(estimates_for_method) == 0:
                    continue

                x_pos = method_idx + 1 + offset
                jitter = (
                    np.linspace(-0.025, 0.025, len(estimates_for_method))
                    if len(estimates_for_method) > 1 else np.array([0.0])
                )
                for estimate, x_jitter in zip(estimates_for_method, jitter):
                    ax.scatter(
                        x_pos + x_jitter,
                        estimate.wavelength,
                        marker=marker,
                        color=color,
                        s=38 if boundary == "pbc" else 46,
                        linewidths=1.25,
                        alpha=1.0,
                        zorder=3,
                    )

        ax.set_title(rf"$a={target_a:.4f}\,\mathrm{{fm}}$", pad=5)
        ax.set_xticks(np.arange(1, len(methods) + 1))
        ax.set_xticklabels(method_labels)
        ax.set_xlim(0.55, len(methods) + 0.45)
        ax.grid(True, axis="y", color="0.88", lw=0.5)
        ax.tick_params(direction="in", labelbottom=True)

    axes[0].set_ylabel(r"$\lambda$ [method, $N_{\rm MC}$]")
    handles = [
        plt.Line2D([0], [0], color=colors["pbc"], marker=markers["pbc"],
                   lw=0, markersize=5.5, label="PBC"),
        plt.Line2D([0], [0], color=colors["obc"], marker=markers["obc"],
                   lw=0, markersize=6.2, mew=1.35, label="OBC"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.52, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    _save(fig, "wavelengths_su2.png")


def _plot_peak_amplitudes(
    rows: list[tuple[str, float, float, int, int, str, WavelengthEstimate | None]],
) -> None:
    rows = _valid_wavelength_rows(rows)
    amplitude_rows = [
        row for row in rows
        if row[5] == "peak_spacing" and row[6] is not None and np.isfinite(row[6].amplitude)
    ]
    targets = sorted({row[1] for row in amplitude_rows}, reverse=True)
    if not targets:
        return

    colors = {"pbc": PBC_COLOR, "obc": OBC_COLOR}
    markers = {"pbc": "o", "obc": "x"}
    x_pos = {"pbc": 1.0, "obc": 2.0}

    fig, axes = plt.subplots(1, len(targets), figsize=(5.4, 3.1), sharey=True)
    if len(targets) == 1:
        axes = [axes]

    for ax, target_a in zip(axes, targets):
        for boundary in ["pbc", "obc"]:
            estimates = [
                estimate
                for row_boundary, row_a, _, _, _, _, estimate in amplitude_rows
                if row_boundary == boundary and np.isclose(row_a, target_a)
            ]
            if not estimates:
                continue
            jitter = np.linspace(-0.04, 0.04, len(estimates)) if len(estimates) > 1 else np.array([0.0])
            for estimate, x_jitter in zip(estimates, jitter):
                ax.scatter(
                    x_pos[boundary] + x_jitter,
                    estimate.amplitude,
                    marker=markers[boundary],
                    color=colors[boundary],
                    s=38 if boundary == "pbc" else 46,
                    linewidths=1.25,
                    zorder=3,
                )

        ax.set_title(rf"$a={target_a:.4f}\,\mathrm{{fm}}$", pad=5)
        ax.set_xticks([1.0, 2.0])
        ax.set_xticklabels(["PBC", "OBC"])
        ax.set_xlim(0.55, 2.45)
        ax.grid(True, axis="y", color="0.88", lw=0.5)
        ax.tick_params(direction="in", labelbottom=True)

    axes[0].set_ylabel(r"$\langle \Delta Q\rangle_{\rm peak-trough}$")
    fig.tight_layout()
    _save(fig, "amplitudes_su2.png")


def _plot_fft_spectra(
    rows: list[tuple[str, float, float, int, int, str, tuple[np.ndarray, np.ndarray, float] | None]],
) -> None:
    rows = [row for row in rows if row[6] is not None]
    targets = sorted({row[1] for row in rows}, reverse=True)
    if not targets:
        return

    colors = {"pbc": PBC_COLOR, "obc": OBC_COLOR}
    linestyles = {"pbc": "-", "obc": "-"}

    fig, axes = plt.subplots(1, len(targets), figsize=(6.2, 3.1), sharey=True)
    if len(targets) == 1:
        axes = [axes]

    for ax, target_a in zip(axes, targets):
        for boundary in ["pbc", "obc"]:
            spectra = [
                spectrum
                for row_boundary, row_a, _, _, _, _, spectrum in rows
                if row_boundary == boundary and np.isclose(row_a, target_a) and spectrum is not None
            ]
            for spectrum_idx, (periods, power, _) in enumerate(spectra):
                label = boundary.upper() if spectrum_idx == 0 else None
                ax.plot(
                    periods,
                    power,
                    color=colors[boundary],
                    ls=linestyles[boundary],
                    lw=1.0,
                    alpha=0.78,
                    label=label,
                )

        ax.set_title(rf"$a={target_a:.4f}\,\mathrm{{fm}}$", pad=5)
        ax.set_xlabel(r"$\lambda_{\rm FFT}$ [$N_{\rm MC}$]")
        ax.grid(True, color="0.88", lw=0.5)
        ax.tick_params(direction="in")

    axes[0].set_ylabel(r"$P_{\rm FFT}/P_{\rm max}$")
    handles = [
        plt.Line2D([0], [0], color=colors["pbc"], ls=linestyles["pbc"], lw=1.0, label="PBC"),
        plt.Line2D([0], [0], color=colors["obc"], ls=linestyles["obc"], lw=1.0, label="OBC"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.52, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    _save(fig, "fft_spectra_su2.png")


def plot_open_gaussian_vs_periodic_gaussian_rounded(runs_open: list[RunData],
                                                    runs_periodic: list[RunData]) -> None:
    targets = FINE_TARGETS
    row_sigmas = FINE_DISPLAY_SIGMAS
    row_colors = FINE_SIGMA_COLORS

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

    wavelength_rows = []
    fft_spectrum_rows = []
    fig, axes = plt.subplots(len(rows), 2, figsize=(6.9, 6.0), sharey=True)

    for row_idx, (target_a, _, open_segments, _, periodic_segments) in enumerate(rows):
        ax_obc = axes[row_idx, 0]
        ax_pbc = axes[row_idx, 1]
        sigma = row_sigmas[row_idx]
        smooth_color = row_colors[row_idx]
        open_total = sum(len(segment) for segment in open_segments)
        periodic_total = sum(len(segment) for segment in periodic_segments)

        offset = 0
        for segment_idx, q_segment in enumerate(open_segments):
            x_segment = offset + np.arange(len(q_segment), dtype=float)
            ax_obc.plot(x_segment, q_segment, color="0.40", lw=0.55, alpha=0.28)
            ax_obc.scatter(x_segment, q_segment, s=3.5, color="0.30", alpha=0.20, zorder=2)

            local_x = np.arange(len(q_segment), dtype=float)
            q_smooth = np.round(gaussian_smooth(local_x, q_segment, sigma))
            ax_obc.plot(x_segment, q_smooth, color=smooth_color, lw=1.1, zorder=3)

            estimates = _wavelength_estimates_for_segment(q_segment, sigma)
            for method, estimate in estimates.items():
                wavelength_rows.append(("obc", target_a, sigma, segment_idx, len(q_segment), method, estimate))
            fft_spectrum_rows.append((
                "obc",
                target_a,
                sigma,
                segment_idx,
                len(q_segment),
                "fft",
                _fft_spectrum_for_segment(q_segment, sigma),
            ))

            offset += len(q_segment)
            if segment_idx < len(open_segments) - 1:
                ax_obc.axvline(offset - 0.5, color=OBC_COLOR, ls=":", lw=0.75, alpha=0.55)

        offset = 0
        for segment_idx, q_segment in enumerate(periodic_segments):
            x_segment = offset + np.arange(len(q_segment), dtype=float)
            ax_pbc.plot(x_segment, q_segment, color=PBC_COLOR, lw=0.8, alpha=0.42)
            ax_pbc.scatter(x_segment, q_segment, s=3.5, color=PBC_COLOR, alpha=0.28, zorder=2)

            local_x = np.arange(len(q_segment), dtype=float)
            q_smooth = np.round(gaussian_smooth(local_x, q_segment, sigma))
            ax_pbc.plot(x_segment, q_smooth, color=PBC_DARK, lw=1.15, zorder=3)

            estimates = _wavelength_estimates_for_segment(q_segment, sigma)
            for method, estimate in estimates.items():
                wavelength_rows.append(("pbc", target_a, sigma, segment_idx, len(q_segment), method, estimate))
            fft_spectrum_rows.append((
                "pbc",
                target_a,
                sigma,
                segment_idx,
                len(q_segment),
                "fft",
                _fft_spectrum_for_segment(q_segment, sigma),
            ))

            offset += len(q_segment)
            if segment_idx < len(periodic_segments) - 1:
                ax_pbc.axvline(offset - 0.5, color=PBC_COLOR, ls=":", lw=0.75, alpha=0.55)

        for ax in (ax_obc, ax_pbc):
            for q_int in range(q_min, q_max + 1):
                ax.axhline(q_int, color="0.87", ls="--", lw=0.4, zorder=0)
            ax.set_ylim(q_min - 0.5, q_max + 0.5)
            ax.tick_params(direction="in")
            ax.grid(False)
        ax_obc.set_xlim(0, OVERVIEW_NMC_MAX)
        ax_pbc.set_xlim(0, OVERVIEW_NMC_MAX)

        ax_obc.set_ylabel(r"$Q$", fontsize=8)
        ax_obc.text(0.03, 0.90, chr(ord("A") + 2 * row_idx),
                    transform=ax_obc.transAxes, va="top", ha="left", fontweight="bold")
        ax_pbc.text(0.03, 0.90, chr(ord("A") + 2 * row_idx + 1),
                    transform=ax_pbc.transAxes, va="top", ha="left", fontweight="bold")

        title = rf"$a={target_a:.4f}\,\mathrm{{fm}}$, $\sigma={sigma:g}$"
        ax_obc.set_title(title, pad=3)
        ax_pbc.set_title(title, pad=3)

    axes[0, 0].text(0.5, 1.26, "OBC Gaussian rounded", transform=axes[0, 0].transAxes,
                    ha="center", va="bottom")
    axes[0, 1].text(0.5, 1.26, "PBC Gaussian rounded", transform=axes[0, 1].transAxes,
                    ha="center", va="bottom")
    for ax in axes[-1, :]:
        ax.set_xlabel(MC_LABEL, fontsize=8)

    fig.legend(
        handles=[
            plt.Line2D([0], [0], color="0.45", marker="o", lw=0.55,
                       markersize=3.0, alpha=0.55),
            *[
                plt.Line2D([0], [0], color=color, lw=1.1)
                for color in row_colors
            ],
            plt.Line2D([0], [0], color=PBC_COLOR, marker="o", lw=0.75,
                       markersize=3.5, alpha=0.60),
            plt.Line2D([0], [0], color=PBC_DARK, lw=1.15),
        ],
        labels=[
            "OBC original",
            *[rf"$\sigma={sigma:g}$" for sigma in row_sigmas],
            "PBC original",
            "PBC smoothed",
        ],
        loc="upper center",
        ncol=6,
        frameon=False,
        bbox_to_anchor=(0.52, 1.02),
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "Q_vs_mctime_obc_gaussian_pbc_gaussian_rounded_fine_su2.png")
    _write_pbc_wavelength_estimates(wavelength_rows)
    _plot_pbc_wavelength_options(wavelength_rows)
    _plot_peak_amplitudes(wavelength_rows)
    _plot_fft_spectra(fft_spectrum_rows)


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
    max_total_len = 0
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
        max_total_len = max(max_total_len, offset)

    for q_int in range(shared_min, shared_max + 1):
        ax.axhline(q_int, color="0.82", ls="--", lw=0.45, zorder=0)

    ax.set_xlabel(MC_LABEL)
    ax.set_ylabel(r"$Q$")
    ax.set_title(r"$160\times 80^3$, cut 40", pad=4)
    ax.set_ylim(shared_min - 0.5, shared_max + 0.5)
    if max_total_len >= 0.75 * OVERVIEW_NMC_MAX:
        ax.set_xlim(0, OVERVIEW_NMC_MAX)
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
    plot_thermalization(runs_periodic + runs_open)
    plot_q_history_concatenated(runs_open, "open", "Q_vs_mctime_concatenated_su2_open.png")
    plot_q_history_concatenated(runs_periodic, "periodic", "Q_vs_mctime_concatenated_su2_periodic.png")
    plot_open_fine_spacing_overlay(runs_open)
    plot_open_gaussian_smoothing(runs_open, target_a=0.0191, a_tag="a0191")
    _scan_sigma_for_wavelengths(runs_open, runs_periodic)
    plot_open_gaussian_vs_periodic_gaussian_rounded(runs_open, runs_periodic)
    plot_t160_history(runs_periodic + runs_open)


if __name__ == "__main__":
    main()
