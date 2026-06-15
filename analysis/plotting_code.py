#!/usr/bin/env python3
"""
Plotting functions for SU(2)/SU(3) topological charge analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from typing import List

from calculations import (
    RunData, AnalysisResult, gaussian,
    lattice_spacing_su2, lattice_spacing_su3,
    CHI_T_FOURTH_ROOT_SU2, CHI_T_FOURTH_ROOT_SU3,
    find_optimal_alpha, is_t160_l80
)

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'figure.dpi': 150,
})


def _smooth(y: np.ndarray, window: int = None) -> np.ndarray:
    """Rolling mean for a guiding line. Window defaults to ~10% of length, min 3."""
    w = max(3, len(y) // 10) if window is None else window
    return np.convolve(y, np.ones(w) / w, mode='same')


def _is_qtarget_run(run_data: RunData) -> bool:
    return os.path.basename(run_data.run_dir).startswith("qtarget_")


def _highlight_title(title: str, obj) -> str:
    return title + ", T160_L80" if is_t160_l80(obj) else title


def _mark_highlight_panel(ax, obj) -> None:
    if is_t160_l80(obj):
        ax.text(0.98, 0.98, "T160_L80", transform=ax.transAxes,
                fontsize=8, fontweight='bold', color='crimson',
                va='top', ha='right')


def _is_t160_l80_key(key: tuple) -> bool:
    return len(key) >= 3 and int(key[1]) == 160 and int(key[2]) == 80


def _boundary_abbrev(boundary: str) -> str:
    return "pbc" if boundary == "periodic" else "obc" if boundary == "open" else boundary


def _boundary_color(boundary: str) -> str:
    return "steelblue" if boundary == "periodic" else "darkorange" if boundary == "open" else "gray"


def _concat_groups(runs: List[RunData]) -> list[tuple[tuple, list[RunData]]]:
    """Group normal runs whose independent seeds should be concatenated."""
    groups = {}
    for run_data in runs:
        if _is_qtarget_run(run_data):
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
        (key, sorted(group, key=lambda x: x.seed))
        for key, group in sorted(groups.items(), key=lambda item: (item[0][0], item[0][3], item[0][1], item[0][2]))
    ]


def _concat_label(key: tuple, group: List[RunData], gauge_group: str) -> str:
    beta, T, L, boundary, exclude = key
    spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    title = f'$a = {spacing(beta):.4f}$ fm, {T}x{L}^3'
    if _is_t160_l80_key(key):
        title += ', T160_L80'
    if exclude > 0:
        title += f', cut {exclude}'
    title += '\n' + f'{len(group)} seed' + ('' if len(group) == 1 else 's')
    return title


def plot_Q_vs_mctime_grid(runs: List[RunData], output_dir: str, gauge_group: str, boundary: str):
    """Plot Q vs MC time grid for one boundary type."""
    n = len(runs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    n_cols = 2 if n >= 2 else 1
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 3 * n_rows))
    axes = np.array(axes).flatten()

    sorted_runs = sorted(runs, key=lambda x: x.beta)

    # Global y-range across production runs so panels are comparable. qtarget
    # diagnostic runs can have intentionally extreme Q values and should not
    # flatten ordinary ensembles into an apparently integer-valued line.
    range_runs = [r for r in sorted_runs if not os.path.basename(r.run_dir).startswith("qtarget_")]
    if not range_runs:
        range_runs = sorted_runs
    raw_min = min(r.Q_rescaled.min() for r in range_runs)
    raw_max = max(r.Q_rescaled.max() for r in range_runs)
    global_min = int(np.floor(raw_min))
    global_max = int(np.ceil(raw_max))
    pad = 0.5

    for idx, run_data in enumerate(sorted_runs):
        ax = axes[idx]

        Q = run_data.Q_rescaled
        mc_time = np.arange(len(Q))
        panel_min, panel_max = global_min, global_max
        if Q.min() < global_min or Q.max() > global_max:
            panel_min = int(np.floor(Q.min()))
            panel_max = int(np.ceil(Q.max()))

        ax.plot(mc_time, Q.astype(float), color='steelblue', lw=0.8, alpha=0.25, zorder=2)
        ax.scatter(mc_time, Q, s=3, color='steelblue', alpha=0.8, zorder=3)

        for q_int in range(panel_min, panel_max + 1):
            ax.axhline(y=q_int, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        ax.set_ylim(panel_min - pad, panel_max + pad)

        a = lattice_spacing(run_data.beta)
        title = f'$a = {a:.4f}$ fm'
        if run_data.exclude_boundary_slices > 0:
            title += f', cut {run_data.exclude_boundary_slices}'
        title = _highlight_title(title, run_data)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('MC time', fontsize=9)
        ax.set_ylabel(r'$Q_{\rm re}$', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')
        _mark_highlight_panel(ax, run_data)

    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    filepath = os.path.join(output_dir, f"Q_vs_mctime_{gauge_group}_{boundary}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_Q_vs_mctime_concatenated(runs: List[RunData], output_dir: str, gauge_group: str):
    """Plot concatenated-seed Q histories split by boundary condition."""
    groups = _concat_groups(runs)
    if not groups:
        return

    legacy_path = os.path.join(output_dir, f"Q_vs_mctime_concatenated_{gauge_group}.png")
    if os.path.exists(legacy_path):
        os.remove(legacy_path)

    pad = 0.5

    def render_group_set(group_set: list[tuple[tuple, list[RunData]]], filename_template: str) -> None:
        for boundary in ("periodic", "open"):
            filepath = os.path.join(output_dir, filename_template.format(boundary=boundary))
            boundary_groups = [(key, group) for key, group in group_set if key[3] == boundary]
            if not boundary_groups:
                if os.path.exists(filepath):
                    os.remove(filepath)
                continue

            q_arrays = [
                np.concatenate([run.Q_rescaled.astype(float) for run in group])
                for _, group in boundary_groups
            ]
            shared_min = int(np.floor(min(q.min() for q in q_arrays)))
            shared_max = int(np.ceil(max(q.max() for q in q_arrays)))

            n = len(boundary_groups)
            n_cols = 2 if n >= 2 else 1
            n_rows = (n + n_cols - 1) // n_cols

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 3 * n_rows))
            axes = np.array(axes).flatten()

            for idx, ((key, group), q_concat) in enumerate(zip(boundary_groups, q_arrays)):
                ax = axes[idx]
                x = np.arange(len(q_concat))

                ax.plot(x, q_concat, color='steelblue', lw=0.8, alpha=0.3, zorder=2)
                ax.scatter(x, q_concat, s=4, color='steelblue', alpha=0.8, zorder=3)

                offset = 0
                for run_data in group[:-1]:
                    offset += len(run_data.Q_rescaled)
                    ax.axvline(offset - 0.5, color='black', linestyle=':', linewidth=1.0, alpha=0.65)

                for q_int in range(shared_min, shared_max + 1):
                    ax.axhline(y=q_int, color='gray', linestyle='--', linewidth=0.5, alpha=0.45)

                ax.set_ylim(shared_min - pad, shared_max + pad)
                ax.set_title(_concat_label(key, group, gauge_group), fontsize=10)
                ax.set_xlabel('MC sample (concatenated seeds)', fontsize=9)
                ax.set_ylabel(r'$Q_{\rm re}$', fontsize=9)
                ax.tick_params(labelsize=8)
                ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                        fontsize=11, va='top', ha='left')

            for j in range(n, len(axes)):
                axes[j].set_visible(False)

            plt.tight_layout()
            plt.savefig(filepath, dpi=200)
            plt.close()
            print(f"Saved: {os.path.basename(filepath)}")

    def render_t160_group_set(group_set: list[tuple[tuple, list[RunData]]]) -> None:
        for boundary in ("periodic", "open"):
            old_path = os.path.join(output_dir, f"Q_vs_mctime_T160_L80_{gauge_group}_{boundary}.png")
            if os.path.exists(old_path):
                os.remove(old_path)

        filepath = os.path.join(output_dir, f"Q_vs_mctime_T160_L80_{gauge_group}.png")
        if not group_set:
            if os.path.exists(filepath):
                os.remove(filepath)
            return

        ordered_groups = sorted(group_set, key=lambda item: 0 if item[0][3] == "periodic" else 1)
        q_arrays = [
            np.concatenate([run.Q_rescaled.astype(float) for run in group])
            for _, group in ordered_groups
        ]
        shared_min = int(np.floor(min(q.min() for q in q_arrays)))
        shared_max = int(np.ceil(max(q.max() for q in q_arrays)))

        first_key = ordered_groups[0][0]
        beta, T, L = first_key[:3]
        spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
        fig, ax = plt.subplots(figsize=(8, 4))

        for key, group, q_concat in zip([item[0] for item in ordered_groups],
                                        [item[1] for item in ordered_groups],
                                        q_arrays):
            color = _boundary_color(key[3])
            label = f"{_boundary_abbrev(key[3])}, {len(group)} seed" + ('' if len(group) == 1 else 's')
            x = np.arange(len(q_concat))

            ax.plot(x, q_concat, color=color, lw=0.9, alpha=0.45, zorder=2)
            ax.scatter(x, q_concat, s=5, color=color, alpha=0.8, zorder=3, label=label)

            offset = 0
            for run_data in group[:-1]:
                offset += len(run_data.Q_rescaled)
                ax.axvline(offset - 0.5, color=color, linestyle=':', linewidth=1.0, alpha=0.5)

        for q_int in range(shared_min, shared_max + 1):
            ax.axhline(y=q_int, color='gray', linestyle='--', linewidth=0.5, alpha=0.45)

        ax.set_ylim(shared_min - pad, shared_max + pad)
        ax.set_title(f'$a = {spacing(beta):.4f}$ fm, {T}x{L}^3, T160_L80', fontsize=10)
        ax.set_xlabel('MC sample (concatenated seeds)', fontsize=9)
        ax.set_ylabel(r'$Q_{\rm re}$', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(frameon=False, fontsize=9)

        plt.tight_layout()
        plt.savefig(filepath, dpi=200)
        plt.close()
        print(f"Saved: {os.path.basename(filepath)}")

    normal_groups = [(key, group) for key, group in groups if not _is_t160_l80_key(key)]
    t160_groups = [(key, group) for key, group in groups if _is_t160_l80_key(key)]

    render_group_set(normal_groups, f"Q_vs_mctime_concatenated_{gauge_group}_{{boundary}}.png")
    render_t160_group_set(t160_groups)


def plot_histograms_grid(runs: List[RunData], output_dir: str, gauge_group: str, boundary: str):
    """Plot histograms grid for one boundary type."""
    n = len(runs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    n_cols = min(3, n)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows))
    if n == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    sorted_runs = sorted(runs, key=lambda x: x.beta)

    # Keep the shared range representative of production ensembles. qtarget
    # diagnostic runs can have intentionally extreme Q values and should not
    # compress every ordinary panel into a thin line around Q=0.
    range_runs = [r for r in sorted_runs if not os.path.basename(r.run_dir).startswith("qtarget_")]
    if not range_runs:
        range_runs = sorted_runs

    def q_values_for_plot(run_data: RunData) -> np.ndarray:
        return run_data.Q_rescaled.astype(float)

    q_mins = [q_values_for_plot(r).min() for r in range_runs]
    q_maxs = [q_values_for_plot(r).max() for r in range_runs]
    shared_min = int(np.floor(min(q_mins))) - 1
    shared_max = int(np.ceil(max(q_maxs))) + 1

    have_rounded_label = False
    have_continuous_label = False
    have_fit_label = False

    for idx, run_data in enumerate(sorted_runs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]

        Q_cont = run_data.alpha * run_data.Q_raw
        continuous_q = (run_data.boundary == 'open') or (getattr(run_data, "exclude_boundary_slices", 0) > 0)
        Q_plot = q_values_for_plot(run_data)
        panel_min, panel_max = shared_min, shared_max
        if Q_plot.min() < shared_min or Q_plot.max() > shared_max:
            panel_min = int(np.floor(Q_plot.min())) - 1
            panel_max = int(np.ceil(Q_plot.max())) + 1

        int_bins = np.arange(panel_min - 0.5, panel_max + 1.5, 1)
        fine_bins = np.linspace(panel_min - 0.5, panel_max + 0.5, 40)

        if not continuous_q:
            Q_rounded = run_data.Q_rescaled
            rounded_label = None if have_rounded_label else r'$Q_{\rm rounded}$'
            have_rounded_label = True
            counts, bin_edges, _ = ax.hist(Q_rounded, bins=int_bins, density=False,
                color='steelblue', edgecolor='steelblue', linewidth=0.5, alpha=0.3,
                label=rounded_label)

            # Show unrounded alpha*Q_raw as rug ticks. A second count histogram
            # with different bin widths would put incompatible quantities on the
            # same y-axis and distort the apparent Gaussian fit.
            cont_label = None if have_continuous_label else r'$\alpha \hat{Q}$ (unrounded)'
            have_continuous_label = True
            ax.plot(Q_cont, np.full_like(Q_cont, 0.02), '|', color='steelblue',
                    alpha=0.45, markersize=4, transform=ax.get_xaxis_transform(),
                    label=cont_label)

            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            fit_x_min, fit_x_max = panel_min - 1, panel_max + 1
            fit_Q_for_mean = Q_rounded
        else:  # continuous Q, no rounding (open BC or temporal subvolume)
            cont_label = None if have_continuous_label else r'$\alpha \hat{Q}$'
            have_continuous_label = True
            counts, bin_edges, _ = ax.hist(Q_plot, bins=fine_bins, density=False,
                color='steelblue', edgecolor='steelblue', linewidth=1.2, alpha=0.7,
                label=cont_label)

            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            fit_x_min, fit_x_max = panel_min - 1, panel_max + 1
            fit_Q_for_mean = Q_plot

        try:
            sigma0 = max(np.std(fit_Q_for_mean), 0.5)
            occupied_bins = np.sum(counts > 0)
            if len(fit_Q_for_mean) >= 20 and occupied_bins >= 3 and sigma0 > 0.1:
                popt, _ = curve_fit(gaussian, bin_centers, counts,
                    p0=[np.mean(fit_Q_for_mean), sigma0, counts.max()],
                    bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                x_fit = np.linspace(fit_x_min, fit_x_max, 200)
                fit_label = None if have_fit_label else 'Gaussian fit'
                have_fit_label = True
                ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5, label=fit_label)
        except:
            pass

        ax.set_xlim(panel_min - 0.5, panel_max + 0.5)

        a = lattice_spacing(run_data.beta)
        title = f'$a = {a:.4f}$ fm'
        if run_data.exclude_boundary_slices > 0:
            title += f', cut {run_data.exclude_boundary_slices}'
        title = _highlight_title(title, run_data)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r'$Q$', fontsize=9)
        ax.set_ylabel(r'$N$', fontsize=9)
        ax.tick_params(labelsize=8)

        # Add panel label (A, B, C...)
        panel_label = chr(ord('A') + idx)
        ax.text(0.02, 0.98, panel_label, transform=ax.transAxes, fontsize=11,
                va='top', ha='left')
        _mark_highlight_panel(ax, run_data)
    
    for idx in range(n, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)
    
    handles, labels = [], []
    for ax in axes.flatten():
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        for handle, label in zip(ax_handles, ax_labels):
            if label and label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.78, 0.98),
                   frameon=False, fontsize=10)
    
    plt.tight_layout(rect=[0, 0, 0.77, 1])
    filepath = os.path.join(output_dir, f"histograms_{gauge_group}_{boundary}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_histograms_concatenated(runs: List[RunData], output_dir: str, gauge_group: str):
    """Plot concatenated-seed histograms split by boundary condition."""
    groups = _concat_groups(runs)
    if not groups:
        return

    legacy_path = os.path.join(output_dir, f"histograms_concatenated_{gauge_group}.png")
    if os.path.exists(legacy_path):
        os.remove(legacy_path)

    def render_group_set(group_set: list[tuple[tuple, list[RunData]]], filename_template: str) -> None:
        for boundary in ("periodic", "open"):
            filepath = os.path.join(output_dir, filename_template.format(boundary=boundary))
            boundary_groups = [(key, group) for key, group in group_set if key[3] == boundary]
            if not boundary_groups:
                if os.path.exists(filepath):
                    os.remove(filepath)
                continue

            q_arrays = [
                np.concatenate([run.Q_rescaled.astype(float) for run in group])
                for _, group in boundary_groups
            ]
            shared_min = int(np.floor(min(q.min() for q in q_arrays))) - 1
            shared_max = int(np.ceil(max(q.max() for q in q_arrays))) + 1

            n = len(boundary_groups)
            n_cols = min(3, n)
            n_rows = (n + n_cols - 1) // n_cols

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows))
            if n == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)

            have_fit_label = False
            have_q_label = False

            for idx, ((key, group), q_concat) in enumerate(zip(boundary_groups, q_arrays)):
                row, col = idx // n_cols, idx % n_cols
                ax = axes[row, col]

                exclude = key[4]
                continuous_q = (boundary == 'open') or (exclude > 0)
                panel_min, panel_max = shared_min, shared_max

                if continuous_q:
                    bins = np.linspace(panel_min - 0.5, panel_max + 0.5, 40)
                    q_label = None if have_q_label else r'$\alpha \hat{Q}$'
                else:
                    bins = np.arange(panel_min - 0.5, panel_max + 1.5, 1)
                    q_label = None if have_q_label else r'$Q_{\rm rounded}$'
                have_q_label = True

                counts, bin_edges, _ = ax.hist(q_concat, bins=bins, density=False,
                    color='steelblue', edgecolor='steelblue', linewidth=0.7,
                    alpha=0.45 if not continuous_q else 0.7, label=q_label)

                try:
                    sigma0 = max(np.std(q_concat), 0.5)
                    occupied_bins = np.sum(counts > 0)
                    if len(q_concat) >= 20 and occupied_bins >= 3 and sigma0 > 0.1:
                        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                        popt, _ = curve_fit(gaussian, bin_centers, counts,
                            p0=[np.mean(q_concat), sigma0, counts.max()],
                            bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                        x_fit = np.linspace(panel_min - 1, panel_max + 1, 200)
                        fit_label = None if have_fit_label else 'Gaussian fit'
                        have_fit_label = True
                        ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5, label=fit_label)
                except Exception:
                    pass

                ax.set_xlim(panel_min - 0.5, panel_max + 0.5)
                ax.set_title(_concat_label(key, group, gauge_group), fontsize=10)
                ax.set_xlabel(r'$Q$', fontsize=9)
                ax.set_ylabel(r'$N$', fontsize=9)
                ax.tick_params(labelsize=8)
                ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                        fontsize=11, va='top', ha='left')

            for idx in range(n, n_rows * n_cols):
                row, col = idx // n_cols, idx % n_cols
                axes[row, col].set_visible(False)

            handles, labels = [], []
            for ax in axes.flatten():
                ax_handles, ax_labels = ax.get_legend_handles_labels()
                for handle, label in zip(ax_handles, ax_labels):
                    if label and label not in labels:
                        handles.append(handle)
                        labels.append(label)
            if handles:
                fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.78, 0.98),
                           frameon=False, fontsize=10)

            plt.tight_layout(rect=[0, 0, 0.77, 1])
            plt.savefig(filepath, dpi=200)
            plt.close()
            print(f"Saved: {os.path.basename(filepath)}")

    def render_t160_group_set(group_set: list[tuple[tuple, list[RunData]]]) -> None:
        for boundary in ("periodic", "open"):
            old_path = os.path.join(output_dir, f"histograms_T160_L80_{gauge_group}_{boundary}.png")
            if os.path.exists(old_path):
                os.remove(old_path)

        filepath = os.path.join(output_dir, f"histograms_T160_L80_{gauge_group}.png")
        if not group_set:
            if os.path.exists(filepath):
                os.remove(filepath)
            return

        ordered_groups = sorted(group_set, key=lambda item: 0 if item[0][3] == "periodic" else 1)
        q_arrays = [
            np.concatenate([run.Q_rescaled.astype(float) for run in group])
            for _, group in ordered_groups
        ]
        shared_min = int(np.floor(min(q.min() for q in q_arrays))) - 1
        shared_max = int(np.ceil(max(q.max() for q in q_arrays))) + 1

        first_key = ordered_groups[0][0]
        beta, T, L = first_key[:3]
        spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
        fig, ax = plt.subplots(figsize=(6.5, 4))

        for key, group, q_concat in zip([item[0] for item in ordered_groups],
                                        [item[1] for item in ordered_groups],
                                        q_arrays):
            boundary = key[3]
            exclude = key[4]
            continuous_q = (boundary == 'open') or (exclude > 0)
            panel_min, panel_max = shared_min, shared_max
            color = _boundary_color(boundary)
            bc_label = _boundary_abbrev(boundary)

            if continuous_q:
                bins = np.linspace(panel_min - 0.5, panel_max + 0.5, 40)
                q_label = rf'{bc_label} $\alpha \hat{{Q}}$'
            else:
                bins = np.arange(panel_min - 0.5, panel_max + 1.5, 1)
                q_label = rf'{bc_label} $Q_{{\rm rounded}}$'

            counts, bin_edges, _ = ax.hist(q_concat, bins=bins, density=False,
                color=color, edgecolor=color, linewidth=0.8,
                alpha=0.35 if not continuous_q else 0.5, label=q_label)

            try:
                sigma0 = max(np.std(q_concat), 0.5)
                occupied_bins = np.sum(counts > 0)
                if len(q_concat) >= 20 and occupied_bins >= 3 and sigma0 > 0.1:
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    popt, _ = curve_fit(gaussian, bin_centers, counts,
                        p0=[np.mean(q_concat), sigma0, counts.max()],
                        bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                    x_fit = np.linspace(panel_min - 1, panel_max + 1, 200)
                    ax.plot(x_fit, gaussian(x_fit, *popt), color=color,
                            linestyle='-', linewidth=1.5, label=f'{bc_label} Gaussian fit')
            except Exception:
                pass

        ax.set_xlim(shared_min - 0.5, shared_max + 0.5)
        ax.set_title(f'$a = {spacing(beta):.4f}$ fm, {T}x{L}^3, T160_L80', fontsize=10)
        ax.set_xlabel(r'$Q$', fontsize=9)
        ax.set_ylabel(r'$N$', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(frameon=False, fontsize=9)

        plt.tight_layout()
        plt.savefig(filepath, dpi=200)
        plt.close()
        print(f"Saved: {os.path.basename(filepath)}")

    normal_groups = [(key, group) for key, group in groups if not _is_t160_l80_key(key)]
    t160_groups = [(key, group) for key, group in groups if _is_t160_l80_key(key)]

    render_group_set(normal_groups, f"histograms_concatenated_{gauge_group}_{{boundary}}.png")
    render_t160_group_set(t160_groups)


def plot_susceptibility_combined(results_periodic: List[AnalysisResult], 
                                  results_open: List[AnalysisResult],
                                  output_dir: str, gauge_group: str):
    """Plot susceptibility vs lattice spacing."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ref_value = CHI_T_FOURTH_ROOT_SU2 if gauge_group == "su2" else CHI_T_FOURTH_ROOT_SU3
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    point_counter = 1
    have_highlight_label = False
    for results, label, marker, color in [
        (results_periodic, 'Periodic', 'o', 'steelblue'),
        (results_open, 'Open', 's', 'darkorange')
    ]:
        if len(results) == 0:
            continue
        
        betas = np.array([r.beta for r in results])
        chi_fourth_a = np.array([r.chi_t_fourth_root_a for r in results])
        chi_fourth_a_err = np.array([r.chi_t_fourth_root_a_err for r in results])
        a_values = np.array([lattice_spacing(b) for b in betas])

        idx = np.argsort(a_values)[::-1]
        a_values, chi_fourth_a, chi_fourth_a_err = a_values[idx], chi_fourth_a[idx], chi_fourth_a_err[idx]
        results_sorted = [results[i] for i in idx]

        chi_fourth_MeV = chi_fourth_a / a_values
        chi_fourth_MeV_err = chi_fourth_a_err / a_values   # Gamma-method error on χ_t^{1/4}

        highlight = np.array([is_t160_l80(r) for r in results_sorted])
        normal = ~highlight
        if np.any(normal):
            ax.errorbar(a_values[normal], chi_fourth_MeV[normal],
                        yerr=chi_fourth_MeV_err[normal], fmt=marker,
                        markersize=7, linestyle='none', capsize=3,
                        color=color, label=label)
        if np.any(highlight):
            highlight_label = None if have_highlight_label else "T160_L80"
            have_highlight_label = True
            ax.errorbar(a_values[highlight], chi_fourth_MeV[highlight],
                        yerr=chi_fourth_MeV_err[highlight], fmt='^',
                        markersize=8, linestyle='none', capsize=3,
                        color='crimson', label=highlight_label, zorder=5)

        # Add numeric labels to each point
        for i, (x, y) in enumerate(zip(a_values, chi_fourth_MeV)):
            ax.annotate(str(point_counter), (x, y), textcoords='offset points',
                        xytext=(3, 3), fontsize=7, color=color)
            point_counter += 1
    
    ax.axhline(y=ref_value, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    ax.set_xlabel(r'$a$ (fm)')
    ax.set_ylabel(r'$\chi_t^{1/4}$ (MeV)')
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=False)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"susceptibility_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_tau_int_combined(results_periodic: List[AnalysisResult],
                           results_open: List[AnalysisResult],
                           output_dir: str, gauge_group: str):
    """Plot tau_int vs lattice spacing."""
    fig, ax = plt.subplots(figsize=(8, 5))

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    have_highlight_label = False
    for results, label, marker, color in [
        (results_periodic, 'Periodic', 'o', 'steelblue'),
        (results_open, 'Open', 's', 'darkorange')
    ]:
        if len(results) == 0:
            continue

        betas = np.array([r.beta for r in results])
        tau_int = np.array([r.tau_int for r in results])
        dtau_int = np.array([r.dtau_int for r in results])
        a_values = np.array([lattice_spacing(b) for b in betas])

        idx = np.argsort(a_values)[::-1]
        a_values, tau_int, dtau_int = a_values[idx], tau_int[idx], dtau_int[idx]
        results_sorted = [results[i] for i in idx]

        highlight = np.array([is_t160_l80(r) for r in results_sorted])
        normal = ~highlight
        if np.any(normal):
            ax.errorbar(a_values[normal], tau_int[normal], yerr=dtau_int[normal],
                        fmt=marker, markersize=7, linestyle='none',
                        capsize=3, color=color, label=label)
        if np.any(highlight):
            highlight_label = None if have_highlight_label else "T160_L80"
            have_highlight_label = True
            ax.errorbar(a_values[highlight], tau_int[highlight],
                        yerr=dtau_int[highlight], fmt='^', markersize=8,
                        linestyle='none', capsize=3, color='crimson',
                        label=highlight_label, zorder=5)

    ax.axhline(0.5, linestyle='--', color='gray', linewidth=1, alpha=0.7,
               zorder=1, label=r'$\tau_{\rm int}=0.5$')

    ax.set_xlabel(r'$a$ (fm)')
    ax.set_ylabel(r'$\tau_{\rm int}(Q^2)$')
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    filepath = os.path.join(output_dir, f"tau_int_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_thermalization_comparison(runs: List[RunData], output_dir: str, gauge_group: str):
    """Plot all thermalization curves (plaquette vs MC time) with detected thermalization points."""
    runs_with_plaq = [r for r in runs if len(r.plaq_data) > 0]
    if len(runs_with_plaq) == 0:
        return
    
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    cmap = plt.cm.viridis
    colors = [cmap(i / max(1, len(runs_with_plaq) - 1)) for i in range(len(runs_with_plaq))]
    
    runs_sorted = sorted(runs_with_plaq, key=lambda r: r.beta)
    
    for idx, run in enumerate(runs_sorted):
        color = colors[idx]
        a = lattice_spacing(run.beta)
        label = f'$a = {a:.4f}$ fm ({run.boundary[:4]})'
        if is_t160_l80(run):
            label += ', T160_L80'
        lw = 1.5 if is_t160_l80(run) else 0.8
        alpha = 0.95 if is_t160_l80(run) else 0.8
        
        # Plot plaquette vs MC time
        ax.plot(run.mc_time, run.plaq_data, '-', color=color, linewidth=lw,
                alpha=alpha, label=label)
        
        # Mark thermalization point
        if run.start_conf > 0 and run.start_conf < run.mc_time[-1]:
            ax.axvline(x=run.start_conf, color=color, linestyle='--', alpha=0.8, linewidth=1.5)
            marker = '^' if is_t160_l80(run) else 'o'
            size = 70 if is_t160_l80(run) else 50
            marker_color = 'crimson' if is_t160_l80(run) else color
            ax.scatter([run.start_conf], [run.plaq_data[np.searchsorted(run.mc_time, run.start_conf)]],
                       color=marker_color, s=size, zorder=5, marker=marker,
                       edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('MC Sweeps')
    ax.set_ylabel(r'$\langle P_{\mu\nu} \rangle$')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"thermalization_comparison_{gauge_group}.png")
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")
