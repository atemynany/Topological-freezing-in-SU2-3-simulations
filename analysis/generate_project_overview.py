#!/usr/bin/env python3
"""Generate a static HTML overview of project inputs, results, and figures."""

from __future__ import annotations

import argparse
import datetime as dt
import html
import os
from pathlib import Path
from typing import Iterable

script_dir = Path(__file__).resolve().parent
base_dir = script_dir.parent

import sys

sys.path.insert(0, str(script_dir))

from calculations import analyze_run, load_run_data
from result_paths import find_run_dirs, resolve_results_dirs


KEY_COLUMNS = [
    "T",
    "L",
    "boundary",
    "start_type",
    "num_sweeps",
    "save_interval",
    "beta",
    "smear_steps",
    "smear_interval",
    "smear_alpha",
    "exclude_boundary_slices",
]


def rel(path: Path) -> str:
    return os.path.relpath(path, base_dir)


def href(path: Path, from_file: Path) -> str:
    return os.path.relpath(path, from_file.parent)


def esc(value: object) -> str:
    return html.escape(str(value))


def parse_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                values[parts[0]] = " ".join(parts[1:])
    return values


def first_lines(path: Path, limit: int = 30) -> list[str]:
    lines = []
    with path.open(errors="replace") as f:
        for idx, line in enumerate(f):
            if idx >= limit:
                lines.append("...")
                break
            lines.append(line.rstrip("\n"))
    return lines


def collect_inputs() -> list[dict[str, object]]:
    rows = []
    for path in sorted((base_dir / "input").rglob("*.txt")):
        values = parse_key_values(path)
        rows.append({
            "path": path,
            "group": rel(path.parent),
            "values": values,
        })
    return rows


def collect_outputs() -> dict[str, list[Path]]:
    files = sorted(path for path in (base_dir / "output").rglob("*") if path.is_file())
    grouped: dict[str, list[Path]] = {}
    for path in files:
        group = rel(path.parent)
        grouped.setdefault(group, []).append(path)
    return grouped


def collect_tools() -> list[dict[str, object]]:
    paths = []
    for pattern in ["*.sh", "analysis/*.py", "scripts/*.sh", "cluster/*.sh", "*.md"]:
        paths.extend(base_dir.glob(pattern))

    rows = []
    for path in sorted(set(paths)):
        if not path.is_file():
            continue
        rel_path = rel(path)
        if rel_path.startswith(".git/"):
            continue
        category = rel_path.split(os.sep, 1)[0] if os.sep in rel_path else "project"
        rows.append({
            "path": path,
            "category": category,
            "description": describe_tool(path),
        })
    return rows


def collect_markdown_docs() -> list[dict[str, object]]:
    docs = []
    for path in sorted(base_dir.rglob("*.md")):
        rel_path = rel(path)
        parts = Path(rel_path).parts
        if any(part.startswith(".") for part in parts):
            continue
        if parts[0] in {"build", "cmake-build-debug", "cmake-build-release"}:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        docs.append({
            "path": path,
            "title": markdown_title(text, path),
            "outline": markdown_outline(text),
            "text": text,
        })
    return docs


def markdown_title(text: str, path: Path) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("# "):
            return stripped.lstrip("# ").strip()
    return rel(path)


def markdown_outline(text: str, limit: int = 16) -> list[str]:
    outline = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("## "):
            outline.append(stripped.lstrip("# ").strip())
            if len(outline) >= limit:
                break
    return outline


def describe_tool(path: Path) -> str:
    """Extract a compact human-readable description from a script/doc file."""
    try:
        lines = first_lines(path, limit=24)
    except UnicodeDecodeError:
        return ""

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped == "#!/usr/bin/env python3" or stripped == "#!/bin/bash":
            continue
        if stripped.startswith('"""') and stripped.strip('"'):
            return stripped.strip('" ')
        if stripped.startswith("#") and not set(stripped) <= {"#", " ", "="}:
            return stripped.lstrip("# ").strip()
        if path.suffix == ".md" and stripped.startswith("#"):
            return stripped.lstrip("# ").strip()
    return ""


def classify_result_files(run_dir: Path) -> list[str]:
    names = []
    for name in [
        "topcharge.dat",
        "topcharge_timeslice.dat",
        "plaquette.dat",
        "action_density.dat",
        "action_density_bulk_ex40.dat",
        "input.txt",
        "run.log",
        "run_info.txt",
    ]:
        if list(run_dir.rglob(name)):
            names.append(name)
    return names


def collect_runs(gauge_group: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    result_dirs = resolve_results_dirs(str(base_dir), None)
    run_paths = find_run_dirs(result_dirs, gauge_group)
    inventory = []
    summaries = []

    for run_path in run_paths:
        path = Path(run_path)
        run_data = load_run_data(run_path, gauge_group)
        if run_data is None:
            inventory.append({
                "path": path,
                "status": "not loaded",
                "files": classify_result_files(path),
            })
            continue

        is_qtarget = path.name.startswith("qtarget_")
        inventory.append({
            "path": path,
            "status": "loaded",
            "source": "qtarget hot candidates" if is_qtarget else path.parent.name,
            "beta": run_data.beta,
            "boundary": run_data.boundary,
            "T": run_data.T,
            "L": run_data.L,
            "seed": "combined" if is_qtarget else run_data.seed,
            "n_conf": len(run_data.Q_raw),
            "files": classify_result_files(path),
        })

        try:
            result = analyze_run(run_data)
        except Exception as exc:
            summaries.append({
                "path": path,
                "error": str(exc),
            })
            continue

        summaries.append({
            "path": path,
            "beta": result.beta,
            "boundary": result.boundary,
            "T": run_data.T,
            "L": run_data.L,
            "seed": "combined" if is_qtarget else run_data.seed,
            "n_conf": result.n_conf,
            "Q2": result.Q2_mean,
            "chi": result.chi_t_fourth_root_a,
            "tau": result.tau_int,
            "alpha": result.alpha,
            "source": "qtarget" if is_qtarget else path.parent.name,
        })

    return inventory, summaries


def image_title(path: Path) -> str:
    stem = path.stem.replace("_", " ")
    parent = path.parent.name
    if parent.startswith("smear"):
        return f"{stem} ({parent})"
    return stem


def image_description(path: Path) -> str:
    name = path.stem
    parent = path.parent.name
    if name.startswith("Q_vs_mctime"):
        return "Topological charge history versus Monte Carlo time."
    if name.startswith("histograms"):
        return "Distribution of topological charge values, useful for checking tunnelling and sector coverage."
    if name.startswith("susceptibility"):
        return "Topological susceptibility summary versus lattice spacing."
    if name.startswith("tau_int"):
        return "Integrated autocorrelation time summary versus lattice spacing."
    if name.startswith("thermalization"):
        return "Plaquette evolution used to inspect equilibration and thermalization."
    if name.startswith("action_density"):
        return "Action density versus sweep, including qtarget/hot-list data when available."
    if name.startswith("timeslice_density"):
        return "Topological charge density as a function of time slice."
    if name.startswith("timeslice_mctime"):
        return "Timeslice charge evolution over Monte Carlo time."
    if name.startswith("timeslice_edge_bulk"):
        return "Comparison of boundary/edge and bulk timeslice behavior."
    if name.startswith("timeslice_susceptibility"):
        return "Timeslice-resolved susceptibility estimate."
    if name.startswith("timeslice_tauint"):
        return "Timeslice-resolved autocorrelation estimate."
    if parent.startswith("smear"):
        return f"Per-smear output from {parent}."
    return "Generated analysis figure."


def render_table(headers: Iterable[str], rows: Iterable[Iterable[object]]) -> str:
    header_html = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for row in rows:
        body.append('<tr class="search-row">' + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return f'<table class="interactive-table"><thead><tr>{header_html}</tr></thead><tbody>{"".join(body)}</tbody></table>'


def render_page(gauge_group: str, output_path: Path) -> str:
    docs = collect_markdown_docs()
    inputs = collect_inputs()
    outputs = collect_outputs()
    tools = collect_tools()
    inventory, summaries = collect_runs(gauge_group)
    generated_at = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    pngs = [path for paths in outputs.values() for path in paths if path.suffix.lower() == ".png"]
    text_outputs = [
        path for paths in outputs.values()
        for path in paths
        if path.suffix.lower() in {".txt", ".dat"}
    ]

    input_rows = []
    for item in inputs:
        path = item["path"]
        values = item["values"]
        input_rows.append([
            f'<a href="{esc(href(path, output_path))}">{esc(rel(path))}</a>',
            esc(item["group"]),
            *[esc(values.get(key, "")) for key in KEY_COLUMNS],
        ])

    inventory_rows = []
    for item in inventory:
        path = item["path"]
        inventory_rows.append([
            f'<a href="{esc(href(path, output_path))}">{esc(rel(path))}</a>',
            esc(item.get("source", "")),
            esc(item.get("beta", "")),
            esc(item.get("boundary", "")),
            esc(item.get("T", "")),
            esc(item.get("L", "")),
            esc(item.get("seed", "")),
            esc(item.get("n_conf", "")),
            esc(", ".join(item.get("files", []))),
            esc(item.get("status", "")),
        ])

    summary_rows = []
    for item in sorted(
        summaries,
        key=lambda x: (str(x.get("boundary", "")), float(x.get("beta", 0) or 0), str(x.get("seed", ""))),
    ):
        path = item["path"]
        if "error" in item:
            summary_rows.append([
                f'<a href="{esc(href(path, output_path))}">{esc(rel(path))}</a>',
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                esc(item["error"]),
            ])
            continue
        summary_rows.append([
            f'<a href="{esc(href(path, output_path))}">{esc(Path(path).name)}</a>',
            esc(item["source"]),
            f'{item["beta"]:.2f}',
            esc(item["boundary"]),
            f'{item["T"]}x{item["L"]}^3',
            esc(item["seed"]),
            esc(item["n_conf"]),
            f'{item["Q2"]:.4g}',
            f'{item["tau"]:.3g}',
            f'{item["alpha"]:.3g}',
        ])

    tool_rows = []
    for item in tools:
        path = item["path"]
        tool_rows.append([
            f'<a href="{esc(href(path, output_path))}">{esc(rel(path))}</a>',
            esc(item["category"]),
            esc(item["description"]),
        ])

    doc_cards = []
    for item in docs:
        path = item["path"]
        outline_items = "".join(f"<li>{esc(entry)}</li>" for entry in item["outline"])
        outline = f"<ul>{outline_items}</ul>" if outline_items else "<p>No section outline found.</p>"
        doc_cards.append(
            f"""
            <details class="doc-card search-card" open>
              <summary>
                <span>{esc(item["title"])}</span>
                <a href="{esc(href(path, output_path))}">{esc(rel(path))}</a>
              </summary>
              <div class="doc-outline">{outline}</div>
              <pre class="doc-pre">{esc(item["text"])}</pre>
            </details>
            """
        )

    figure_cards = []
    for path in pngs:
        figure_cards.append(
            f"""
            <figure class="search-card" data-kind="{esc(path.parent.name)}">
              <a href="{esc(href(path, output_path))}" class="figure-link">
                <img src="{esc(href(path, output_path))}" alt="{esc(image_title(path))}">
              </a>
              <figcaption>
                <b>{esc(image_title(path))}</b>
                <span>{esc(image_description(path))}</span>
                <code>{esc(rel(path))}</code>
              </figcaption>
            </figure>
            """
        )

    text_cards = []
    for path in text_outputs:
        preview = "\n".join(first_lines(path, limit=18))
        text_cards.append(
            f"""
            <details class="search-card">
              <summary><a href="{esc(href(path, output_path))}">{esc(rel(path))}</a></summary>
              <pre>{esc(preview)}</pre>
            </details>
            """
        )

    html_page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SU23 Freezing Project Overview</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f7f7f4;
      --panel: #ffffff;
      --ink: #202124;
      --muted: #62676f;
      --line: #d8d9d2;
      --accent: #25636b;
      --accent-soft: #e2eff0;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font: 14px/1.45 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    header {{
      padding: 32px max(28px, 6vw) 24px;
      background: #1e3538;
      color: white;
    }}
    header h1 {{
      margin: 0 0 8px;
      font-size: clamp(28px, 4vw, 44px);
      letter-spacing: 0;
    }}
    header p {{
      max-width: 980px;
      margin: 6px 0 0;
      color: #dce8e9;
    }}
    main {{
      padding: 24px max(18px, 4vw) 48px;
    }}
    nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-bottom: 22px;
    }}
    nav a, .pill {{
      display: inline-block;
      padding: 6px 10px;
      border-radius: 6px;
      background: var(--accent-soft);
      color: var(--accent);
      text-decoration: none;
      font-weight: 600;
    }}
    section {{
      margin: 24px 0;
      padding-top: 8px;
    }}
    h2 {{
      margin: 0 0 12px;
      font-size: 24px;
      letter-spacing: 0;
    }}
    h3 {{
      margin: 18px 0 8px;
      font-size: 17px;
    }}
    .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 16px;
      overflow-x: auto;
    }}
    .toolbar {{
      display: grid;
      grid-template-columns: minmax(220px, 1fr) auto auto;
      gap: 10px;
      align-items: center;
      margin: 14px 0 20px;
      position: sticky;
      top: 0;
      z-index: 5;
      background: color-mix(in srgb, var(--bg) 88%, transparent);
      backdrop-filter: blur(8px);
      padding: 10px 0;
    }}
    .toolbar input, .toolbar select, .toolbar button {{
      border: 1px solid var(--line);
      border-radius: 6px;
      background: white;
      color: var(--ink);
      min-height: 36px;
      padding: 0 10px;
      font: inherit;
    }}
    .toolbar button {{
      cursor: pointer;
      background: var(--accent);
      color: white;
      border-color: var(--accent);
      font-weight: 700;
    }}
    .guide {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
      gap: 12px;
    }}
    .guide article {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
    }}
    .guide h3 {{
      margin-top: 0;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}
    .stat {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
    }}
    .stat b {{
      display: block;
      font-size: 24px;
    }}
    table {{
      border-collapse: collapse;
      min-width: 920px;
      width: 100%;
    }}
    th, td {{
      border-bottom: 1px solid var(--line);
      padding: 7px 8px;
      text-align: left;
      vertical-align: top;
      white-space: nowrap;
    }}
    th {{
      background: #f0f1ec;
      position: sticky;
      top: 0;
      z-index: 1;
      cursor: pointer;
      user-select: none;
    }}
    a {{
      color: var(--accent);
    }}
    code, pre {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    }}
    pre {{
      margin: 10px 0 0;
      padding: 12px;
      overflow-x: auto;
      background: #f4f5f2;
      border-radius: 6px;
      border: 1px solid var(--line);
      max-height: 360px;
    }}
    .figures {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 14px;
    }}
    figure {{
      margin: 0;
      padding: 10px;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
    }}
    figure img {{
      width: 100%;
      height: 220px;
      object-fit: contain;
      background: white;
      border: 1px solid #ecece8;
    }}
    figcaption {{
      margin-top: 8px;
      color: var(--muted);
      font-size: 12px;
      overflow-wrap: anywhere;
      display: grid;
      gap: 4px;
    }}
    figcaption b {{
      color: var(--ink);
      font-size: 13px;
    }}
    figcaption code {{
      color: #455;
      white-space: normal;
    }}
    details {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 10px 12px;
      margin: 8px 0;
    }}
    summary {{
      cursor: pointer;
      font-weight: 600;
    }}
    .note {{
      color: var(--muted);
      max-width: 980px;
    }}
    .doc-card {{
      margin: 12px 0;
      padding: 14px;
    }}
    .doc-card summary {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      align-items: baseline;
      justify-content: space-between;
    }}
    .doc-card summary span {{
      font-size: 17px;
      color: var(--ink);
    }}
    .doc-outline {{
      margin-top: 12px;
      padding: 10px 12px;
      background: #f4f8f8;
      border: 1px solid var(--line);
      border-radius: 6px;
    }}
    .doc-outline ul {{
      margin: 0;
      padding-left: 18px;
      columns: 2;
    }}
    .doc-pre {{
      white-space: pre-wrap;
      max-height: 520px;
    }}
    .background-lead {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
      gap: 12px;
      margin-bottom: 14px;
    }}
    .background-lead article {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
    }}
    .background-lead h3 {{
      margin-top: 0;
    }}
    .hidden-by-filter {{
      display: none !important;
    }}
    .empty-note {{
      display: none;
      color: var(--muted);
      padding: 12px 0;
    }}
    @media (max-width: 760px) {{
      .toolbar {{
        grid-template-columns: 1fr;
        position: static;
      }}
      figure img {{
        height: 180px;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>SU23 Freezing Project Overview</h1>
    <p>Generated at {esc(generated_at)} from the files currently present in this workspace.</p>
    <p>This page summarizes input configurations, synced result roots, loaded analysis runs, generated figures, and text outputs.</p>
  </header>
  <main>
    <nav>
      <a href="#background">Background</a>
      <a href="#workflow">Workflow</a>
      <a href="#guide">Guide</a>
      <a href="#runs">Runs</a>
      <a href="#inputs">Inputs</a>
      <a href="#tools">Tools</a>
      <a href="#figures">Figures</a>
      <a href="#outputs">Text Outputs</a>
      <a href="#regenerate">Regenerate</a>
    </nav>
    <div class="toolbar">
      <input id="globalSearch" type="search" placeholder="Search runs, figures, inputs, scripts...">
      <select id="quickFilter" aria-label="Quick filter">
        <option value="">All content</option>
        <option value="qtarget">qtarget / hot ensemble</option>
        <option value="periodic">periodic</option>
        <option value="open">open</option>
        <option value="b3.15">beta 3.15</option>
        <option value="smear">per-smear</option>
        <option value="action_density">action density</option>
        <option value="timeslice">timeslice</option>
      </select>
      <button type="button" id="clearFilters">Clear</button>
    </div>

    <section id="background">
      <h2>Project Background From Documentation</h2>
      <div class="background-lead">
        <article>
          <h3>Physics Aim</h3>
          <p>This project studies topological freezing in pure-gauge SU(2) and SU(3) lattice gauge theory by measuring topological susceptibility and integrated autocorrelation times as the lattice spacing is reduced.</p>
        </article>
        <article>
          <h3>Core Observable</h3>
          <p>Topological charge is measured after APE smearing with a clover discretisation. The analysis rescales near-integer smeared charges using an alpha estimator before computing <code>&lt;Q^2&gt;</code>, susceptibility, and autocorrelation times.</p>
        </article>
        <article>
          <h3>Boundary Comparison</h3>
          <p>Periodic boundaries preserve topological sectors and can freeze at fine lattice spacing. Open temporal boundaries allow topological charge flow through the boundary, reducing freezing but requiring boundary-slice exclusions.</p>
        </article>
        <article>
          <h3>Production Pipeline</h3>
          <p>The code supports two-phase heatbath plus topcharge scans and a fused SU(2) heatbath/topcharge workflow. Cluster output is synced into <code>results_home</code> and <code>results_work</code>, then analyzed locally.</p>
        </article>
      </div>
      <p class="note">The full Markdown documentation is embedded below. Open or collapse each document as needed; the global search also searches inside these blocks.</p>
      {''.join(doc_cards)}
    </section>

    <section id="workflow">
      <h2>Workflow</h2>
      <div class="grid">
        <div class="stat"><b>{len(docs)}</b>Markdown documentation files</div>
        <div class="stat"><b>{len(inputs)}</b>input text files</div>
        <div class="stat"><b>{len(inventory)}</b>{gauge_group.upper()} analysis runs loaded or discovered</div>
        <div class="stat"><b>{len(pngs)}</b>PNG figures</div>
        <div class="stat"><b>{len(text_outputs)}</b>text or data outputs</div>
        <div class="stat"><b>{len(tools)}</b>scripts and docs indexed</div>
      </div>
      <div class="panel">
        <p>The synced data roots are <code>data/results_home</code> and <code>data/results_work</code>. The sync script copies only <code>.dat</code> files from both cluster locations into these roots.</p>
        <p>Standard runs use directories named like <code>T160_L80_b3.15_open_seed857952899_su2</code>. The qtarget hot ensemble is represented by <code>data/results_work/qtarget_T80_L80_b3.15_periodic_20260527_145728_su2</code>; its <code>hot_candidates/try_*</code> measurements are combined into one analysis ensemble.</p>
        <p>The main analysis loads topological charge, plaquette, topcharge-timeslice, and action-density files where available, then writes plots and summaries under <code>output/figures_analysis</code> or a chosen <code>--output-dir</code>.</p>
      </div>
    </section>

    <section id="guide">
      <h2>What Is Displayed Here?</h2>
      <div class="guide">
        <article>
          <h3>Run Tables</h3>
          <p>The run inventory tells you what data folders were discovered, which files were present, and whether the loader accepted them. The analysis summary recomputes the central observables from the current files.</p>
        </article>
        <article>
          <h3>Topological Charge Plots</h3>
          <p><code>Q_vs_mctime</code> shows Monte Carlo history. <code>histograms</code> shows sector coverage and whether the periodic data clusters near integer topological charge after rescaling.</p>
        </article>
        <article>
          <h3>Continuum Summaries</h3>
          <p><code>susceptibility</code> and <code>tau_int</code> compare ensembles versus lattice spacing. These are the high-level freezing and scaling plots.</p>
        </article>
        <article>
          <h3>Thermalization And Action Density</h3>
          <p><code>thermalization_comparison</code> uses plaquette data. <code>action_density</code> shows action density histories, including the qtarget/hot-list ensemble when present.</p>
        </article>
        <article>
          <h3>Timeslice Analysis</h3>
          <p>The timeslice plots break topological charge density and autocorrelation estimates down along the temporal direction, including edge/bulk comparisons for open boundaries.</p>
        </article>
        <article>
          <h3>Per-Smear Outputs</h3>
          <p>The <code>smear5</code>, <code>smear10</code>, and <code>smear15</code> folders show the same pipeline repeated at fixed APE smearing levels.</p>
        </article>
      </div>
    </section>

    <section id="runs">
      <h2>Run Inventory</h2>
      <div class="panel">
        {render_table(["path", "source", "beta", "boundary", "T", "L", "seed", "N conf", "files found", "status"], inventory_rows)}
      </div>

      <h3>Analysis Summary</h3>
      <p class="note">Values are recomputed from the currently loaded result files when this overview is generated.</p>
      <div class="panel">
        {render_table(["run", "source", "beta", "boundary", "volume", "seed", "N conf", "<Q^2>", "tau_int", "alpha"], summary_rows)}
      </div>
    </section>

    <section id="inputs">
      <h2>Input Configurations</h2>
      <div class="panel">
        {render_table(["file", "group", *KEY_COLUMNS], input_rows)}
      </div>
    </section>

    <section id="tools">
      <h2>Scripts And Project Tools</h2>
      <p class="note">This index links to the local scripts and docs that generate, sync, analyze, and document the project.</p>
      <div class="panel">
        {render_table(["file", "category", "description"], tool_rows)}
      </div>
    </section>

    <section id="figures">
      <h2>Generated Figures</h2>
      <div class="figures">
        {''.join(figure_cards)}
      </div>
    </section>

    <section id="outputs">
      <h2>Text Outputs And Data Previews</h2>
      <p class="note">Each item links to the full file and shows the first few lines.</p>
      {''.join(text_cards)}
    </section>

    <section id="regenerate">
      <h2>Regenerate</h2>
      <div class="panel">
        <h3>Sync cluster results</h3>
        <pre>./cluster/sync_results.sh
./cluster/sync_results.sh --refresh</pre>
        <h3>Run main analysis</h3>
        <pre>MPLCONFIGDIR=/tmp/matplotlib python3 analysis/analysis.py --su2
MPLCONFIGDIR=/tmp/matplotlib python3 analysis/analysis.py --su2 --output-dir output/action_density</pre>
        <h3>Run per-smear analysis</h3>
        <pre>MPLCONFIGDIR=/tmp/matplotlib python3 analysis/analysis_per_smear.py --su2</pre>
        <h3>Rebuild this overview</h3>
        <pre>python3 analysis/generate_project_overview.py</pre>
      </div>
    </section>
  </main>
  <script>
    const searchInput = document.getElementById('globalSearch');
    const quickFilter = document.getElementById('quickFilter');
    const clearFilters = document.getElementById('clearFilters');

    function normalizedFilter() {{
      return [searchInput.value, quickFilter.value]
        .join(' ')
        .trim()
        .toLowerCase();
    }}

    function applyFilters() {{
      const terms = normalizedFilter().split(/\\s+/).filter(Boolean);
      const filterTargets = document.querySelectorAll('.search-row, .search-card');

      filterTargets.forEach((node) => {{
        const haystack = node.textContent.toLowerCase();
        const visible = terms.every((term) => haystack.includes(term));
        node.classList.toggle('hidden-by-filter', !visible);
      }});

      document.querySelectorAll('section').forEach((section) => {{
        const emptyNote = section.querySelector('.empty-note');
        if (!emptyNote) return;
        const visible = section.querySelector('.search-row:not(.hidden-by-filter), .search-card:not(.hidden-by-filter)');
        emptyNote.style.display = visible ? 'none' : 'block';
      }});
    }}

    function sortableValue(text) {{
      const cleaned = text.trim().replace(/,/g, '');
      const number = Number(cleaned);
      return Number.isFinite(number) && cleaned !== '' ? number : cleaned.toLowerCase();
    }}

    document.querySelectorAll('.interactive-table th').forEach((th, columnIndex) => {{
      th.title = 'Click to sort';
      th.addEventListener('click', () => {{
        const table = th.closest('table');
        const body = table.querySelector('tbody');
        const direction = th.dataset.sortDirection === 'asc' ? 'desc' : 'asc';
        table.querySelectorAll('th').forEach((other) => delete other.dataset.sortDirection);
        th.dataset.sortDirection = direction;

        const rows = Array.from(body.querySelectorAll('tr'));
        rows.sort((a, b) => {{
          const av = sortableValue(a.children[columnIndex]?.textContent || '');
          const bv = sortableValue(b.children[columnIndex]?.textContent || '');
          if (av < bv) return direction === 'asc' ? -1 : 1;
          if (av > bv) return direction === 'asc' ? 1 : -1;
          return 0;
        }});
        rows.forEach((row) => body.appendChild(row));
      }});
    }});

    document.querySelectorAll('.figure-link').forEach((link) => {{
      link.addEventListener('click', (event) => {{
        if (!event.altKey) return;
        event.preventDefault();
        const win = window.open('', '_blank');
        if (!win) return;
        win.document.write(`<title>${{link.textContent || 'Figure'}}</title><img src="${{link.href}}" style="max-width:100%;height:auto">`);
      }});
    }});

    searchInput.addEventListener('input', applyFilters);
    quickFilter.addEventListener('change', applyFilters);
    clearFilters.addEventListener('click', () => {{
      searchInput.value = '';
      quickFilter.value = '';
      applyFilters();
    }});
  </script>
</body>
</html>
"""
    return html_page


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gauge", choices=["su2", "su3"], default="su2")
    parser.add_argument("--output", default=str(base_dir / "output" / "project_overview.html"))
    args = parser.parse_args()

    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = base_dir / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_page(args.gauge, output_path), encoding="utf-8")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
