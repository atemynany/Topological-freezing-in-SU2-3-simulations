"""Helpers for locating synced result directories."""

import glob
import os
from typing import Iterable, List, Optional


def default_results_dirs(base_dir: str) -> List[str]:
    """Default local result roots produced by cluster/sync_results.sh."""
    return [
        os.path.join(base_dir, "data", "results_home"),
        os.path.join(base_dir, "data", "results_work"),
    ]


def resolve_results_dirs(base_dir: str, results_dirs: Optional[Iterable[str]]) -> List[str]:
    """Use explicit result roots, or the home/work sync roots by default."""
    if results_dirs:
        return [os.path.abspath(path) for path in results_dirs]
    return [path for path in default_results_dirs(base_dir) if os.path.isdir(path)]


def find_run_dirs(results_dirs: Iterable[str], gauge_group: str) -> List[str]:
    """Find run directories across result roots without double-counting names."""
    seen = set()
    run_dirs = []

    for results_dir in results_dirs:
        patterns = [
            os.path.join(results_dir, f"T*_L*_b*_*_seed*_{gauge_group}"),
            os.path.join(results_dir, f"qtarget_*_{gauge_group}"),
        ]
        for pattern in patterns:
            for run_dir in sorted(glob.glob(pattern)):
                run_name = os.path.basename(run_dir)
                if run_name in seen:
                    continue
                seen.add(run_name)
                run_dirs.append(run_dir)

    return run_dirs
