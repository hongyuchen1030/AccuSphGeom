#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path("/global/homes/h/hyvchen/AccuSphGeom")
OUTPUT_DIR = REPO_ROOT / "tests" / "performance_test" / "output"
INPUT_CSV = OUTPUT_DIR / "gca_constlat_saturation_sweep.csv"
SWEEP_PLOT = OUTPUT_DIR / "gca_constlat_saturation_sweep.png"
RATIO_PLOT = OUTPUT_DIR / "gca_constlat_saturation_best_ratio.png"
REPORT_PATH = OUTPUT_DIR / "gca_constlat_saturation_report.txt"


def method_label(method: str, width: int) -> str:
    if method == "pure_fp":
        return "pure_fp" if width == 1 else f"pure_fp_w{width}"
    if method == "accux":
        return "accux" if width == 1 else f"accux_w{width}"
    return f"{method}_w{width}"


def method_color(method: str, width: int) -> str:
    colors = {
        ("pure_fp", 1): "#1f77b4",
        ("pure_fp", 2): "#377eb8",
        ("pure_fp", 4): "#6baed6",
        ("pure_fp", 8): "#9ecae1",
        ("accux", 1): "#2ca02c",
        ("accux", 2): "#4daf4a",
        ("accux", 4): "#8bc34a",
        ("accux", 8): "#b2df8a",
    }
    return colors.get((method, width), "#777777")


def detect_saturation(data_sizes: np.ndarray, times: np.ndarray) -> Optional[int]:
    if len(data_sizes) < 3:
        return None

    rel_changes: List[float] = []
    for idx in range(1, len(times)):
        rel_changes.append(abs(times[idx] - times[idx - 1]) / times[idx])

    for idx in range(1, len(rel_changes)):
        if rel_changes[idx - 1] < 0.05 and rel_changes[idx] < 0.05:
            return int(data_sizes[idx + 1])

    return None


def detect_ratio_stabilization(
    data_sizes: np.ndarray, ratios: np.ndarray
) -> Optional[int]:
    if len(data_sizes) < 3:
        return None

    rel_changes: List[float] = []
    for idx in range(1, len(ratios)):
        rel_changes.append(abs(ratios[idx] - ratios[idx - 1]) / ratios[idx])

    for idx in range(1, len(rel_changes)):
        if rel_changes[idx - 1] < 0.05 and rel_changes[idx] < 0.05:
            return int(data_sizes[idx + 1])

    return None


def main() -> None:
    df = pd.read_csv(INPUT_CSV)

    required = {
        "data_size",
        "method",
        "threadsNum",
        "vec_width",
        "min_time",
        "median_time",
        "mean_time",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df["data_size"] = df["data_size"].astype(int)
    df["method"] = df["method"].astype(str)
    df["threadsNum"] = df["threadsNum"].astype(int)
    df["vec_width"] = df["vec_width"].astype(int)
    df["min_time"] = df["min_time"].astype(float)
    df["median_time"] = df["median_time"].astype(float)
    df["mean_time"] = df["mean_time"].astype(float)

    fig, ax = plt.subplots(figsize=(13, 7))
    for (method, threads_num, vec_width), group in df.groupby(
        ["method", "threadsNum", "vec_width"], sort=True
    ):
        group = group.sort_values("data_size")
        ax.plot(
            group["data_size"],
            group["median_time"],
            marker="o",
            linewidth=1.6,
            markersize=3.5,
            color=method_color(method, vec_width),
            label=f"{method_label(method, vec_width)}_t{threads_num}",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("data_size")
    ax.set_ylabel("median_time (s/point)")
    ax.set_title("GCA constant-latitude saturation sweep")
    ax.grid(True, which="both", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(SWEEP_PLOT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    best_rows: List[Dict[str, float]] = []
    report_lines: List[str] = []
    report_lines.append("GCA constant-latitude saturation report")
    report_lines.append("====================================")
    report_lines.append("")

    for threads_num in sorted(df["threadsNum"].unique()):
        thread_df = df[df["threadsNum"] == threads_num]
        data_sizes = np.array(sorted(thread_df["data_size"].unique()), dtype=int)

        best_pure_df = (
            thread_df[thread_df["method"] == "pure_fp"]
            .sort_values(["data_size", "median_time", "vec_width"])
            .groupby("data_size", as_index=False)
            .first()
        )
        best_accux_df = (
            thread_df[thread_df["method"] == "accux"]
            .sort_values(["data_size", "median_time", "vec_width"])
            .groupby("data_size", as_index=False)
            .first()
        )

        merged = best_pure_df.merge(
            best_accux_df,
            on="data_size",
            suffixes=("_pure", "_accux"),
            how="inner",
        )
        merged["ratio"] = merged["median_time_accux"] / merged["median_time_pure"]

        report_lines.append(f"threads={threads_num}")
        report_lines.append("----------------")
        for _, row in merged.iterrows():
            report_lines.append(
                f"data_size={int(row['data_size'])} "
                f"best_pure_fp=width{int(row['vec_width_pure'])} median={row['median_time_pure']:.6e} "
                f"best_accux=width{int(row['vec_width_accux'])} median={row['median_time_accux']:.6e} "
                f"ratio={row['ratio']:.6f}"
            )
            best_rows.append(
                {
                    "threadsNum": threads_num,
                    "data_size": int(row["data_size"]),
                    "best_pure_width": int(row["vec_width_pure"]),
                    "best_pure_time": float(row["median_time_pure"]),
                    "best_accux_width": int(row["vec_width_accux"]),
                    "best_accux_time": float(row["median_time_accux"]),
                    "ratio": float(row["ratio"]),
                }
            )

        pure_sat = detect_saturation(
            merged["data_size"].to_numpy(), merged["median_time_pure"].to_numpy()
        )
        accux_sat = detect_saturation(
            merged["data_size"].to_numpy(), merged["median_time_accux"].to_numpy()
        )
        ratio_sat = detect_ratio_stabilization(
            merged["data_size"].to_numpy(), merged["ratio"].to_numpy()
        )

        if pure_sat is None:
            report_lines.append("pure_fp best kernel saturation: not detected")
        else:
            report_lines.append(f"pure_fp best kernel saturation: data_size={pure_sat}")

        if accux_sat is None:
            report_lines.append("accux best kernel saturation: not detected")
        else:
            report_lines.append(f"accux best kernel saturation: data_size={accux_sat}")

        if ratio_sat is None:
            report_lines.append("best-kernel ratio stabilization: not detected")
        else:
            report_lines.append(f"best-kernel ratio stabilization: data_size={ratio_sat}")

        ratio_series = merged["ratio"].to_numpy()
        if len(ratio_series) >= 2:
            diffs = np.diff(ratio_series)
            if np.all(diffs <= 0.0):
                report_lines.append("ratio trend: non-increasing with data_size")
            elif np.all(diffs >= 0.0):
                report_lines.append("ratio trend: non-decreasing with data_size")
            else:
                report_lines.append("ratio trend: mixed")
        else:
            report_lines.append("ratio trend: insufficient data")

        if accux_sat is None or ratio_sat is None:
            report_lines.append("thread-count saturation status: unresolved by current criterion")
        else:
            report_lines.append("thread-count saturation status: saturated by current criterion")

        report_lines.append("")

    best_df = pd.DataFrame(best_rows)

    fig, ax = plt.subplots(figsize=(10, 6))
    for threads_num, group in best_df.groupby("threadsNum", sort=True):
        group = group.sort_values("data_size")
        ax.plot(
            group["data_size"],
            group["ratio"],
            marker="o",
            linewidth=1.8,
            markersize=4,
            label=f"threads={threads_num}",
        )

    ax.set_xscale("log")
    ax.set_xlabel("data_size")
    ax.set_ylabel("best_accux / best_pure_fp ratio")
    ax.set_title("Best-kernel ratio vs data_size")
    ax.grid(True, which="both", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(RATIO_PLOT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    if 16 in set(best_df["threadsNum"].unique()):
        thread16 = best_df[best_df["threadsNum"] == 16].sort_values("data_size")
        ratio_sat_16 = detect_ratio_stabilization(
            thread16["data_size"].to_numpy(), thread16["ratio"].to_numpy()
        )
        if ratio_sat_16 is None:
            report_lines.append("16-thread run is still unsaturated: ratio stabilization not detected")
        else:
            report_lines.append(
                f"16-thread run saturation status: ratio stabilized by data_size={ratio_sat_16}"
            )

    REPORT_PATH.write_text("\n".join(report_lines) + "\n", encoding="ascii")

    print(f"Wrote sweep plot to {SWEEP_PLOT}")
    print(f"Wrote best-ratio plot to {RATIO_PLOT}")
    print(f"Wrote saturation report to {REPORT_PATH}")


if __name__ == "__main__":
    main()
