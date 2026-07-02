#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_SUMMARY_CSV = Path(
    "tests/performance_test/output/gca_constlat_SIMDPack_timing_summary.csv"
)
DEFAULT_REPEATS_CSV = Path(
    "tests/performance_test/output/gca_constlat_SIMDPack_timing_repeats.csv"
)
DEFAULT_OLD_CSV = Path(
    "tests/performance_test/output/gca_constlat_SIMDPack_timing.csv"
)


def method_label(method: str, width: int) -> str:
    if method == "pure_fp":
        if width == 1:
            return "FP64 Full API Scalar"
        return f"FP64 Full API SIMD Vector Width {width}"

    if method == "accux":
        if width == 1:
            return "AccuX Constant-Lat Full API Scalar"
        return f"AccuX Constant-Lat Full API SIMD Vector Width {width}"

    return f"{method} width {width}"


def method_color(method: str, width: int):
    # Green shades for AccuX, blue shades for FP64.
    # This intentionally uses fixed colors to mimic the paper-style figure.
    colors = {
        ("accux", 1): "#2ca02c",
        ("accux", 2): "#4daf4a",
        ("accux", 4): "#8bc34a",
        ("accux", 8): "#b2df8a",
        ("pure_fp", 1): "#1f77b4",
        ("pure_fp", 2): "#377eb8",
        ("pure_fp", 4): "#6baed6",
        ("pure_fp", 8): "#9ecae1",
    }
    return colors.get((method, width), "#777777")


def resolve_csv_path(csv_arg: Optional[str]) -> Path:
    if csv_arg is not None:
        csv_path = Path(csv_arg)

        preferred_paths = []
        if csv_path.name == DEFAULT_OLD_CSV.name:
            preferred_paths.extend(
                [
                    csv_path.with_name(DEFAULT_SUMMARY_CSV.name),
                    csv_path.with_name(DEFAULT_REPEATS_CSV.name),
                ]
            )

        for preferred_path in preferred_paths:
            if preferred_path.exists():
                return preferred_path

        if csv_path.exists():
            return csv_path

        raise FileNotFoundError(f"Timing CSV not found: {csv_path}")

    for candidate in (DEFAULT_SUMMARY_CSV, DEFAULT_REPEATS_CSV, DEFAULT_OLD_CSV):
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "No timing CSV found. Looked for summary, repeats, and legacy timing CSVs."
    )


def normalize_timing_frame(df: pd.DataFrame) -> pd.DataFrame:
    required_base = {"method", "threadsNum", "vec_width"}
    missing_base = required_base - set(df.columns)
    if missing_base:
        raise ValueError(f"Missing required columns in CSV: {sorted(missing_base)}")

    if "median_time" in df.columns:
        df = df.copy()
        df["time"] = df["median_time"].astype(float)
        return df

    if "time" in df.columns:
        return df

    raise ValueError("CSV must contain either 'median_time' or 'time'.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        default=None,
        help=(
            "Input timing CSV. Defaults to the summary CSV, then repeats CSV, "
            "then the legacy timing CSV."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output figure path, usually .png or .pdf.",
    )
    parser.add_argument(
        "--title",
        default="GCA Constant-Latitude Full-API Performance",
        help="Plot title.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output DPI for PNG.",
    )
    args = parser.parse_args()

    csv_path = resolve_csv_path(args.csv)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = normalize_timing_frame(pd.read_csv(csv_path))

    df["method"] = df["method"].astype(str)
    df["threadsNum"] = df["threadsNum"].astype(int)
    df["vec_width"] = df["vec_width"].astype(int)
    df["time"] = df["time"].astype(float)

    threads = sorted(df["threadsNum"].unique())

    ordered_series = [
        ("accux", 1),
        ("accux", 2),
        ("accux", 4),
        ("accux", 8),
        ("pure_fp", 1),
        ("pure_fp", 2),
        ("pure_fp", 4),
        ("pure_fp", 8),
    ]

    available_series = []
    for method, width in ordered_series:
        sub = df[(df["method"] == method) & (df["vec_width"] == width)]
        if not sub.empty:
            available_series.append((method, width))

    if not available_series:
        raise ValueError("No known pure_fp or accux rows found in CSV.")

    pivot = {}
    for method, width in available_series:
        sub = df[(df["method"] == method) & (df["vec_width"] == width)]
        value_by_thread = dict(zip(sub["threadsNum"], sub["time"]))
        pivot[(method, width)] = [value_by_thread.get(t, np.nan) for t in threads]

    # Ratio labels are relative to the fastest available result within each thread group.
    # This mirrors the paper-style "x" annotations without using accuracy metrics.
    best_by_thread = {}
    for j, t in enumerate(threads):
        values = [
            pivot[key][j]
            for key in available_series
            if not np.isnan(pivot[key][j]) and pivot[key][j] > 0.0
        ]
        best_by_thread[t] = min(values) if values else np.nan

    fig, ax = plt.subplots(figsize=(12, 5.8))

    x = np.arange(len(threads))
    group_width = 0.82
    bar_width = group_width / len(available_series)
    start_offset = -0.5 * group_width + 0.5 * bar_width

    positive_times = df[df["time"] > 0.0]["time"].to_numpy()
    ymin = positive_times.min() * 0.45
    ymax = positive_times.max() * 4.0

    for k, (method, width) in enumerate(available_series):
        times = np.asarray(pivot[(method, width)], dtype=float)
        xpos = x + start_offset + k * bar_width

        bars = ax.bar(
            xpos,
            times,
            width=bar_width * 0.92,
            label=method_label(method, width),
            color=method_color(method, width),
            edgecolor="white",
            linewidth=0.4,
        )

        for j, bar in enumerate(bars):
            value = times[j]
            t = threads[j]
            baseline = best_by_thread[t]

            if np.isnan(value) or value <= 0.0 or np.isnan(baseline):
                continue

            ratio = value / baseline
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                value * 1.13,
                f"{ratio:.2f}x",
                ha="center",
                va="bottom",
                rotation=45,
                fontsize=8,
            )

    ax.set_yscale("log")
    ax.set_ylim(ymin, ymax)

    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in threads], fontsize=12)

    ax.set_xlabel("Number of Threads", fontsize=14)
    ax.set_ylabel("Execution Time Per Point (s)", fontsize=14)
    ax.set_title(args.title, fontsize=14)

    ax.grid(True, which="major", axis="y", linestyle=":", linewidth=0.7, alpha=0.7)
    ax.grid(True, which="minor", axis="y", linestyle=":", linewidth=0.4, alpha=0.35)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
        frameon=False,
        fontsize=9,
    )

    fig.subplots_adjust(bottom=0.30)
    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    main()
