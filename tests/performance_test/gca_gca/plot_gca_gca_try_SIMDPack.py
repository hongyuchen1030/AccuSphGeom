#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_SUMMARY_CSV = Path(
    "tests/performance_test/gca_gca/output/gca_gca_SIMDPack_timing_summary.csv"
)
DEFAULT_REPEATS_CSV = Path(
    "tests/performance_test/gca_gca/output/gca_gca_SIMDPack_timing_repeats.csv"
)
DEFAULT_OLD_CSV = Path(
    "tests/performance_test/gca_gca/output/gca_gca_SIMDPack_timing.csv"
)


def method_label(method: str, width: int) -> str:
    if method == "pure_fp":
        if width == 1:
            return "FP64 GCA-GCA Scalar"
        return f"FP64 GCA-GCA SIMD Vector Width {width}"

    if method == "accux":
        if width == 1:
            return "AccuX GCA-GCA Scalar"
        return f"AccuX GCA-GCA SIMD Vector Width {width}"

    if method == "pure_fp_try":
        if width == 1:
            return "FP64 GCA-GCA Try-API Scalar"
        return f"FP64 GCA-GCA Try-API SIMD Vector Width {width}"

    if method == "accux_try":
        if width == 1:
            return "AccuX GCA-GCA Try-API Scalar"
        return f"AccuX GCA-GCA Try-API SIMD Vector Width {width}"

    return f"{method} width {width}"


def method_color(method: str, width: int):
    colors = {
        ("accux", 1): "#2ca02c",
        ("accux", 2): "#4daf4a",
        ("accux", 4): "#8bc34a",
        ("accux", 8): "#b2df8a",
        ("accux_try", 1): "#2ca02c",
        ("accux_try", 2): "#4daf4a",
        ("accux_try", 4): "#8bc34a",
        ("accux_try", 8): "#b2df8a",
        ("pure_fp", 1): "#1f77b4",
        ("pure_fp", 2): "#377eb8",
        ("pure_fp", 4): "#6baed6",
        ("pure_fp", 8): "#9ecae1",
        ("pure_fp_try", 1): "#1f77b4",
        ("pure_fp_try", 2): "#377eb8",
        ("pure_fp_try", 4): "#6baed6",
        ("pure_fp_try", 8): "#9ecae1",
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

    df = df.copy()

    if "median_time" in df.columns:
        df["time"] = df["median_time"].astype(float)
        return df[["method", "threadsNum", "vec_width", "time"]]

    if "time" in df.columns:
        df["time"] = df["time"].astype(float)

        if "repeat" in df.columns:
            grouped = (
                df.groupby(["method", "threadsNum", "vec_width"], as_index=False)[
                    "time"
                ]
                .median()
                .copy()
            )
            return grouped

        return df[["method", "threadsNum", "vec_width", "time"]]

    raise ValueError("CSV must contain either 'median_time' or 'time'.")


def output_with_suffix(output_path: Path, suffix: str) -> Path:
    return output_path.with_name(f"{output_path.stem}{suffix}{output_path.suffix}")


def plot_method_block(
    df: pd.DataFrame,
    output_path: Path,
    title: str,
    ordered_series: list[tuple[str, int]],
    dpi: int,
) -> bool:
    available_series = []
    for method, width in ordered_series:
        sub = df[(df["method"] == method) & (df["vec_width"] == width)]
        if not sub.empty:
            available_series.append((method, width))

    if not available_series:
        return False

    block_df = df[
        df.apply(
            lambda row: (row["method"], int(row["vec_width"])) in available_series,
            axis=1,
        )
    ].copy()

    threads = sorted(block_df["threadsNum"].unique())

    pivot = {}
    for method, width in available_series:
        sub = block_df[
            (block_df["method"] == method) & (block_df["vec_width"] == width)
        ]
        value_by_thread = dict(zip(sub["threadsNum"], sub["time"]))
        pivot[(method, width)] = [value_by_thread.get(t, np.nan) for t in threads]

    positive_times = block_df[block_df["time"] > 0.0]["time"].to_numpy()
    if positive_times.size == 0:
        raise ValueError("No positive timing values found in the selected series.")

    global_best = positive_times.min()

    fig, ax = plt.subplots(figsize=(12, 5.8))

    x = np.arange(len(threads))
    group_width = 0.82
    bar_width = group_width / len(available_series)
    start_offset = -0.5 * group_width + 0.5 * bar_width

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

            if np.isnan(value) or value <= 0.0:
                continue

            ratio = value / global_best
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
    ax.set_title(title, fontsize=14)

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
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote plot to {output_path}")
    return True


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
        help="Output figure path for the try_gca_gca_intersection plot.",
    )
    parser.add_argument(
        default="GCA-GCA Try-API Performance",
        dest="title",
        help="Plot title for try_gca_gca_intersection.",
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

    try_series = [
        ("accux_try", 1),
        ("accux_try", 2),
        ("accux_try", 4),
        ("accux_try", 8),
        ("pure_fp_try", 1),
        ("pure_fp_try", 2),
        ("pure_fp_try", 4),
        ("pure_fp_try", 8),
    ]

    wrote_try = plot_method_block(
        df=df,
        output_path=output_path,
        title=args.title,
        ordered_series=try_series,
        dpi=args.dpi,
    )

    if not wrote_try:
        raise ValueError(
            "No known rows found. Expected pure_fp_try/accux_try."
        )


if __name__ == "__main__":
    main()
