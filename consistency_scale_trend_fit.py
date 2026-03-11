from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def fit_linear(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    A = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ beta
    sse = float(((y - yhat) ** 2).sum())
    return beta, yhat, sse


def fit_log(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    lx = np.log(x)
    A = np.column_stack([np.ones_like(lx), lx])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ beta
    sse = float(((y - yhat) ** 2).sum())
    return beta, yhat, sse


def fit_inverse(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    invx = 1.0 / x
    A = np.column_stack([np.ones_like(invx), invx])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ beta
    sse = float(((y - yhat) ** 2).sum())
    return beta, yhat, sse


def fit_power(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float] | None:
    if np.any(y <= 0):
        return None
    lx = np.log(x)
    ly = np.log(y)
    A = np.column_stack([np.ones_like(lx), lx])
    beta, *_ = np.linalg.lstsq(A, ly, rcond=None)
    lyhat = A @ beta
    yhat = np.exp(lyhat)
    sse = float(((y - yhat) ** 2).sum())
    return beta, yhat, sse


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fit trend for minimal consistency scale c_N.")
    parser.add_argument(
        "--config",
        default="config_width_height_consistency_scale_scan.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    report_path = out_dir / "width_height_consistency_scale_scan_report.csv"

    report_df = pd.read_csv(report_path)
    report_df["scale_factor"] = report_df["scale_label"].astype(float)

    threshold_rows = []
    for n, sub in report_df.groupby("n", sort=True):
        cross = sub[sub["status"] == "crossing"].sort_values("scale_factor")
        first = cross.iloc[0] if not cross.empty else None
        threshold_rows.append(
            {
                "n": int(n),
                "c_n_min": float(first["scale_factor"]) if first is not None else math.nan,
                "gamma_c_at_threshold": float(first["gamma_c_est"]) if first is not None else math.nan,
            }
        )
    threshold_df = pd.DataFrame(threshold_rows)
    fit_df = threshold_df.dropna().copy()

    x = fit_df["n"].to_numpy(dtype=float)
    y = fit_df["c_n_min"].to_numpy(dtype=float)

    fit_rows = []
    beta_lin, yhat_lin, sse_lin = fit_linear(x, y)
    fit_rows.append({"model": "linear", "sse": sse_lin, "param_0": float(beta_lin[0]), "param_1": float(beta_lin[1])})

    beta_log, yhat_log, sse_log = fit_log(x, y)
    fit_rows.append({"model": "log", "sse": sse_log, "param_0": float(beta_log[0]), "param_1": float(beta_log[1])})

    beta_inv, yhat_inv, sse_inv = fit_inverse(x, y)
    fit_rows.append({"model": "plateau_inverse_n", "sse": sse_inv, "param_0": float(beta_inv[0]), "param_1": float(beta_inv[1])})

    power_fit = fit_power(x, y)
    if power_fit is not None:
        beta_pow, yhat_pow, sse_pow = power_fit
        fit_rows.append({"model": "power", "sse": sse_pow, "param_0": float(beta_pow[0]), "param_1": float(beta_pow[1])})
    else:
        yhat_pow = None

    fit_summary = pd.DataFrame(fit_rows).sort_values("sse")
    fit_summary.to_csv(out_dir / "consistency_scale_trend_fit_summary.csv", index=False, encoding="utf-8-sig")
    threshold_df.to_csv(out_dir / "consistency_scale_thresholds.csv", index=False, encoding="utf-8-sig")

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.scatter(x, y, color="#1f77b4", s=55, label="Observed c_N")

    x_plot = np.linspace(float(x.min()), float(x.max()), 300)
    for model, color in [("linear", "#ff7f0e"), ("log", "#2ca02c"), ("plateau_inverse_n", "#d62728")]:
        row = fit_summary[fit_summary["model"] == model].iloc[0]
        if model == "linear":
            y_plot = row["param_0"] + row["param_1"] * x_plot
        elif model == "log":
            y_plot = row["param_0"] + row["param_1"] * np.log(x_plot)
        else:
            y_plot = row["param_0"] + row["param_1"] / x_plot
        ax.plot(x_plot, y_plot, color=color, linewidth=1.8, label=f"{model} fit")

    ax.set_xlabel("N")
    ax.set_ylabel("Minimal restoring scale factor c_N")
    best = fit_summary.iloc[0]
    ax.set_title(f"Minimal consistency scale c_N vs N\nBest simple fit: {best['model']} (SSE={best['sse']:.4f})")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "consistency_scale_trend_fit.png", dpi=160, bbox_inches="tight")

    print((out_dir / "consistency_scale_thresholds.csv").as_posix())
    print((out_dir / "consistency_scale_trend_fit_summary.csv").as_posix())
    print((out_dir / "consistency_scale_trend_fit.png").as_posix())
    print()
    print(threshold_df.to_string(index=False))
    print()
    print(fit_summary.to_string(index=False))
