import numpy as np
import pandas as pd
from pathlib import Path


BASELINE_DIR = Path("/data1/sandesh/arka/Joint-Heston-Model-Calibration/bates/baseline_real")
OUTPUT_DIR = BASELINE_DIR / "Output"

SINGLE_STATE_ERRORS = OUTPUT_DIR / "single_state_errors.csv"
MULTI_STATE_ERRORS = OUTPUT_DIR / "multi_state_errors.csv"
STOCK_PATH = BASELINE_DIR / "S_path.csv"
REPORT_PATH = OUTPUT_DIR / "Q_errors_report.txt"

ERROR_COLUMNS = [
    "date",
    "hour",
    "strike",
    "maturity",
    "true_price",
    "computed_price",
    "abs_error",
]


def read_bates_errors(path):
    # Bates error rows include hour, but the current CSV header omits it.
    df = pd.read_csv(path, names=ERROR_COLUMNS, skiprows=1)
    df["date"] = pd.to_datetime(df["date"])
    df["hour"] = df["hour"].astype(int)

    numeric_columns = [
        "strike",
        "maturity",
        "true_price",
        "computed_price",
        "abs_error",
    ]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric)
    return df


def get_moneyness_metrics(df, stock_df):
    state_df = df.merge(
        stock_df[["date", "hour", "stock_price"]],
        on=["date", "hour"],
        how="left",
    )

    missing_stock = state_df["stock_price"].isna().sum()
    if missing_stock:
        raise ValueError(f"Missing stock_price for {missing_stock} option error rows")

    state_df["Moneyness_Val"] = state_df["strike"] / state_df["stock_price"]
    conditions = [
        state_df["Moneyness_Val"] < 0.97,
        state_df["Moneyness_Val"].between(0.97, 1.03),
        state_df["Moneyness_Val"] > 1.03,
    ]
    choices = ["ITM", "ATM", "OTM"]
    state_df["Moneyness"] = np.select(conditions, choices, default="Unknown")

    squared_error = (state_df["true_price"] - state_df["computed_price"]) ** 2
    state_df["squared_error"] = squared_error

    moneyness_metrics = (
        state_df.groupby("Moneyness")
        .agg(
            MAE=("abs_error", "mean"),
            **{"Median AE": ("abs_error", "median")},
            RMSE=("squared_error", lambda x: np.sqrt(x.mean())),
            count=("abs_error", "size"),
        )
    )

    overall = {
        "MAE": state_df["abs_error"].mean(),
        "Median AE": state_df["abs_error"].median(),
        "RMSE": np.sqrt(squared_error.mean()),
        "Max AE": state_df["abs_error"].max(),
    }

    time_metrics = (
        state_df.groupby(["date", "hour"])
        .agg(
            MAE=("abs_error", "mean"),
            RMSE=("squared_error", lambda x: np.sqrt(x.mean())),
        )
    )

    time_metrics_summary = {
        "Average Date-Hour MAE": time_metrics["MAE"].mean(),
        "Median Date-Hour MAE": time_metrics["MAE"].median(),
        "Average Date-Hour RMSE": time_metrics["RMSE"].mean(),
    }

    return moneyness_metrics, overall, time_metrics_summary


def format_report_section(title, moneyness_metrics, overall, time_metrics_summary):
    overall_df = pd.Series(overall).to_frame("value")
    time_summary_df = pd.Series(time_metrics_summary).to_frame("value")

    return "\n".join([
        title,
        "=" * len(title),
        "",
        "Overall",
        overall_df.to_string(float_format=lambda x: f"{x:.6f}"),
        "",
        "Moneyness Metrics",
        moneyness_metrics.to_string(float_format=lambda x: f"{x:.6f}"),
        "",
        "Date-Hour Metrics Summary",
        time_summary_df.to_string(float_format=lambda x: f"{x:.6f}"),
        "",
    ])


stock_df = pd.read_csv(STOCK_PATH)
stock_df = stock_df.rename(columns={"Underlying Value": "stock_price"})
stock_df["date"] = pd.to_datetime(stock_df["date"])
stock_df["hour"] = stock_df["hour"].astype(int)

single_state_df = read_bates_errors(SINGLE_STATE_ERRORS)
single_moneyness_metrics, single_overall, single_time_metrics_summary = get_moneyness_metrics(
    single_state_df,
    stock_df,
)

multi_state_df = read_bates_errors(MULTI_STATE_ERRORS)
multi_moneyness_metrics, multi_overall, multi_time_metrics_summary = get_moneyness_metrics(
    multi_state_df,
    stock_df,
)

report = "\n".join([
    format_report_section(
        "Bates Single State Errors",
        single_moneyness_metrics,
        single_overall,
        single_time_metrics_summary,
    ),
    format_report_section(
        "Bates Multi State Errors",
        multi_moneyness_metrics,
        multi_overall,
        multi_time_metrics_summary,
    ),
])

REPORT_PATH.write_text(report)
print(f"Wrote report to {REPORT_PATH}")
