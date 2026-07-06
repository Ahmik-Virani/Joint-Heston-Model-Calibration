import numpy as np
import pandas as pd
from pathlib import Path

single_state_errors = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output/single_state_errors.csv"
multi_state_errors = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output/multi_state_errors.csv"
stock_path = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/S_path.csv"
output_dir = Path("/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output")
report_path = output_dir / "Q_errors_report.txt"

def getMoneyness(df, stock_df):
    state_df = df.merge(stock_df[["date", "stock_price"]], on="date", how="left")
    state_df["Moneyness_Val"] = state_df["strike"] / state_df["stock_price"]
    conditions = [
        state_df["Moneyness_Val"] < 0.97,
        state_df["Moneyness_Val"].between(0.97, 1.03),
        state_df["Moneyness_Val"] > 1.03,
    ]
    choices = ["ITM", "ATM", "OTM"]
    state_df["Moneyness"] = np.select(
        conditions,
        choices,
        default="Unknown",
    )

    moneyness_metrics = (
        state_df.groupby("Moneyness")
        .apply(lambda x: pd.Series({
            "MAE": x["abs_error"].mean(),
            "Median AE": x["abs_error"].median(),
            "RMSE": np.sqrt(((x["true_price"] - x["computed_price"]) ** 2).mean()),
            "count": len(x),
        }))
    )

    overall = {
        "MAE": state_df["abs_error"].mean(),
        "Median AE": state_df["abs_error"].median(),
        "RMSE": np.sqrt(((state_df["true_price"] - state_df["computed_price"]) ** 2).mean()),
        "Max AE": state_df["abs_error"].max(),
    }
    
    daily_metrics = (
        state_df
        .groupby("date")
        .apply(lambda x: pd.Series({
            "MAE": x["abs_error"].mean(),
            "RMSE": np.sqrt(((x["true_price"] - x["computed_price"]) ** 2).mean()),
        }))
    )
    
    daily_metrics_summary = {
        "Average Daily MAE": daily_metrics["MAE"].mean(),
        "Median Daily MAE": daily_metrics["MAE"].median(),
        "Average Daily RMSE": daily_metrics["RMSE"].mean(),
    }

    return moneyness_metrics, overall, daily_metrics_summary


def format_report_section(title, moneyness_metrics, overall, daily_metrics_summary):
    overall_df = pd.Series(overall).to_frame("value")
    daily_summary_df = pd.Series(daily_metrics_summary).to_frame("value")

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
        "Daily Metrics Summary",
        daily_summary_df.to_string(float_format=lambda x: f"{x:.6f}"),
        "",
    ])


stock_df = pd.read_csv(stock_path)
stock_df = stock_df.rename(columns={"Date": "date", "Underlying Value": "stock_price"})
stock_df["date"] = pd.to_datetime(stock_df["date"])

df = pd.read_csv(single_state_errors)
df["date"] = pd.to_datetime(df["date"])
single_moneyness_metrics, single_overall, single_daily_metrics_summary = getMoneyness(df, stock_df)

df = pd.read_csv(multi_state_errors)
df["date"] = pd.to_datetime(df["date"])
multi_moneyness_metrics, multi_overall, multi_daily_metrics_summary = getMoneyness(df, stock_df)

report = "\n".join([
    format_report_section(
        "Single State Errors",
        single_moneyness_metrics,
        single_overall,
        single_daily_metrics_summary,
    ),
    format_report_section(
        "Multi State Errors",
        multi_moneyness_metrics,
        multi_overall,
        multi_daily_metrics_summary,
    ),
])

report_path.write_text(report)
print(f"Wrote report to {report_path}")
