######################################################
### imports
######################################################

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
from scipy.stats import mannwhitneyu

######################################################
### variables
######################################################

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
SUMMARY_OUTPUT_COLUMNS = [
    "sample",
    "total_count",
    "population",
    "count",
    "percentage",
]
OUTPUT_DIR = ROOT / "outputs"

######################################################
### functions
######################################################

def get_connection():
    return sqlite3.connect(DB_PATH)

def get_summary_table():
    query = """
    SELECT
        project,
        subject,
        condition,
        age,
        sex,
        treatment,
        response,
        sample,
        sample_type,
        time_from_treatment_start,
        b_cell,
        cd8_t_cell,
        cd4_t_cell,
        nk_cell,
        monocyte,
        b_cell + cd8_t_cell + cd4_t_cell + nk_cell + monocyte AS total_count
    FROM cell_counts;
    """

    with get_connection() as conn:
        wide = pd.read_sql_query(query, conn)

    summary = wide.melt(
        id_vars=[
            "project",
            "subject",
            "condition",
            "age",
            "sex",
            "treatment",
            "response",
            "sample",
            "sample_type",
            "time_from_treatment_start",
            "total_count",
        ],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )

    summary["percentage"] = summary["count"] / summary["total_count"] * 100

    return summary[
        [
            "sample",
            "total_count",
            "population",
            "count",
            "percentage",
            "project",
            "subject",
            "condition",
            "age",
            "sex",
            "treatment",
            "response",
            "sample_type",
            "time_from_treatment_start",
        ]
    ]

### melanoma summary
def get_miraclib_pbmc_melanoma_summary():
    summary = get_summary_table()

    return summary[
        (summary["condition"] == "melanoma")
        & (summary["treatment"] == "miraclib")
        & (summary["sample_type"] == "PBMC")
        & (summary["response"].isin(["yes", "no"]))
    ].copy()

### responder comparison
def compare_responders_vs_nonresponders():
    data = get_miraclib_pbmc_melanoma_summary()

    results = []

    for population in POPULATIONS:
        responders = data[
            (data["population"] == population) & (data["response"] == "yes")
        ]["percentage"]

        nonresponders = data[
            (data["population"] == population) & (data["response"] == "no")
        ]["percentage"]

        test = mannwhitneyu(
            responders,
            nonresponders,
            alternative="two-sided",
            method="auto",
        )

        results.append(
            {
                "population": population,
                "n_responders": len(responders),
                "n_nonresponders": len(nonresponders),
                "median_responders": responders.median(),
                "median_nonresponders": nonresponders.median(),
                "u_statistic": test.statistic,
                "p_value": test.pvalue,
                "significant_at_0.05": test.pvalue < 0.05,
            }
        )

    return pd.DataFrame(results)


def build_responder_comparison_boxplot():
    data = get_miraclib_pbmc_melanoma_summary()

    return px.box(
        data,
        x="population",
        y="percentage",
        color="response",
        points="outliers",
        labels={
            "population": "Cell population",
            "percentage": "Relative frequency (%)",
            "response": "Response",
        },
        title="Relative frequencies by response group",
    )


## baselines
def get_baseline_miraclib_melanoma_pbmc_samples():
    query = """
    SELECT *
    FROM cell_counts
    WHERE condition = 'melanoma'
      AND treatment = 'miraclib'
      AND sample_type = 'PBMC'
      AND time_from_treatment_start = 0;
    """

    with get_connection() as conn:
        return pd.read_sql_query(query, conn)


def get_baseline_subset_counts():
    baseline = get_baseline_miraclib_melanoma_pbmc_samples()

    samples_by_project = (
        baseline.groupby("project")["sample"]
        .nunique()
        .reset_index(name="sample_count")
    )

    subjects_by_response = (
        baseline.drop_duplicates("subject")
        .groupby("response")["subject"]
        .nunique()
        .reset_index(name="subject_count")
    )

    subjects_by_sex = (
        baseline.drop_duplicates("subject")
        .groupby("sex")["subject"]
        .nunique()
        .reset_index(name="subject_count")
    )

    return samples_by_project, subjects_by_response, subjects_by_sex


def get_average_b_cells_melanoma_male_responders_baseline() -> float:
    """Average B cell count for melanoma male responders at baseline (time=0)."""
    baseline = get_baseline_miraclib_melanoma_pbmc_samples()
    subset = baseline[(baseline["sex"] == "M") & (baseline["response"] == "yes")]

    if subset.empty:
        return 0.0

    return round(subset["b_cell"].mean(), 2)


def write_part2_outputs(output_dir: Path = OUTPUT_DIR) -> Path:
    summary = get_summary_table()
    output_path = output_dir / "cell_population_summary.csv"
    summary[SUMMARY_OUTPUT_COLUMNS].to_csv(output_path, index=False)
    print(f"Wrote Part 2 summary: {output_path}")
    return output_path


def write_part3_outputs(output_dir: Path = OUTPUT_DIR) -> None:
    stats = compare_responders_vs_nonresponders()
    stats_path = output_dir / "responder_comparison_stats.csv"
    stats.to_csv(stats_path, index=False)
    print(f"Wrote Part 3 stats: {stats_path}")

    plot_path = output_dir / "responder_comparison_boxplot.html"
    build_responder_comparison_boxplot().write_html(plot_path)
    print(f"Wrote Part 3 boxplot: {plot_path}")

    significant = stats[stats["significant_at_0.05"]]
    significant_path = output_dir / "significant_populations.csv"
    significant.to_csv(significant_path, index=False)
    print(f"Wrote Part 3 significant populations: {significant_path}")

    if significant.empty:
        print("No cell populations are significant at p < 0.05.")
    else:
        populations = ", ".join(significant["population"].tolist())
        print(f"Significant cell populations at p < 0.05: {populations}")


def write_part4_outputs(output_dir: Path = OUTPUT_DIR) -> None:
    samples_by_project, subjects_by_response, subjects_by_sex = get_baseline_subset_counts()
    avg_b_cells = get_average_b_cells_melanoma_male_responders_baseline()

    samples_by_project.to_csv(output_dir / "baseline_samples_by_project.csv", index=False)
    subjects_by_response.to_csv(output_dir / "baseline_subjects_by_response.csv", index=False)
    subjects_by_sex.to_csv(output_dir / "baseline_subjects_by_sex.csv", index=False)

    avg_path = output_dir / "average_b_cells_melanoma_male_responders_baseline.csv"
    pd.DataFrame(
        [{"average_b_cells_melanoma_male_responders_baseline": avg_b_cells}]
    ).to_csv(avg_path, index=False)

    print(f"Wrote Part 4 outputs to: {output_dir}")
    print(f"Average B cells (melanoma, male, responder, baseline): {avg_b_cells:.2f}")


def run_pipeline(output_dir: Path = OUTPUT_DIR) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    write_part2_outputs(output_dir)
    write_part3_outputs(output_dir)
    write_part4_outputs(output_dir)

    print(f"Pipeline complete. Outputs written to: {output_dir}")