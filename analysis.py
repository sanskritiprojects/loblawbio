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
POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}
POPULATION_COLORS = {
    "b_cell": "#1B4965",
    "cd8_t_cell": "#62B6CB",
    "cd4_t_cell": "#EE6C4D",
    "nk_cell": "#5A189A",
    "monocyte": "#2A9D8F",
}
RESPONSE_COLORS = {
    "yes": "#2A9D8F",
    "no": "#E76F51",
}
RESPONSE_LABELS = {
    "yes": "Responder",
    "no": "Non-responder",
}
SEX_COLORS = {
    "M": "#264653",
    "F": "#E9C46A",
}
PROJECT_COLORS = {
    "prj1": "#1B4965",
    "prj3": "#62B6CB",
}
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


def _label_populations(df: pd.DataFrame, column: str = "population") -> pd.DataFrame:
    labeled = df.copy()
    labeled[column] = labeled[column].map(POPULATION_LABELS).fillna(labeled[column])
    return labeled


def _apply_chart_theme(fig):
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        font={"family": "Arial, sans-serif", "color": "#293241", "size": 13},
        legend={"title": ""},
        margin={"t": 60, "r": 20, "b": 40, "l": 40},
    )
    fig.update_xaxes(showgrid=False, linecolor="#D9D9D9")
    fig.update_yaxes(gridcolor="#EDEDED", linecolor="#D9D9D9")
    return fig

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


def build_mean_frequency_barplot(summary: pd.DataFrame | None = None):
    if summary is None:
        summary = get_summary_table()

    mean_freq = (
        summary.groupby("population", observed=True)["percentage"]
        .mean()
        .reindex(POPULATIONS)
        .reset_index()
    )
    mean_freq = _label_populations(mean_freq)

    fig = px.bar(
        mean_freq,
        x="population",
        y="percentage",
        labels={
            "population": "Cell population",
            "percentage": "Mean relative frequency (%)",
        },
        title="Mean relative frequency by cell population",
    )
    fig.update_traces(
        marker_color=[POPULATION_COLORS[p] for p in POPULATIONS],
        showlegend=False,
    )
    return _apply_chart_theme(fig)


def build_population_distribution_violin(summary: pd.DataFrame | None = None):
    if summary is None:
        summary = get_summary_table()

    plot_data = _label_populations(summary.copy())
    plot_data["population"] = pd.Categorical(
        plot_data["population"],
        categories=[POPULATION_LABELS[p] for p in POPULATIONS],
        ordered=True,
    )

    fig = px.violin(
        plot_data,
        x="population",
        y="percentage",
        color="population",
        box=True,
        points=False,
        labels={
            "population": "Cell population",
            "percentage": "Relative frequency (%)",
        },
        title="Distribution of relative frequencies across samples",
        color_discrete_map={
            POPULATION_LABELS[p]: POPULATION_COLORS[p] for p in POPULATIONS
        },
    )
    fig.update_layout(showlegend=False)
    return _apply_chart_theme(fig)


def build_responder_comparison_boxplot():
    data = get_miraclib_pbmc_melanoma_summary().copy()
    data["population_label"] = pd.Categorical(
        data["population"].map(POPULATION_LABELS),
        categories=[POPULATION_LABELS[p] for p in POPULATIONS],
        ordered=True,
    )
    data["response_label"] = data["response"].map(RESPONSE_LABELS)

    fig = px.box(
        data,
        x="population_label",
        y="percentage",
        color="response_label",
        points="outliers",
        category_orders={
            "population_label": [POPULATION_LABELS[p] for p in POPULATIONS],
            "response_label": list(RESPONSE_LABELS.values()),
        },
        labels={
            "population_label": "Cell population",
            "percentage": "Relative frequency (%)",
            "response_label": "Response",
        },
        title="Relative frequencies by response group",
        color_discrete_map={
            RESPONSE_LABELS["yes"]: RESPONSE_COLORS["yes"],
            RESPONSE_LABELS["no"]: RESPONSE_COLORS["no"],
        },
    )
    return _apply_chart_theme(fig)


def build_responder_median_comparison_barplot(stats: pd.DataFrame | None = None):
    if stats is None:
        stats = compare_responders_vs_nonresponders()

    plot_data = stats[["population", "median_responders", "median_nonresponders"]].copy()
    plot_data = plot_data.melt(
        id_vars="population",
        value_vars=["median_responders", "median_nonresponders"],
        var_name="group",
        value_name="median_percentage",
    )
    plot_data["response_label"] = plot_data["group"].map(
        {
            "median_responders": RESPONSE_LABELS["yes"],
            "median_nonresponders": RESPONSE_LABELS["no"],
        }
    )
    plot_data["population_label"] = pd.Categorical(
        plot_data["population"].map(POPULATION_LABELS),
        categories=[POPULATION_LABELS[p] for p in POPULATIONS],
        ordered=True,
    )

    fig = px.bar(
        plot_data,
        x="population_label",
        y="median_percentage",
        color="response_label",
        barmode="group",
        category_orders={
            "population_label": [POPULATION_LABELS[p] for p in POPULATIONS],
            "response_label": list(RESPONSE_LABELS.values()),
        },
        labels={
            "population_label": "Cell population",
            "median_percentage": "Median relative frequency (%)",
            "response_label": "Response",
        },
        title="Median relative frequency: responders vs non-responders",
        color_discrete_map={
            RESPONSE_LABELS["yes"]: RESPONSE_COLORS["yes"],
            RESPONSE_LABELS["no"]: RESPONSE_COLORS["no"],
        },
    )
    return _apply_chart_theme(fig)


def build_pvalue_barplot(stats: pd.DataFrame | None = None):
    if stats is None:
        stats = compare_responders_vs_nonresponders()

    plot_data = stats.copy()
    plot_data["population_label"] = pd.Categorical(
        plot_data["population"].map(POPULATION_LABELS),
        categories=[POPULATION_LABELS[p] for p in POPULATIONS],
        ordered=True,
    )
    plot_data["significance"] = plot_data["significant_at_0.05"].map(
        {True: "Significant", False: "Not significant"}
    )

    fig = px.bar(
        plot_data,
        x="population_label",
        y="p_value",
        color="significance",
        labels={
            "population_label": "Cell population",
            "p_value": "p-value",
            "significance": "Result",
        },
        title="Statistical significance by cell population",
        color_discrete_map={
            "Significant": RESPONSE_COLORS["yes"],
            "Not significant": "#BDBDBD",
        },
    )
    fig.add_hline(
        y=0.05,
        line_dash="dash",
        line_color="#E76F51",
        annotation_text="p = 0.05",
        annotation_position="top right",
    )
    return _apply_chart_theme(fig)


def build_baseline_project_barplot(samples_by_project: pd.DataFrame):
    fig = px.bar(
        samples_by_project,
        x="project",
        y="sample_count",
        text="sample_count",
        labels={"project": "Project", "sample_count": "Sample count"},
        title="Baseline samples by project",
        color="project",
        color_discrete_map=PROJECT_COLORS,
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False)
    return _apply_chart_theme(fig)


def build_baseline_response_barplot(subjects_by_response: pd.DataFrame):
    plot_data = subjects_by_response.copy()
    plot_data["response_label"] = plot_data["response"].map(RESPONSE_LABELS)

    fig = px.bar(
        plot_data,
        x="response_label",
        y="subject_count",
        text="subject_count",
        labels={"response_label": "Response", "subject_count": "Subject count"},
        title="Baseline subjects by response",
        color="response_label",
        color_discrete_map={
            RESPONSE_LABELS["yes"]: RESPONSE_COLORS["yes"],
            RESPONSE_LABELS["no"]: RESPONSE_COLORS["no"],
        },
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False)
    return _apply_chart_theme(fig)


def build_baseline_sex_barplot(subjects_by_sex: pd.DataFrame):
    fig = px.bar(
        subjects_by_sex,
        x="sex",
        y="subject_count",
        text="subject_count",
        labels={"sex": "Sex", "subject_count": "Subject count"},
        title="Baseline subjects by sex",
        color="sex",
        color_discrete_map=SEX_COLORS,
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False)
    return _apply_chart_theme(fig)


def build_baseline_male_b_cells_by_response(baseline: pd.DataFrame):
    male_baseline = baseline[baseline["sex"] == "M"].copy()
    plot_data = (
        male_baseline.groupby("response", observed=True)["b_cell"]
        .mean()
        .reindex(["yes", "no"])
        .reset_index()
    )
    plot_data["response_label"] = plot_data["response"].map(RESPONSE_LABELS)

    fig = px.bar(
        plot_data,
        x="response_label",
        y="b_cell",
        text=plot_data["b_cell"].map(lambda value: f"{value:.2f}"),
        labels={"response_label": "Response", "b_cell": "Mean B cell count"},
        title="Mean B cell count at baseline (melanoma males)",
        color="response_label",
        color_discrete_map={
            RESPONSE_LABELS["yes"]: RESPONSE_COLORS["yes"],
            RESPONSE_LABELS["no"]: RESPONSE_COLORS["no"],
        },
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False)
    return _apply_chart_theme(fig)


def build_baseline_population_composition(baseline: pd.DataFrame):
    totals = baseline[POPULATIONS].sum()
    composition = (totals / totals.sum() * 100).reset_index()
    composition.columns = ["population", "percentage"]
    composition = _label_populations(composition)

    fig = px.pie(
        composition,
        names="population",
        values="percentage",
        title="Baseline cell population composition",
        color="population",
        color_discrete_map={
            POPULATION_LABELS[p]: POPULATION_COLORS[p] for p in POPULATIONS
        },
    )
    return _apply_chart_theme(fig)


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