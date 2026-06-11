######################################################
### imports
######################################################

import sqlite3
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu

######################################################
### variables
######################################################

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

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
