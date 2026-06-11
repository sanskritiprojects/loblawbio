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

