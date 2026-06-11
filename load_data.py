######################################################
### imports
######################################################

import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
TABLE_NAME = "cell_counts"

REQUIRED_COLUMNS = [
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
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


SCHEMA_SQL = f"""
DROP TABLE IF EXISTS {TABLE_NAME};

CREATE TABLE {TABLE_NAME} (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    project TEXT NOT NULL,
    subject TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL,
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT,
    sample TEXT NOT NULL UNIQUE,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    b_cell INTEGER NOT NULL,
    cd8_t_cell INTEGER NOT NULL,
    cd4_t_cell INTEGER NOT NULL,
    nk_cell INTEGER NOT NULL,
    monocyte INTEGER NOT NULL
);

CREATE INDEX idx_cell_counts_sample
ON {TABLE_NAME}(sample);

CREATE INDEX idx_cell_counts_subject
ON {TABLE_NAME}(subject);

CREATE INDEX idx_cell_counts_analysis_filters
ON {TABLE_NAME}(
    condition,
    treatment,
    sample_type,
    time_from_treatment_start,
    response
);
"""


######################################################
### main functions
######################################################
# Create sqllite database connection
def get_database_connection(db_path: Path) -> sqlite3.Connection:
    return sqlite3.connect(db_path)

# create the database schema
def initialize_database(connection: sqlite3.Connection) -> None:
    connection.executescript(SCHEMA_SQL)
    connection.commit()

# read cell-count csv file into a pandas DataFrame
def read_csv(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"Could not find CSV file: {csv_path}")

    return pd.read_csv(csv_path)

######################################################
# main workflow
######################################################
def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only required columns and make sure numeric columns
    have the expected integer type.
    """
    df = df[REQUIRED_COLUMNS].copy()

    integer_columns = [
        "age",
        "time_from_treatment_start",
        "b_cell",
        "cd8_t_cell",
        "cd4_t_cell",
        "nk_cell",
        "monocyte",
    ]

    for column in integer_columns:
        df[column] = df[column].astype(int)

    text_columns = [
        "project",
        "subject",
        "condition",
        "sex",
        "treatment",
        "response",
        "sample",
        "sample_type",
    ]

    for column in text_columns:
        df[column] = df[column].astype("string")

    return df

# load df 
def load_dataframe(
    connection: sqlite3.Connection,
    df: pd.DataFrame,
    table_name: str = TABLE_NAME,) -> None:
    """Load the cleaned DataFrame into the SQLite table."""
    df.to_sql(
        table_name,
        connection,
        if_exists="append",
        index=False,
        )
    connection.commit()

# main function to run the workflow
def load_data(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> None:
    raw_df = read_csv(csv_path)
    clean_df = clean_data(raw_df)

    with get_database_connection(db_path) as connection:
        initialize_database(connection)
        load_dataframe(connection, clean_df)

    print(f"Created database: {db_path}")
    print(f"Loaded {len(clean_df)} rows into table: {TABLE_NAME}")

def main() -> None:
    load_data()


if __name__ == "__main__":
    main()