# loblawbio

## To run

pip install -r requirements.txt
python load_data.py
streamlit run dashboard.py

## What you get
A streamlit dashboard with sorted csvs and a couple of key figures.
1) Viewing relative immune cell population frequencies by sample.
1) Comparing melanoma PBMC miraclib responders versus non-responders.
1) Visualizing responder/non-responder differences with boxplots.
1) Reporting statistical test results for each immune cell population.
1) Summarizing baseline melanoma PBMC miraclib samples by project, response, and sex.

## Methodology reasoning
Responder and non-responder relative frequencies were compared separately for each immune cell population using a two-sided Mann-Whitney U test. Used this test because the two response groups are independent and the relative frequency values are continuous percentages that may not necessarily follow a normal distribution!