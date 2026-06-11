# loblawbio

## To run
#### manually
pip install -r requirements.txt
python load_data.py
streamlit run dashboard.py

##### or...

make setup
make pipeline
make dashboard

You will get the dashboard at http://localhost:8501

## What you get
A streamlit dashboard with sorted csvs and a couple of key figures.
1) Viewing relative immune cell population frequencies by sample.
1) Comparing melanoma PBMC miraclib responders versus non-responders.
1) Visualizing responder/non-responder differences with boxplots.
1) Reporting statistical test results for each immune cell population.
1) Summarizing baseline melanoma PBMC miraclib samples by project, response, and sex.

## Methodology reasoning
Responder and non-responder relative frequencies were compared separately for each immune cell population using a two-sided Mann-Whitney U test. Used this test because the two response groups are independent and the relative frequency values are continuous percentages that may not necessarily follow a normal distribution!

## Scaling up

For hundreds of projects, thousands of samples, and more complex analytics, the schema could be normalized into multiple related tables :)

Currently using a denormalized table, which works for this, but would be good to normalize at scale.