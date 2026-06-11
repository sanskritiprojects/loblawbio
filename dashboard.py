######################################################
### imports
######################################################

import plotly.express as px
import streamlit as st

from analysis import (
    get_summary_table,
    get_miraclib_pbmc_melanoma_summary,
    compare_responders_vs_nonresponders,
    get_baseline_miraclib_melanoma_pbmc_samples,
    get_baseline_subset_counts,
)


st.set_page_config(
    page_title="Immune Cell Trial Dashboard",
    layout="wide",
)

st.title("Immune Cell Population Analysis")
st.write("Clinical trial dashboard for Loblaw Bio")


tab1, tab2, tab3 = st.tabs(
    [
        "Part 2: Data Overview",
        "Part 3: Responder Analysis",
        "Part 4: Baseline Subset",
    ]
)


with tab1:
    st.header("Relative frequency of each cell type in each sample")

    summary = get_summary_table()

    display_columns = [
        "sample",
        "total_count",
        "population",
        "count",
        "percentage",
    ]

    st.dataframe(summary[display_columns], use_container_width=True)

    st.download_button(
        label="Download summary table as CSV",
        data=summary[display_columns].to_csv(index=False),
        file_name="cell_population_summary.csv",
        mime="text/csv",
    )


with tab2:
    st.header("Melanoma PBMC samples treated with miraclib")

    data = get_miraclib_pbmc_melanoma_summary()
    stats = compare_responders_vs_nonresponders()

    fig = px.box(
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

    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Statistical test results")
    st.dataframe(stats, use_container_width=True)

    significant = stats[stats["significant_at_0.05"]]

    if len(significant) == 0:
        st.info("No cell populations are significant at p < 0.05.")
    else:
        st.success("Significant cell populations at p < 0.05:")
        st.dataframe(significant, use_container_width=True)


with tab3:
    st.header("Baseline melanoma PBMC samples treated with miraclib")

    baseline = get_baseline_miraclib_melanoma_pbmc_samples()
    samples_by_project, subjects_by_response, subjects_by_sex = get_baseline_subset_counts()

    st.subheader("Matching baseline samples")
    st.dataframe(baseline, use_container_width=True)

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Samples by project")
        st.dataframe(samples_by_project, use_container_width=True)

    with col2:
        st.subheader("Subjects by response")
        st.dataframe(subjects_by_response, use_container_width=True)

    with col3:
        st.subheader("Subjects by sex")
        st.dataframe(subjects_by_sex, use_container_width=True)