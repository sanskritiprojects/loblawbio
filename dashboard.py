######################################################
### imports
######################################################

import streamlit as st

from analysis import (
    SUMMARY_OUTPUT_COLUMNS,
    build_baseline_male_b_cells_by_response,
    build_baseline_population_composition,
    build_baseline_project_barplot,
    build_baseline_response_barplot,
    build_baseline_sex_barplot,
    build_mean_frequency_barplot,
    build_population_distribution_violin,
    build_pvalue_barplot,
    build_responder_comparison_boxplot,
    build_responder_median_comparison_barplot,
    compare_responders_vs_nonresponders,
    get_average_b_cells_melanoma_male_responders_baseline,
    get_baseline_miraclib_melanoma_pbmc_samples,
    get_baseline_subset_counts,
    get_summary_table,
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

    st.dataframe(summary[SUMMARY_OUTPUT_COLUMNS], use_container_width=True)

    st.download_button(
        label="Download summary table as CSV",
        data=summary[SUMMARY_OUTPUT_COLUMNS].to_csv(index=False),
        file_name="cell_population_summary.csv",
        mime="text/csv",
    )

    st.subheader("Summary figures")
    col1, col2 = st.columns(2)

    with col1:
        st.plotly_chart(
            build_mean_frequency_barplot(summary),
            use_container_width=True,
        )

    with col2:
        st.plotly_chart(
            build_population_distribution_violin(summary),
            use_container_width=True,
        )


with tab2:
    st.header("Melanoma PBMC samples treated with miraclib")

    stats = compare_responders_vs_nonresponders()

    st.plotly_chart(
        build_responder_comparison_boxplot(),
        use_container_width=True,
    )

    col1, col2 = st.columns(2)

    with col1:
        st.plotly_chart(
            build_responder_median_comparison_barplot(stats),
            use_container_width=True,
        )

    with col2:
        st.plotly_chart(
            build_pvalue_barplot(stats),
            use_container_width=True,
        )

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

    st.subheader("Baseline summary figures")
    fig_col1, fig_col2 = st.columns([1, 1])

    with fig_col1:
        st.plotly_chart(
            build_baseline_population_composition(baseline),
            use_container_width=True,
        )

    with fig_col2:
        st.plotly_chart(
            build_baseline_male_b_cells_by_response(baseline),
            use_container_width=True,
        )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Samples by project")
        st.plotly_chart(
            build_baseline_project_barplot(samples_by_project),
            use_container_width=True,
        )
        st.dataframe(samples_by_project, use_container_width=True)

    with col2:
        st.subheader("Subjects by response")
        st.plotly_chart(
            build_baseline_response_barplot(subjects_by_response),
            use_container_width=True,
        )
        st.dataframe(subjects_by_response, use_container_width=True)

    with col3:
        st.subheader("Subjects by sex")
        st.plotly_chart(
            build_baseline_sex_barplot(subjects_by_sex),
            use_container_width=True,
        )
        st.dataframe(subjects_by_sex, use_container_width=True)

    avg_b_cells = get_average_b_cells_melanoma_male_responders_baseline()
    st.subheader("Melanoma males: average B cells at baseline (responders)")
    st.metric("Average B cell count", f"{avg_b_cells:.2f}")
