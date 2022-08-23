import streamlit as st
import numpy as np
import plotly.express as px
import pandas as pd


def show_expression(countData, sampleData, annot_dict={},
                    sample_col="sampleID", default_annot='FLYBASE', add_annot='SYMBOL'):
    df = countData.copy()
    df.index.name = 'Gene'
    df = df.reset_index()
    if annot_dict:
        gene_name = st.radio('Choose gene annotation', [default_annot, add_annot])
    else:
        gene_name = default_annot
    if gene_name == default_annot:
        gene_options = df['Gene'].unique()
    else:
        gene_options = list(annot_dict.keys())
        #gene_options = [a for a,b in annot_dict.items() if b in df.Gene.values]
    genes = st.multiselect("Choose gene(s) of interest", gene_options, key='gois')
    if gene_name == add_annot:
        genes = set([annot_dict[g] for g in genes])
    st.write(genes)
    df = df.apply(lambda x: np.log2(x + 0.5) if np.issubdtype(x.dtype, np.number) else x)
    sampleDataAb = sampleData.reset_index()
    df = df.reset_index()
    c1, c2 = st.columns(2)
    compare_by = c1.selectbox('Compare by', sampleDataAb.columns)
    categories = c1.multiselect(f'Categories of {compare_by} to display',
                                ['All'] + list(sampleDataAb[compare_by].unique()), key='cats')
    filter_by = c2.selectbox("Filter by",  sampleDataAb.columns)
    filter_out = c2.selectbox(f'Which category of {filter_by} to keep?',
                               [None] + list(sampleDataAb[filter_by].unique()))

    if 'All' in categories:
        categories = list(sampleDataAb[compare_by].unique())
    if genes:
        if len(genes) * len(categories) > 40:
            st.write('Too many genes/categories to display, consider choosing fewer genes')
            st.stop()
        c3, c4 = st.columns(2)
        tpm_label = 'log2 (TPM)'
        gene_df = df[df['Gene'].isin(genes)]
        sample_df = sampleDataAb[sampleDataAb[compare_by].isin(categories)]

        if filter_out:
            sample_df = sample_df[sample_df[filter_by] == filter_out]
        gene_df2 = (gene_df.melt(id_vars=['Gene'], value_name=tpm_label, var_name=sample_col)
                    .merge(sample_df, how='inner', on=sample_col))
        groupby = st.radio('Group by', ['Gene', compare_by])
        color_by = [c for c in ['Gene', compare_by] if c != groupby][0]
        fig = px.box(gene_df2, x=groupby, y=tpm_label, color=color_by,
                     hover_data=['Gene'] + list(sampleData.columns), points='all')
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=16))
        fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
        st.plotly_chart(fig, use_container_width=True)
    return genes


def show_expression_heatmap():
    pass

