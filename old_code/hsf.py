import streamlit as st
from processing.load_data import load_data, load_from_bucket, load_annotations
from graphs.show_pca import show_pca
from graphs.show_gene_expression import show_expression
from graphs.show_dge import show_volcano, link_to_string, download_filtered_hits
from check_password import check_password
import plotly.express as px
page_name ='HSF RNAseq'
st.set_page_config(layout="wide")

def app():
    st.title('HSF')
    results, tpms, vsd, sd = load_from_bucket('jagannathan', experiment_name='hsf', gene_name='')
    #results, tpms, vsd, sd = load_data('/Users/ansintsova/git_repos/fly_rnaseq/data/hsf/results', '')
    st.header("PCA")
    with st.expander('Show PCA'):
        show_pca(vsd, sd)

    st.header("Expression")
    with st.expander('Show Gene Expression'):
        annotations = load_annotations()['FLYBASE']
        annotations = {v:k for k,v in annotations.items()}
        tpms = tpms.reset_index()
        tpms = tpms.melt(id_vars='index', var_name='sampleID', value_name='Normalised Counts')
        tpms['SYMBOL'] = tpms['index'].map(annotations)
        tpms = tpms.rename({'index': 'FLYBASE'}, axis=1)
        tpms = tpms.merge(sd.reset_index(), left_on='sampleID', right_on='sampleID')
        gene_name = st.radio('Choose gene annotation', ['FLYBASE', 'SYMBOL'])
        # if 'fly_genes' not in st.session_state.keys():
        #     st.session_state['fly_genes'] = list(tpms.FLYBASE.unique())
        #
        # if 'sym_genes' not in st.session_state.keys():
        #     st.session_state['sym_genes'] = list(tpms.SYMBOL.unique())
        if gene_name == 'FLYBASE':
            genes = list(tpms.FLYBASE.unique())
            GENE_COL = 'FLYBASE'
        else:
            genes = list(tpms.SYMBOL.unique())
            GENE_COL = 'SYMBOL'

        genes_to_show = st.multiselect('Choose gene of interest', genes)
        c1, c2 = st.columns(2)
        compare_by = c1.selectbox('Compare by', sd.columns)
        categories = c1.multiselect(f'Categories of {compare_by} to display',
                                    ['All'] + list(tpms[compare_by].unique()), key='cats')
        filter_by = c2.selectbox("Filter by", sd.columns)
        filter_out = c2.selectbox(f'Which category of {filter_by} to keep?',
                                  [None] + list(tpms[filter_by].unique()))

        if 'All' in categories:
            categories = list(tpms[compare_by].unique())
        if genes_to_show:
            if len(genes_to_show) * len(categories) > 40:
                st.write('Too many genes/categories to display, consider choosing fewer genes')
                st.stop()
            c3, c4 = st.columns(2)
            tpm_label = 'Normalised Counts'
            gene_df = tpms[(tpms[GENE_COL].isin(genes_to_show)) &(tpms[compare_by].isin(categories))]
            if filter_out:
                gene_df = gene_df[gene_df[filter_by] == filter_out]
            groupby = st.radio('Group by', [GENE_COL, compare_by])
            color_by = [c for c in [GENE_COL, compare_by] if c != groupby][0]
            fig = px.box(gene_df, x=groupby, y=tpm_label, color=color_by, log_y=True,
                         hover_data=[GENE_COL] + list(sd.columns), points='all')
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=16))
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
            st.plotly_chart(fig, use_container_width=True)


    st.header("DGE Results")
    with st.expander('Show DGE Results'):

        hits_df, fig = show_volcano(results, gene_name='Gene', genes_to_highlight=list(genes_to_show),
                               annotation_col='SYMBOL')
        c1, c2 = st.columns([3,1])
        c1.plotly_chart(fig, use_container_width=True)
        c2.markdown("### STRING PPI Network")
        link_to_string(hits_df, c2)
        c2.markdown("### Download results")
        download_filtered_hits(hits_df, c2)


if check_password():
    app()