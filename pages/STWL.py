import streamlit as st
from processing.load_data import load_data, load_from_bucket, load_annotations
from graphs.dge_dataset import DgeDataSet
from graphs.show_pca import show_pca
from graphs.show_gene_expression import show_expression
from graphs.show_dge import show_volcano, link_to_string, download_filtered_hits
from check_password import check_password
import plotly.express as px

page_name ='STWL OE'
st.set_page_config(layout="wide")


def app():
    st.title('STWL OE')
    dge = DgeDataSet(bucket='jagannathan', experiment_name='stwloe')
    st.header("PCA")
    with st.expander('Show PCA'):
        show_pca(dge.vsd, dge.sd)
    st.header("Expression")
    with st.expander('Show Gene Expression'):
        dge.load_annotations()
        gene_name = st.radio('Choose gene annotation', ['FLYBASE', 'SYMBOL'])
        if gene_name == 'FLYBASE':
            genes = list(dge.tpms.FLYBASE.unique())
            GENE_COL = 'FLYBASE'
        else:
            genes = list(dge.tpms.SYMBOL.unique())
            GENE_COL = 'SYMBOL'

        genes_to_show = st.multiselect('Choose gene of interest', genes)
        c1, c2 = st.columns(2)
        compare_by = c1.selectbox('Compare by', dge.sd.columns)
        categories = c1.multiselect(f'Categories of {compare_by} to display',
                                    ['All'] + list(dge.tpms[compare_by].unique()), default= 'All', key='cats')
        filter_by = c2.selectbox("Filter by", [None] + list(dge.sd.columns))

        filter_out = c2.selectbox(f'Which category of {filter_by} to keep?',
                                       list(dge.tpms[filter_by].unique())) if filter_by else None
        if 'All' in categories:
            categories = list(dge.tpms[compare_by].unique())
        if genes_to_show:
            if len(genes_to_show) * len(categories) > 40:
                st.write('Too many genes/categories to display, consider choosing fewer genes')
                st.stop()
            c3, c4 = st.columns(2)
            tpm_label = 'Normalised Counts'
            gene_df = dge.tpms[(dge.tpms[GENE_COL].isin(genes_to_show)) & (dge.tpms[compare_by].isin(categories))]
            if filter_out:
                gene_df = gene_df[gene_df[filter_by] == filter_out]
            groupby = st.radio('Group by', [GENE_COL, compare_by])
            color_by = [c for c in [GENE_COL, compare_by] if c != groupby][0]
            fig = px.box(gene_df, x=groupby, y=tpm_label, color=color_by, log_y=True,
                         hover_data=[GENE_COL] + list(dge.sd.columns), points='all')
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=16))
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
            st.plotly_chart(fig, use_container_width=True)

    st.header("DGE Results")
    with st.expander('Show DGE Results'):

        hits_df, fig = show_volcano(dge.results,  gene_name='Gene', genes_to_highlight=list(genes),
                               annotation_col='SYMBOL')

        c1, c2 = st.columns([3, 1])
        c1.plotly_chart(fig, use_container_width=True)
        c2.markdown("### STRING PPI Network")
        link_to_string(hits_df, c2)
        c2.markdown("### Download results")
        download_filtered_hits(hits_df, c2)

#
# if check_password():
#     app()

app()