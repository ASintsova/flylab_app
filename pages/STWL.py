import streamlit as st
from processing.load_data import load_data, load_from_bucket, load_annotations
from graphs.show_pca import show_pca
from graphs.show_gene_expression import show_expression
from graphs.show_dge import show_volcano, link_to_string, download_filtered_hits
from check_password import check_password

page_name ='STWL OE'


def app():
    st.title('STWL')

    results, tpms, vsd, sd = load_from_bucket('jagannathan', experiment_name='stwloe', gene_name='')
    st.header("PCA")
    with st.expander('Show PCA'):
        show_pca(vsd, sd)

    st.header("Expression")
    with st.expander('Show Gene Expression'):
        annotations = load_annotations()['FLYBASE']
        st.write(annotations['stwl'])
        genes = show_expression(tpms, sd, annotations, sample_col='sampleID')

    st.header("DGE Results")
    with st.expander('Show DGE Results'):
        c1, c2, c3 = st.columns([1, 3, 1])
        hits_df = show_volcano(results, c1, c2, gene_name='Gene', genes_to_highlight=list(genes),
                               annotation_col='SYMBOL')
        c3.markdown("### STRING PPI Network")
        link_to_string(hits_df, c3)
        c3.markdown("### Download results")
        download_filtered_hits(hits_df, c3)


if check_password():
    app()