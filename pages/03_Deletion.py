import streamlit as st
from graphs.page_outline import RnaseqPage
from graphs.dge_dataset import DgeDataSet

st.set_page_config(layout="wide")
page_name ='Deletion'


def app():
    st.title('Deletion')
    dge = DgeDataSet(data_dir='/Users/ansintsova/git_repos/fly_rnaseq/data/deletion/results')
    dge.load_annotations()
    page = RnaseqPage(dge)
    page.pca_layout()
    page.dge_layout()
    page.clustergram()

app()