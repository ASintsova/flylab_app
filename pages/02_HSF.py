import streamlit as st
from graphs.page_outline import RnaseqPage
from graphs.dge_dataset import DgeDataSet
from check_password import check_password


st.set_page_config(layout="wide")
page_name ='HSF'


def app():
    st.title('HSF')
    dge = DgeDataSet(bucket='jagannathan', experiment_name='hsf', sample_id = 'sampleID')
    dge.load_annotations()
    page = RnaseqPage(dge)
    page.pca_layout()
    page.dge_layout()

if check_password():
    app()