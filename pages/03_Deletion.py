import streamlit as st
from graphs.page_outline import RnaseqPage
from graphs.dge_dataset import DgeDataSet
from check_password import check_password

st.set_page_config(layout="wide")
page_name ='Deletion'


def app():
    st.title('Deletion')
    dge = DgeDataSet(bucket='jagannathan', experiment_name='deletion', sample_id='sampleID')
    dge.load_annotations()
    page = RnaseqPage(dge)
    page.pca_layout()
    page.dge_layout()
    #page.clustergram()

if check_password():
    app()