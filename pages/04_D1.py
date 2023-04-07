import streamlit as st
from graphs.page_outline import RnaseqPage
from graphs.dge_dataset import DgeDataSet

st.set_page_config(layout="wide")
page_name ='Deletion'


def app():
    st.title('D1')
    dge = DgeDataSet(bucket='jagannathan', experiment_name='07-04-23_d1', sample_id='sample_id')
    dge.load_annotations()
    page = RnaseqPage(dge)
    page.pca_layout()
    page.dge_layout()
    #page.clustergram()

app()