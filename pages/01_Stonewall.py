import streamlit as st
from graphs.page_outline import RnaseqPage
from graphs.dge_dataset import DgeDataSet

st.set_page_config(layout="wide")
page_name ='Stonewall'


def app():
    st.title('Stonewall')
    dge = DgeDataSet(bucket='jagannathan', experiment_name='stwloe', sample_id = 'sampleID')
    dge.load_annotations()
    page = RnaseqPage(dge)
    page.pca_layout()
    page.dge_layout()

app()