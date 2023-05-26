import streamlit as st
from pathlib import Path
import pandas as pd
from google.oauth2 import service_account
from google.cloud import storage
from io import StringIO
import json


@st.cache_data
def load_data(datadir, gene_name=''):
    result_files = list(Path(datadir).glob("*results*csv"))
    tpms = list(Path(datadir).glob("*tpms*csv"))
    if not tpms:
        tpms = list(Path(datadir).glob("*norm_cnts*csv"))
    vsd = list(Path(datadir).glob("*vsd*csv"))
    sampleData = list(Path(datadir).glob("*metadata.csv"))
    results = pd.concat([pd.read_csv(f, index_col=0).assign(contrast=f.stem.split("_unfiltered")[0])
                         for f in result_files]) if result_files else pd.DataFrame()
    if not results.empty:
        if not gene_name:
            gene_name = 'Gene'
            results.index.name = gene_name
        results = results.reset_index()
    tpms = pd.read_csv(tpms[0], index_col=0) if tpms else None
    vsd = pd.read_csv(vsd[0], index_col=0) if vsd else None
    sd = pd.read_csv(sampleData[0], index_col=0, dtype=str) if sampleData else None
    return results, tpms, vsd, sd





def load_from_bucket(bucket_name, experiment_name, gene_name):
    # Create API client.
    credentials = service_account.Credentials.from_service_account_info(
        st.secrets["gcp_service_account"]
    )
    client = storage.Client(credentials=credentials)

    @st.cache_data(ttl=600)
    def read_file(file_path):
        bucket = client.bucket(bucket_name)
        content = bucket.blob(file_path).download_as_string().decode("utf-8")
        return pd.read_csv(StringIO(content), sep=',', index_col=0)

    file_blobs = list(client.list_blobs(bucket_name))

    result_files = [blob.name for blob in file_blobs if 'results.csv' in str(blob) and experiment_name in str(blob)]
    tpms_file = [blob.name for blob in file_blobs if 'tpms.csv' in str(blob) and experiment_name in str(blob)]
    if not tpms_file:
        tpms_file = [blob.name for blob in file_blobs if 'norm_cnts.csv' in str(blob) and experiment_name in str(blob)]
    vsd_file = [blob.name for blob in file_blobs if 'vsd.csv' in str(blob) and experiment_name in str(blob)]
    sampleData_file = [blob.name for blob in file_blobs if 'metadata.csv' in str(blob) and experiment_name in str(blob)]
    tpms = read_file(tpms_file[0]) if tpms_file else pd.DataFrame()
    vsd = read_file(vsd_file[0]) if vsd_file else pd.DataFrame()
    sd = read_file(sampleData_file[0]) if sampleData_file else pd.DataFrame()

    results = pd.concat([read_file(f).assign(contrast=f.split('-')[-1])
                         for f in result_files]) if result_files else pd.DataFrame()
    if not sd.empty:
        sd = sd.astype(str)

    if not results.empty:
        if not gene_name:
            gene_name = 'Gene'
            results.index.name = gene_name
        results = results.reset_index()

    return results, tpms, vsd, sd


def load_annotations():
    credentials = service_account.Credentials.from_service_account_info(
        st.secrets["gcp_service_account"]
    )
    client = storage.Client(credentials=credentials)

    @st.cache_data(ttl=600)
    def read_file(file_path):
        bucket = client.bucket('jagannathan')
        content = bucket.blob(file_path).download_as_string().decode("utf-8")
        return json.loads(content)

    return read_file('symbol_to_fb.json')

    #pms = read_file(client, bucket_name, tpms_file[0])
    #st.write(tpms.head())
    #content = read_file(client, bucket_name, file_path)

# bucket_name = "vorholt"
# file_path = "AK06/AK06_Sample_metadata.csv"
#
#
# df = pd.read_csv(content, sep=',')
#
# # Print results.
# st.write(df)