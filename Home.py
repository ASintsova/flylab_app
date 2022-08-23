import streamlit as st
from google.oauth2 import service_account
from google.cloud import storage
from io import StringIO
from processing import load_data

import pandas as pd
from check_password import check_password


st.set_page_config(page_title="Jagannathan Lab", layout='wide')


def home_page():
    st.title("RNA-sequencing projects")

if check_password():
    home_page()



