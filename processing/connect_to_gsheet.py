from google.oauth2 import service_account
import pandas as pd
from gsheetsdb import connect
import streamlit as st



def gsheet_to_df(url):
    #  Create a connection object.
    credentials = service_account.Credentials.from_service_account_info(
        st.secrets["gcp_service_account"],
        scopes=[
            "https://www.googleapis.com/auth/spreadsheets",
        ],
    )
    conn = connect(credentials=credentials)
    # Perform SQL query on the Google Sheet.
    # Uses st.cache to only rerun when the query changes or after 10 min.

    @st.cache(ttl=600)
    def run_query(query):
        rows = conn.execute(query, headers=1)
        rows = rows.fetchall()
        return rows

    return pd.DataFrame(run_query(f'SELECT * FROM "{url}"'))