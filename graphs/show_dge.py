import streamlit as st
import plotly.express as px
import numpy as np
import requests
from time import sleep


@st.cache
def convert_df(df):
    return df.to_csv().encode('utf-8')

@st.cache
def get_contrasts(results, contrast_col):
    return list(results[contrast_col].unique())

@st.cache
def get_genes(results, gene_name):
    return list(results[gene_name].unique())

# @st.cache
# def get_symbols(results, gene_name, annotations):
#     return [a for a, b in annotations.items() if b in results[gene_name]]

@st.cache
def get_volcano_df(results, contrast_col, volcano_contrast, pval_col, ):
    volcano_df = results[results[contrast_col] == volcano_contrast].copy()
    volcano_df['log10FDR'] = -1 * np.log10(volcano_df[pval_col])
    volcano_df.loc[volcano_df['log10FDR'] > 20, 'log10FDR'] = 20
    return volcano_df


def show_volcano(results, gene_name='Gene',
                 contrast_col='contrast', lfc_col='log2FoldChange',
                 pval_col='padj',  genes_to_highlight=(), annotation_col='',):

    st.markdown('### Options')

    contrasts = get_contrasts(results, contrast_col)
    genes = get_genes(results, gene_name)
    genes_to_highlight = [g for g in genes_to_highlight if g in genes]
    volcano_contrast = st.selectbox('Select a contrast', contrasts, key='volcano_contrasts')
    volcano_df = get_volcano_df(results, contrast_col, volcano_contrast, pval_col).copy()
    k1, k2, k3 = st.columns(3)
    fdr = k1.number_input('FDR cutoff', value=0.05)
    lfc_th = k2.number_input('Log FC cutoff (absolute)', value=1.0)
    volcano_df['Hit'] = ((abs(volcano_df[lfc_col]) > lfc_th) & (volcano_df[pval_col] < fdr))
    genes_to_highlight = k3.multiselect("Choose gene(s) of interest", genes,
                                             default=genes_to_highlight, key='higenes')
    if annotation_col:
        volcano_df[annotation_col] = volcano_df[annotation_col].fillna('N/A')
    if genes_to_highlight:
        volcano_df['GOI'] = (volcano_df[gene_name].isin(genes_to_highlight).astype(int) * 50 + 10).astype(float)
        size_max = 25
    else:
        volcano_df['GOI'] = 10
        size_max = 10

    on_hover = {lfc_col: True,
                pval_col: ':.3f',
                'log10FDR': False,
                'Hit': False,
                'GOI': False,
                                 }
    if annotation_col:
        on_hover[annotation_col] = True
    fig = px.scatter(volcano_df, x=lfc_col, y='log10FDR', color='Hit', size='GOI',
                     height=700, size_max=size_max,
                     title=volcano_contrast,
                     color_discrete_map={False: '#fff9f0',
                                         True: '#378b84'},
                     category_orders={'Hit': [False, True], 'in Pathway': [False, True]},
                     hover_name=volcano_df[gene_name],
                     hover_data=on_hover)


    fig.add_vline(x=lfc_th, line_width=2, line_dash="dash", line_color="grey")
    fig.add_vline(x=-lfc_th, line_width=2, line_dash="dash", line_color="grey")
    fig.add_hline(y=-1 * np.log10(fdr), line_width=2, line_dash="dash", line_color="grey")
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                      font=dict(size=18))
    fig.update_traces(marker=dict(
                                  line=dict(width=1,
                                            color='DarkSlateGrey'), opacity=0.8),
                      selector=dict(mode='markers'))
    return volcano_df[volcano_df['Hit'] == True], fig


def link_to_string(hits_df, st_col, lfc_col='log2FoldChange', gene_name='Gene'):
    up = st_col.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'))
    if up == 'Upregulated Only':
        hits_df = hits_df[hits_df[lfc_col] > 0]
    elif up == 'Downregulated Only':
        hits_df = hits_df[hits_df[lfc_col] < 0]

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = 'tsv-no-header'
    method = 'get_link'
    if gene_name:
        my_genes = set(hits_df[gene_name].values)
    else:
        my_genes = list(hits_df.index)
    request_url = "/".join([string_api_url, output_format, method])
    #species = st_col.number_input("NCBI species taxid", value=3702, help='Arabidopsis thaliana: 3702')
    species = 7227
    params = {
        "identifiers": "\r".join(my_genes),  # your protein
        "species": species,  # species NCBI identifier
        "network_flavor": "confidence",  # show confidence links
        "caller_identity": "explodata"  # your app name
    }
#
    if st_col.button('Get STRING network'):
        network = requests.post(request_url, data=params)
        network_url = network.text.strip()
        st_col.markdown(f"[Link to STRING network]({network_url})")
        sleep(1)


def download_filtered_hits(hits_df, st_col, contrast_col="contrast"):

    fname_default = hits_df[contrast_col].unique()[0]
    fname = st_col.text_input("File name", value=fname_default)
    fname = fname+".csv"
    st_col.download_button("Download data as csv file", convert_df(hits_df), file_name=fname)


# def compare_between_contrasts():
#     contrasts = results[contrast_col].unique()
#     compContrasts = st.multiselect('Select contrasts to display', contrasts)
#     c1, c2 = st.columns(2)
#     filters = {}
#     for col, con in zip(cycle([c1, c2]), compContrasts):
#         col.write(con)
#         l = col.number_input('LFC cutoff', value=-1.0, step=0.5, key=f'{con}_lfc')
#         f = col.number_input('FDR cutoff', value=0.05, step=0.01, key=f'{con}_fdr')
#         filters[con] = (l, f)
#     comp_genes = []
#
#     for key, value in filters.items():
#         if value[0] > 0:
#             genes = set(vennDf[(vennDf[contrast_col] == key) & (vennDf[lfc] > value[0])
#                                & (vennDf[pval] < value[1])][gene_name].values)
#         else:
#             genes = set(vennDf[(vennDf[contrast_col] == key) & (vennDf[lfc] < value[0])
#                                & (vennDf[pval] < value[1])][gene_name].values)
#         comp_genes.append(genes)
#     if not comp_genes:
#         st.stop()
#     intersect_genes = set.intersection(*comp_genes)
#     vennDf = vennDf[vennDf[gene_name].isin(intersect_genes)].copy()
#     vennDf = vennDf[[gene_name, lfc, pval, contrast_col]].drop_duplicates()
#     # vennDf2 = vennDf.pivot(index='Name', columns='contrast_col', values=['LFC', 'fdr'])
#     st.write(vennDf.shape)