import streamlit as st
from processing import eda
import plotly.express as px
import pandas as pd


def show_pca(countData, sampleData):
    c1, c2 = st.columns((4, 1))
    c2.write('### PCA Options')
    maxComponents = min(countData.shape[0], countData.shape[1])
    numPCs = c2.slider("Select number of Principal Components", min_value=2, max_value=maxComponents, value=2)
    numGenes = c2.slider("Number of genes to use", value=500, max_value=countData.shape[0])
    choose_by = c2.selectbox('Choose genes based on highest', ['variance', 'log2FoldChange (not implemented)'])
    pDf, pc_var = eda.find_PCs(countData, sampleData, numPCs, numGenes, choose_by)
    pcX_labels = [f'PC{i}' for i in range(1, numPCs+1)]
    expVars = [c for c in pDf.columns if c not in pcX_labels]
    pcX = c2.selectbox('X-axis component', pcX_labels)
    pcY = c2.selectbox('Y-axis component', [pc for pc in pcX_labels if pc != pcX])
    pcVar = c2.radio('Color', expVars, key='c')
    pcSym = c2.radio('Symbol', [None] + expVars, key='s')
    fig = px.scatter(pDf, x=pcX, y=pcY, color=pcVar, symbol=pcSym,
                         labels ={pcX: f'{pcX}, {pc_var[pcX]} % Variance',
                                  pcY: f'{pcY}, {pc_var[pcY]} % Variance'},
                         height=700, hover_data=expVars, hover_name=pDf.index)
    fig.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
                          )
    fig.update_traces(marker=dict(size=12,
                                      line=dict(width=2,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
    c1.write(f'### {pcX} vs {pcY}, highlighting {pcVar}')
    c1.plotly_chart(fig, use_container_width=True)
    c3, c4 = st.columns(2)
    pDf_sum = pDf.groupby(pcVar).median()
    varDf = pd.DataFrame.from_dict(pc_var, orient='index').reset_index()
    varDf.columns = ['PC', '% Variance']
    fig2 = px.line(varDf, x='PC', y='% Variance', markers=True,
                   labels={'PC': ''})
    fig2.update_traces(marker=dict(size=12,
                                  line=dict(width=2,
                                            color='DarkSlateGrey')))
    c3.write('### Scree Plot')
    c3.plotly_chart(fig2)
    c4.write(f'### PCs summarized by {pcVar}')
    c4.plotly_chart(px.imshow(pDf_sum), use_container_width=True)