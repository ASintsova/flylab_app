import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA


@st.cache
def find_PCs(countData, sampleData, numPCs=2, numGenes=None, choose_by='variance'):
    """
    :param countData: each column is a sampleID, index is featureID
    :param sampleData:
    :param numPCs:
    :param numGenes:
    :return:
    """
    if numGenes:
        # calculate var for each, pick numGenes top var across samples -> df
        if choose_by == 'variance':
            genes = countData.var(axis=1).sort_values(ascending=False).head(numGenes).index
            df = countData.loc[genes].T
        else:
            pass
            # todo implement log2fc selection
    else:
        df = countData.T
    pca = PCA(n_components=numPCs)
    principalComponents = pca.fit_transform(df)
    pcs = [f'PC{i}' for i in range(1, numPCs+1)]
    pDf = (pd.DataFrame(data=principalComponents, columns=pcs)
           .set_index(df.index))
    pc_var = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, numPCs)}
    pDf2 = pDf.merge(sampleData, left_index=True, right_index=True)
    return pDf2, pc_var