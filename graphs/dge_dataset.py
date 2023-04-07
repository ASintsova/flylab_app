from processing.load_data import load_data, load_from_bucket, load_annotations
import pandas as pd
from sklearn.decomposition import PCA
import streamlit as st

class DgeDataSet:

    def __init__(self, bucket='', experiment_name='', data_dir='',  gene_name='',
                 sample_id = 'sampleID'):
        self.data_dir = data_dir
        self.bucket = bucket
        self.experiment_name = experiment_name
        self.gene_name = gene_name
        self.results, self.tpms, self.vsd, self.sd = self.load_experiment()
        self.sample_id = sample_id

    def load_experiment(self):
        if self.data_dir:
            return load_data(self.data_dir)
        else:
            return load_from_bucket(self.bucket,
                                    experiment_name=self.experiment_name,
                                    gene_name=self.gene_name)

    def load_annotations(self):
        annotations = load_annotations()['FLYBASE']
        annotations = {v:k for k,v in annotations.items()}
        self.tpms = self.tpms.reset_index()
        self.tpms = self.tpms.melt(id_vars='index', var_name=self.sample_id, value_name='Normalised Counts')
        self.tpms['SYMBOL'] = self.tpms['index'].map(annotations)
        self.tpms = self.tpms.rename({'index': 'FLYBASE'}, axis=1)
        self.tpms = self.tpms.merge(self.sd.reset_index(), left_on=self.sample_id, right_on=self.sample_id)

    @st.cache
    def find_pcs(self, num_pcs=2, num_genes=None, choose_by='variance'):
        """
        :param numPCs:
        :param numGenes:
        :return:
        """
        if num_genes:
            # calculate var for each, pick numGenes top var across samples -> df
            if choose_by == 'variance':
                genes = self.vsd.var(axis=1).sort_values(ascending=False).head(num_genes).index
                df = self.vsd.loc[genes].T
            else:
                pass
                # todo implement log2fc selection
        else:
            df = self.vsd.T
        pca = PCA(n_components=num_pcs)
        principal_components = pca.fit_transform(df)
        pcs = [f'PC{i}' for i in range(1, num_pcs + 1)]
        pc_df = (pd.DataFrame(data=principal_components, columns=pcs).set_index(df.index))
        pc_var = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, num_pcs)}
        pc_df = pc_df.merge(self.sd, left_index=True, right_index=True)
        return pc_df, pc_var

