from processing.load_data import load_data, load_from_bucket, load_annotations

class DgeDataSet:

    def __init__(self, bucket='', experiment_name='', data_dir='',  gene_name=''):
        self.data_dir = data_dir
        self.bucket = bucket
        self.experiment_name = experiment_name
        self.gene_name = gene_name
        self.results, self.tpms, self.vsd, self.sd = self.load_experiment()

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
        self.tpms = self.tpms.melt(id_vars='index', var_name='sampleID', value_name='Normalised Counts')
        self.tpms['SYMBOL'] = self.tpms['index'].map(annotations)
        self.tpms = self.tpms.rename({'index': 'FLYBASE'}, axis=1)
        self.tpms = self.tpms.merge(self.sd.reset_index(), left_on='sampleID', right_on='sampleID')
