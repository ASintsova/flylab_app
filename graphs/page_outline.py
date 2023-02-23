import streamlit as st
import plotly.express as px
import pandas as pd
import numpy as np
from streamlit_plotly_events import plotly_events

class RnaseqPage:
    def __init__(self, dge):
        self.dge = dge
        self.count_data = self.dge.vsd
        self.sample_data = self.dge.sd
        self.results = self.dge.results

        self.config = {'contrast_col': 'contrast', 'pval_col': 'padj',
                       'lfc_col': 'log2FoldChange',
                       'gene_name': 'Gene'}

        self.contrast_col = self.config['contrast_col']
        self.pval_col = self.config['pval_col']
        self.lfc_col = self.config['lfc_col']
        self.gene_name = self.config['gene_name']

    @st.cache
    def convert_df(self, df):
        return df.to_csv().encode('utf-8')

    #@st.cache
    def pca_graph(self, pc_df, pc_x, pc_y, pc_var, pc_col, pc_sym, exp_vars):
        pc_df_sum = pc_df.groupby(pc_col).median()
        var_df = pd.DataFrame.from_dict(pc_var, orient='index').reset_index()
        var_df.columns = ['PC', '% Variance']
        fig = px.scatter(pc_df, x=pc_x, y=pc_y, color=pc_col, symbol=pc_sym,
                         labels={pc_x: f'{pc_x}, {pc_var[pc_x]} % Variance',
                                 pc_y: f'{pc_y}, {pc_var[pc_y]} % Variance'},
                         height=700, hover_data=exp_vars, hover_name=pc_df.index)
        fig.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
                          )
        fig.update_traces(marker=dict(size=12,
                                      line=dict(width=2,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))

        fig2 = px.line(var_df, x='PC', y='% Variance', markers=True,
                       labels={'PC': ''})
        fig2.update_traces(marker=dict(size=12,
                                       line=dict(width=2,
                                                 color='DarkSlateGrey')))
        fig3 = px.imshow(pc_df_sum)
        return fig, fig2, fig3

    def pca_layout(self):
        st.header("PCA")
        with st.expander('Show PCA'):
            c1, c2 = st.columns((4, 1))
            c2.write('### PCA Options')
            max_components = min(self.count_data.shape[0], self.sample_data.shape[0])
            num_pcs = c2.slider("Select number of Principal Components", min_value=2, max_value=max_components, value=2)
            num_genes = c2.slider("Number of genes to use", value=500, max_value=self.count_data.shape[0])
            choose_by = c2.selectbox('Choose genes based on highest', ['variance', 'log2FoldChange (not implemented)'])
            pc_df, pc_var = self.dge.find_pcs(num_pcs, num_genes, choose_by)
            pc_x_labels = [f'PC{i}' for i in range(1, num_pcs + 1)]
            exp_vars = [c for c in pc_df.columns if c not in pc_x_labels]
            pc_x = c2.selectbox('X-axis component', pc_x_labels)
            pc_y = c2.selectbox('Y-axis component', [pc for pc in pc_x_labels if pc != pc_x])
            pc_col = c2.radio('Color', exp_vars, key='c')
            pc_sym = c2.radio('Symbol', [None] + exp_vars, key='s')
            fig1, fig2, fig3 = self.pca_graph(pc_df, pc_x, pc_y, pc_var, pc_col, pc_sym, exp_vars)
            c1.write(f'### {pc_x} vs {pc_y}, highlighting {pc_col}')
            c1.plotly_chart(fig1, use_container_width=True)
            c3, c4 = st.columns(2)
            c3.write('### Scree Plot')
            c3.plotly_chart(fig2)
            c4.write(f'### PCs summarized by {pc_col}')
            c4.plotly_chart(fig3, use_container_width=True)

    @st.cache
    def get_genes(self):
        return list(self.results[self.gene_name].unique())

    @st.cache
    def get_contrasts(self):
        return list(self.results[self.contrast_col].unique())


    def get_volcano_df(self, volcano_contrast):
        volcano_df = self.results[self.results[self.contrast_col] == volcano_contrast].copy()
        volcano_df['log10FDR'] = -1 * np.log10(volcano_df[self.pval_col])
        volcano_df.loc[volcano_df['log10FDR'] > 20, 'log10FDR'] = 20
        return volcano_df

    def volcano_graph(self, volcano_df,  size_max, volcano_contrast, lfc_th, fdr, on_hover):
        fig = px.scatter(volcano_df, x=self.lfc_col, y='log10FDR', color='Hit', size='GOI',
                         height=700,
                         size_max=size_max,
                         title=volcano_contrast,
                         color_discrete_map={False: '#fff9f0',
                                         True: '#378b84'},
                         category_orders={'Hit': [False, True], 'in Pathway': [False, True]},
                         hover_name=volcano_df[self.gene_name],
                         hover_data=on_hover)

        fig.add_vline(x=lfc_th, line_width=2, line_dash="dash", line_color="grey")
        fig.add_vline(x=-lfc_th, line_width=2, line_dash="dash", line_color="grey")
        fig.add_hline(y=-1 * np.log10(fdr), line_width=2, line_dash="dash", line_color="grey")
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=14))
        fig.update_traces(marker=dict(
            line=dict(width=0.8,
                      color='DarkSlateGrey'), opacity=0.9),
            selector=dict(mode='markers'))
        return fig

    # def link_to_string(self, hits_df, st_col):
    #     up = st_col.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'))
    #     if up == 'Upregulated Only':
    #         hits_df = hits_df[hits_df[lfc_col] > 0]
    #     elif up == 'Downregulated Only':
    #         hits_df = hits_df[hits_df[lfc_col] < 0]
    #
    #     string_api_url = "https://version-11-5.string-db.org/api"
    #     output_format = 'tsv-no-header'
    #     method = 'get_link'
    #     if gene_name:
    #         my_genes = set(hits_df[gene_name].values)
    #     else:
    #         my_genes = list(hits_df.index)
    #     request_url = "/".join([string_api_url, output_format, method])
    #     # species = st_col.number_input("NCBI species taxid", value=3702, help='Arabidopsis thaliana: 3702')
    #     species = 7227
    #     params = {
    #         "identifiers": "\r".join(my_genes),  # your protein
    #         "species": species,  # species NCBI identifier
    #         "network_flavor": "confidence",  # show confidence links
    #         "caller_identity": "explodata"  # your app name
    #     }
    #     #
    #     if st_col.button('Get STRING network'):
    #         network = requests.post(request_url, data=params)
    #         network_url = network.text.strip()
    #         st_col.markdown(f"[Link to STRING network]({network_url})")
    #         sleep(1)

    def dge_layout(self, annotation_col='SYMBOL'):
        st.header("DGE Results")
        with st.expander('Show DGE Results'):
            st.markdown('### Options')
            contrasts = self.get_contrasts()
            genes = self.get_genes()
            volcano_contrast = st.selectbox('Select a contrast', contrasts, key='volcano_contrasts')
            volcano_df = self.get_volcano_df(volcano_contrast)
            k1, k2, k3 = st.columns(3)
            fdr = k1.number_input('FDR cutoff', value=0.05)
            lfc_th = k2.number_input('Log FC cutoff (absolute)', value=1.0)
            volcano_df['Hit'] = ((abs(volcano_df[self.lfc_col]) > lfc_th) & (volcano_df[self.pval_col] < fdr))
            genes_to_highlight = k3.multiselect("Choose gene(s) of interest", genes, key='higenes')
            if annotation_col:
                volcano_df[annotation_col] = volcano_df[annotation_col].fillna('N/A')
            if genes_to_highlight:
                volcano_df['GOI'] = (volcano_df[self.gene_name].isin(genes_to_highlight).astype(int) * 50 + 10).astype(
                    float)
                size_max = 25
            else:
                volcano_df['GOI'] = 10
                size_max = 10

            on_hover = {self.lfc_col: True,
                        self.pval_col: ':.3f',
                        'log10FDR': False,
                        'Hit': False,
                        'GOI': False,
                        }
            if annotation_col:
                on_hover[annotation_col] = True

            hits_df = volcano_df[volcano_df['Hit'] == True]
            fig = self.volcano_graph(volcano_df,  size_max, volcano_contrast, lfc_th, fdr, on_hover)
            c1, c2 = st.columns([3, 1])
            st.plotly_chart(fig, use_container_width=True)

    def clustergram(self):
        import dash_bio as dashbio
        genes = self.count_data.var(axis=1).sort_values(ascending=False).head(500).index
        test = self.count_data.loc[genes]


        rows = list(test.index.values)
        columns = list(test.columns)
        fig = dashbio.Clustergram(
                    data=test.loc[rows].values,
                    row_labels=rows,
                    column_labels = columns,
                standardize = 'row',
                    color_threshold={
                        'row': 250,
                        'col': 700
                    },
                    height=800,
                    width=700,
                    hidden_labels='row'
                )
        st.plotly_chart(fig)


            # c2.markdown("### STRING PPI Network")
            # link_to_string(hits_df, c2)
            # c2.markdown("### Download results")
            # download_filtered_hits(hits_df, c2)