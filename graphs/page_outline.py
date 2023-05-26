import streamlit as st
import plotly.express as px
import pandas as pd
import numpy as np
import requests
from time import sleep
from datetime import datetime

from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode, JsCode

class RnaseqPage:
    def __init__(self, dge):
        self.dge = dge
        self.count_data = self.dge.vsd
        self.tpm_data = self.dge.tpms
        self.sample_data = self.dge.sd
        self.results = self.dge.results

        self.config = {'contrast_col': 'contrast', 'pval_col': 'padj',
                       'lfc_col': 'log2FoldChange',
                       'gene_name': 'Gene',
                       'primary_annot': 'FLYBASE',
                       'secondary_annot': 'SYMBOL'}

        self.contrast_col = self.config['contrast_col']
        self.pval_col = self.config['pval_col']
        self.lfc_col = self.config['lfc_col']
        self.gene_name = self.config['gene_name']
        self.primary_annot = self.config['primary_annot']
        self.secondary_annot = self.config['secondary_annot']


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


    def get_genes(self):
        return list(self.results[self.gene_name].unique())


    def get_contrasts(self):
        return list(self.results[self.contrast_col].unique())


    def get_volcano_df(self, volcano_contrast):
        volcano_df = self.results[self.results[self.contrast_col] == volcano_contrast].copy()
        volcano_df['log10FDR'] = -1 * np.log10(volcano_df[self.pval_col])
        #volcano_df.loc[volcano_df['log10FDR'] > 20, 'log10FDR'] = 20
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


    def show_expression_graph(self, selected = []):
        if self.secondary_annot:
            annotation_column = st.radio('Choose gene annotation', [self.primary_annot, self.secondary_annot])
        else:
            annotation_column = self.primary_annot
        gene_options = self.tpm_data[annotation_column].unique()
        if selected:
            tpm_genes = selected
            df = self.tpm_data[self.tpm_data[self.primary_annot].isin(selected)].copy()
        else:
            tpm_genes = st.multiselect("Choose gene(s) of interest", gene_options, key='gois')
            df = self.tpm_data[self.tpm_data[annotation_column].isin(tpm_genes)].copy()
        df = df.apply(lambda x: np.log2(x + 0.5) if np.issubdtype(x.dtype, np.number) else x)
        c1, c2 = st.columns(2)
        compare_by = c1.selectbox('Compare by', self.sample_data.columns)
        categories = c1.multiselect(f'Categories of {compare_by} to display',
                                    ['All'] + list(self.sample_data[compare_by].unique()), key='cats')
        if len(self.sample_data.columns) > 1:
            filter_by = c2.selectbox("Filter by", self.sample_data.columns)
            filter_out = c2.selectbox(f'Which category of {filter_by} to keep?',
                                  [None] + list(self.sample_data[filter_by].unique()))
        else:
            filter_by = filter_out = ''
        if 'All' in categories:
            categories = list(self.sample_data[compare_by].unique())
        if tpm_genes:
            if len(tpm_genes) * len(categories) > 40:
                st.write('Too many genes/categories to display, consider choosing fewer genes')
                st.stop()
            c3, c4 = st.columns(2)
            tpm_label = 'log2 (Normalised Count)'
            #df = df[df[compare_by].isin(categories)]

            if filter_out:
                df = df[df[filter_by] != filter_out]
            # gene_df2 = (gene_df.melt(id_vars=['Gene'], value_name=tpm_label, var_name=sample_col)
            #             .merge(sample_df, how='inner', on=sample_col))
            groupby = st.radio('Group by', [annotation_column, compare_by])
            color_by = [c for c in [annotation_column, compare_by] if c != groupby][0]
            fig = px.box(df, x=groupby, y='Normalised Counts', color=color_by,
                         labels = {'Normalised Counts':'log2 Normalised Counts'},
                         hover_data=list(self.tpm_data.columns), points='all')
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=16))
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
            st.plotly_chart(fig, use_container_width=True)


    def link_to_string(self, hits_df, st_col):
        up = st_col.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'))
        if up == 'Upregulated Only':
            hits_df = hits_df[hits_df[self.lfc_col] > 0]
        elif up == 'Downregulated Only':
            hits_df = hits_df[hits_df[self.lfc_col] < 0]

        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = 'tsv-no-header'
        method = 'get_link'
        if self.gene_name:
            my_genes = set(hits_df[self.gene_name].values)
        else:
            my_genes = list(hits_df.index)
        request_url = "/".join([string_api_url, output_format, method])
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


    def convert_df(_self, df):
        return df.to_csv().encode('utf-8')

    def download_filtered_hits(self, hits_df, st_col, contrast_col="contrast"):

        fname_default = hits_df[contrast_col].unique()[0]
        fname = st_col.text_input("File name", value=fname_default)
        fname = fname + ".csv"
        st_col.download_button("Download data as csv file", self.convert_df(hits_df), file_name=fname)

    def dge_layout(self, annotation_col='SYMBOL'):

        ### Nice stuff comes by playing with the arguments the package has (https://streamlit-aggrid.readthedocs.io/en/docs/)
        # @st.cache
        def load_gridOptions(df):
            """Cached function to load grid options"""
            #df['Link'] = df['Gene'].apply(lambda x: f"http://flybase.org/reports/{x}")
            gb = GridOptionsBuilder.from_dataframe(df)
            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)  # Add pagination
            gb.configure_selection('multiple',
                                   use_checkbox=True)  # Enable multiselection, if nested rows add this -->  groupSelectsChildren="Group checkbox select children"
            gb.configure_column("Gene", headerName="Gene", cellRenderer=JsCode(
                '''function(params) {return `<a href=http://flybase.org/reports/${params.value} target="_blank">${params.value}</a>`}'''),
                                width=300)

            # gb.configure_column(
            #     "Link", "Link",
            #     cellRenderer=JsCode("""
            #         class UrlCellRenderer {
            #           init(params) {
            #             this.eGui = document.createElement('a');
            #             this.eGui.innerText = params.value;
            #             this.eGui.setAttribute('href', params.value);
            #             this.eGui.setAttribute('style', "text-decoration:none");
            #             this.eGui.setAttribute('target', "_blank");
            #           }
            #           getGui() {
            #             return this.eGui;
            #           }
            #         }
            #     """)
            # )
            gridOptions = gb.build()
            return gridOptions

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
            st.write(f"Number of hits: {volcano_df['Hit'].sum()}")
            volcano_short = volcano_df[['Gene', 'SYMBOL', 'baseMean', 'log2FoldChange', 'padj', 'contrast']]
            hits_df = volcano_df[volcano_df['Hit'] == True]
            c1, c2 = st.columns(2)
            self.link_to_string(hits_df, c1)
            self.download_filtered_hits(hits_df, c2)
            fig = self.volcano_graph(volcano_df, size_max, volcano_contrast, lfc_th, fdr, on_hover)

            st.plotly_chart(fig, use_container_width=True)
            if 'grid_key' not in st.session_state:
                st.session_state['grid_key'] = datetime.now()
                st.session_state['volcano_df'] = volcano_short
            if st.button(label='Update table'):
                st.session_state['volcano_df'] = volcano_short
                st.session_state['grid_key'] = datetime.now()

            grid_response = AgGrid(
                volcano_short,
                gridOptions=load_gridOptions(volcano_short),
                data_return_mode='AS_INPUT',
                #update_mode='NO_UPDATE',
                update_mode='MODEL_CHANGED',
                # header_checkbox_selection_filtered_only=True,
                fit_columns_on_grid_load=False,
                # theme='blue', #Add theme color to the table
                enable_enterprise_modules=True,
                height=350,
                width='100%',
                key=str(st.session_state['grid_key']),
                reload_data=False,
                allow_unsafe_jscode=True

            )
            selected = grid_response['selected_rows']
            if selected:
                #selected = pd.DataFrame(selected)[['Gene', 'SYMBOL','contrast', self.lfc_col, 'padj']]
                selected = [s['Gene'] for s in selected]
                #selected = [s[0] for s in selected]
                #df2 = selected.set_index('Gene').merge(self.count_data, left_index=True, right_index =True, how='left')  # Pass the selected rows to a new dataframe df
                #df3 = self.results[self.results[self.gene_name].isin(selected[self.gene_name].values)]
                #st.write(df2.reset)
                #st.write(df3)

            self.show_expression_graph(selected)
            st.markdown("## Show clustergram (under construction)")
            self.clustergram()

    def clustergram(self):
        import dash_bio as dashbio
        genes = self.count_data.var(axis=1).sort_values(ascending=False).head(500).index

        test = self.count_data.loc[genes]

        test = (test.merge(self.tpm_data[['FLYBASE', 'SYMBOL']].drop_duplicates(), left_index=True, right_on='FLYBASE', how='left')
                .set_index(['SYMBOL', 'FLYBASE']))
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
                    height=1000,
                    width=1200,
                    hidden_labels='row'
                )
        st.plotly_chart(fig)


            # c2.markdown("### STRING PPI Network")
            # link_to_string(hits_df, c2)
            # c2.markdown("### Download results")
            # download_filtered_hits(hits_df, c2)