import os
import itertools

import epiweeks
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash.dependencies


# -------------------------- PARAMETERS ---------------------------- #

reportable_vocs = ('B.1.1.7', 'B.1.351', 'P.1')
min_unambig = 24000


# -------------------------- LOAD DATA ---------------------------- #

assemblies_tsv = os.path.join('data/assemblies.tsv')
collab_tsv = os.path.join('data/collab_ids.tsv')

df_assemblies = pd.read_csv(assemblies_tsv, sep='\t')
if collab_tsv and os.path.isfile(collab_tsv):
    collab_ids = pd.read_csv(collab_tsv, sep='\t')[list(['external_id', 'collaborator_id'])]
    collab_ids.columns = ['sample', 'collaborator_id']
else:
    collab_ids = pd.DataFrame(columns = ['sample', 'collaborator_id']) 


# -------------------------- TRANSFORM DATA ---------------------------- #

# format dates properly
df_assemblies = df_assemblies.astype({'collection_date':np.datetime64,'run_date':np.datetime64})

# fix missing data in purpose_of_sequencing
df_assemblies.loc[:,'purpose_of_sequencing'] = df_assemblies.loc[:,'purpose_of_sequencing'].fillna('Missing').replace('', 'Missing')

# derived columns: collection_epiweek, run_epiweek
df_assemblies.loc[:,'collection_epiweek'] = list(epiweeks.Week.fromdate(x) if not pd.isna(x) else x for x in df_assemblies.loc[:,'collection_date'])
df_assemblies.loc[:,'run_epiweek'] = list(epiweeks.Week.fromdate(x) if not pd.isna(x) else x for x in df_assemblies.loc[:,'run_date'])
df_assemblies.loc[:,'collection_epiweek_end'] = list(x.enddate().strftime('%Y-%m-%d') if not pd.isna(x) else '' for x in df_assemblies.loc[:,'collection_epiweek'])
df_assemblies.loc[:,'run_epiweek_end'] = list(x.enddate().strftime('%Y-%m-%d') if not pd.isna(x) else '' for x in df_assemblies.loc[:,'run_epiweek'])
df_assemblies.loc[:,('collection_date', 'collection_epiweek', 'collection_epiweek_end', 'run_date', 'run_epiweek', 'run_epiweek_end')]

# derived column: sample_age_at_runtime
df_assemblies.loc[:,'sample_age_at_runtime'] = list(x.days for x in df_assemblies.loc[:,'run_date'] - df_assemblies.loc[:,'collection_date'])

# derived column: genome_status
df_assemblies.loc[:,'genome_status'] = list(
        'failed_sequencing' if df_assemblies.loc[id, 'assembly_length_unambiguous'] < min_unambig
        else 'failed_annotation' if df_assemblies.loc[id, 'vadr_num_alerts'] > 0
        else 'submittable'
        for id in df_assemblies.index)

# derived columns: geo_country, geo_state, geo_locality
df_assemblies.loc[:,'geo_country'] = list(g.split(': ')[0] if not pd.isna(g) else '' for g in df_assemblies.loc[:,'geo_loc_name'])
df_assemblies.loc[:,'geo_state'] = list(g.split(': ')[1].split(', ')[0] if not pd.isna(g) else '' for g in df_assemblies.loc[:,'geo_loc_name'])
df_assemblies.loc[:,'geo_locality'] = list(g.split(': ')[1].split(', ')[1] if not pd.isna(g) and ', ' in g else '' for g in df_assemblies.loc[:,'geo_loc_name'])

# join column: collaborator_id
df_assemblies = df_assemblies.merge(collab_ids, on='sample', how='left', validate='one_to_one')


# get lists
states_all = list(x for x in df_assemblies['geo_state'].unique() if x)
collaborators_all = list(x for x in df_assemblies['collected_by'].unique() if x)
purposes_all = list(x for x in df_assemblies['purpose_of_sequencing'].unique() if x)


# -------------------------- DASH ---------------------------- #

external_stylesheets = [dbc.themes.BOOTSTRAP]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets, assets_folder='assets')
server = app.server

app.config.suppress_callback_exceptions = True


# -------------------------- VISUAL COMPONENTS ---------------------------- #

def get_subset(states, collabs, purpose):
    df = df_assemblies[df_assemblies['geo_state'].isin(states)]
    df = df[df['collected_by'].isin(collabs)]
    df = df[df['purpose_of_sequencing'].isin(purpose)]
    return df

@app.callback(
    dash.dependencies.Output('collab_selector', 'options'),
    dash.dependencies.Output('collab_selector', 'value'),
    dash.dependencies.Input('state_selector', 'value'),
    )
def collab_selector(states):
    df = df_assemblies[df_assemblies['geo_state'].isin(states)]
    collabs = list(x for x in df['collected_by'].unique() if x)
    return list({'label': x, 'value': x} for x in collabs), collabs

@app.callback(
    dash.dependencies.Output('fig_sequencing_success', 'figure'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    )
def fig_sequencing_success(states, collabs, purpose):
    ''' plot sequencing success rates '''
    df = get_subset(states, collabs, purpose)
    return px.histogram(df,
        title='Sequencing success by run date epiweek end',
        x='run_epiweek_end', color='genome_status',
        color_discrete_sequence=px.colors.qualitative.Pastel,
        category_orders={'genome_status':['submittable', 'failed_annotation', 'failed_sequencing']}
        )

@app.callback(
    dash.dependencies.Output('fig_sample_age', 'figure'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    )
def fig_sample_age(states, collabs, purpose):
    ''' plot collection-vs-sequencing dates '''
    df = get_subset(states, collabs, purpose)
    return px.scatter(df,
        title='Sequencing date vs collection date',
        x='collection_date', y='run_date', opacity=0.5,
        color='sample_age_at_runtime',
        color_continuous_scale=px.colors.sequential.Plasma,
        hover_data=['sample', 'collected_by', 'collaborator_id', 'biosample_accession', 'collection_date', 'run_date', 'sample_age_at_runtime', 'purpose_of_sequencing', 'pango_lineage', 'nextclade_clade', 'geo_loc_name']
        )

@app.callback(
    dash.dependencies.Output('fig_root_to_tip', 'figure'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    dash.dependencies.Input('selector_rtt_colorby', 'value'),
    )
def fig_root_to_tip(states, collabs, purpose, colorby):
    ''' Genetic distance as a QC check '''
    df = get_subset(states, collabs, purpose)
    df_good = df.query('genome_status != "failed_sequencing"')
    return px.scatter(df_good,
        title='Genetic distance root-to-tip vs sample collection date',
        x='collection_date', y='dist_to_ref_snps',
        color=colorby, opacity=0.7,
        hover_data=['sample', 'collaborator_id', 'biosample_accession', 'collected_by', 'collection_date', 'run_date', 'sample_age_at_runtime', 'purpose_of_sequencing', 'pango_lineage', 'nextclade_clade', 'geo_loc_name']
        )

@app.callback(
    dash.dependencies.Output('fig_nextclade_over_time', 'figure'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    )
def fig_nextclade_over_time(states, collabs, purpose):
    ''' Nextclade histograms over time '''
    df = get_subset(states, collabs, purpose)
    df_good = df.query('genome_status != "failed_sequencing"')
    return px.histogram(df_good,
        title='Nextclade phylogenetic classifications vs sample collection date',
        x='collection_epiweek_end',
        color='nextclade_clade',
        color_discrete_sequence=px.colors.qualitative.Vivid
        )

@app.callback(
    dash.dependencies.Output('table_vocs', 'data'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    )
def table_vocs(states, collabs, purpose):
    # Report on major VoCs
    df = get_subset(states, collabs, purpose)
    df_good = df.query('genome_status != "failed_sequencing"')
    df_vocs = df_good[df_good['pango_lineage'].isin(reportable_vocs)]
    table = df_vocs.groupby(
        ["collection_epiweek_end", "pango_lineage"], as_index=False, dropna=False
        ).agg(n=('assembly_fasta', 'count')
        ).pivot(index='collection_epiweek_end', columns='pango_lineage', values='n'
        ).reset_index().rename(columns={'index':'collection_epiweek_end'})
    return table.to_dict('records')

@app.callback(
    dash.dependencies.Output('table_seq_by_week', 'data'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    )
def table_numbers_by_week(states, collabs, purpose):
    # Report on major VoCs
    df = get_subset(states, collabs, purpose)
    table = df.groupby(
        ["run_epiweek", "genome_status"], as_index=False, dropna=False
        ).agg(n=('assembly_fasta', 'count')
        ).pivot(index='run_epiweek', columns='genome_status', values='n'
        ).reset_index().rename(columns={'index':'run_epiweek'})
    table.loc[:,'run_epiweek_end'] = list(x.enddate().strftime('%Y-%m-%d') if not pd.isna(x) else '' for x in table.loc[:,'run_epiweek'])
    table['n_submittable'] = table['submittable']
    table['n_good'] = table['submittable'] + table['failed_annotation']
    table['n_attempted'] = table['submittable'] + table['failed_annotation'] + table['failed_sequencing']
    return table.to_dict('records')

@app.callback(
    dash.dependencies.Output('table_vocs_by_sample', 'data'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    dash.dependencies.Input('selector_sample_details', 'value'),
    )
def table_vocs_by_sample(states, collabs, purpose, sample_set):
    # Report on major VoCs
    df = get_subset(states, collabs, purpose)
    df_good = df.query('genome_status != "failed_sequencing"')
    if sample_set == 'vocs':
        df_vocs = df_good[df_good['pango_lineage'].isin(reportable_vocs)]
        return df_vocs.to_dict('records')
    else:
        return df_good.to_dict('records')


@app.callback(
    dash.dependencies.Output('basic_stats_card', 'children'),
    dash.dependencies.Input('state_selector', 'value'),
    dash.dependencies.Input('collab_selector', 'value'),
    dash.dependencies.Input('purpose_selector', 'value'),
    )
def basic_stats_card(states, collabs, purpose):
    df = get_subset(states, collabs, purpose)
    df_good = df.query('genome_status != "failed_sequencing"')
    df_vocs = df_good[df_good['pango_lineage'].isin(reportable_vocs)]

    return [
        html.P("Samples sequenced: {}".format(
            len(df.index)),
            className='card-text'),
        html.P("Genomes assembled: {}".format(
            len(df_good.index)),
            className='card-text'),
        html.P("Genomes released or in process for data release: {}".format(
            len(df_good.query('genome_status == "submittable"'))),
            className='card-text'),
        html.P("Samples with Variants of Concern (VoCs): {} ({})".format(
            len(df_vocs.index), ', '.join(df_vocs['pango_lineage'].unique())),
            className='card-text'),
        html.P("Samples sequenced between {} and {}.".format(
            df['run_date'].min().strftime('%Y-%m-%d'), df['run_date'].max().strftime('%Y-%m-%d')),
            className='card-text'),
        html.P("Samples collected between {} and {}.".format(
            df['collection_date'].min().strftime('%Y-%m-%d'), df['collection_date'].max().strftime('%Y-%m-%d')),
            className='card-text'),
        html.P("Samples originating from: {}".format(
            ', '.join(df['geo_state'].unique())),
            className='card-text'),
        html.P("Samples sent by: {}".format(
            ', '.join(df['collected_by'].unique())),
            className='card-text'),
    ]


# -------------------------- PROJECT DASHBOARD ---------------------------- #


app.layout = html.Div(children=[
    html.H1(
        children=[
            html.Div(
                id='banner',
                className='banner',
                children=[
                    html.Img(src="https://www.broadinstitute.org/files/news/media-images/logos/BroadInstLogoforDigitalRGB.png"),
                ],
            ),
            html.P(
                id='top',
                children='Summary of SARS-CoV-2 Sequencing by the Broad Institute'),
            ]
    ),

    html.Hr(),

    html.Div(id='body', className='container', children=[

        html.Div([html.P(
            """The Broad Institute Viral Genomics group, in partnership with the Broad Genomics Platform and
            Data Sciences Platform, has been engaged in viral sequencing of COVID-19 patients since March 2020.
            This dashboard describes samples sequenced since {} (we will backfill with more historical data
            later).""".format(df_assemblies['run_date'].min().strftime('%B %d, %Y'))
        )]),

        dbc.Card(dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    html.H3('Choose which states and labs to subset this report to')
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    html.Div([
                        html.Label('States of origin'),
                        dcc.Dropdown(
                            id='state_selector',
                            options=list({'label': x, 'value': x} for x in states_all),
                            value=states_all,
                            multi=True
                        ),
                    ]),
                ]),
                dbc.Col([
                    html.Div([
                        html.Label('Sample submitting lab'),
                        dcc.Dropdown(
                            id='collab_selector',
                            multi=True
                        ),
                    ]),
                ], width=8),
            ]),

            dbc.Row([
                dbc.Col([
                    html.Label('Purpose of sequencing'),
                    dcc.Dropdown(
                        id='purpose_selector',
                        options=list({'label': x, 'value': x} for x in purposes_all),
                        value=purposes_all,
                        multi=True
                    ),
                ], width=12),
            ]),

            #dbc.Row([dbc.Col([
            #    html.Label('Sequencing dates'),
            #    dcc.RangeSlider(id='run_date_range',
            #        min=df_assemblies['run_date'].min(), max=df_assemblies['run_date'].max(),
            #        value=[df_assemblies['run_date'].min(), df_assemblies['run_date'].max()]
            #    ),
            #])]),

            #dbc.Row([dbc.Col([
            #    html.Label('Collection dates'),
            #    dcc.RangeSlider(id='collection_date_range',
            #        min=df_assemblies['collection_date'].min(), max=df_assemblies['collection_date'].max(),
            #        value=[df_assemblies['collection_date'].min(), df_assemblies['collection_date'].max()]
            #    ),
            #])]),
        ])),

        html.Hr(),

        html.H3('Summary of selected data'),
        html.Div(id='basic_stats_card', className='container'),

        html.Br(),

        dbc.Card(dbc.CardBody([
            dcc.Graph(id="fig_sequencing_success"),
            html.P('''This describes the total number of patient samples sequenced in this data set,
            plotted by the date of the sequencing run in our lab (by CDC epiweek).'''),
            html.P('''"Submittable" genomes pass all QC checks and are quickly released to public genome repositories.
            "Failed sequencing" are samples that failed to produce at least {} unambiguous base pairs of
            viral genome. Raw data from these samples are submitted to NCBI's SRA database, but the genomes
            are not used for any analyses. "Failed annotation" are samples that produced a sufficiently complete
            genome, but did not pass NCBI's VADR quality checks.'''.format(min_unambig))
        ])),

        html.Br(),

        dbc.Card(dbc.CardBody([
            dcc.Graph(id="fig_sample_age"),
            html.P('''This plot describes the "timeliness" of the sequencing run for the purpose of
            real-time surveillance of circulating lineages and variants of interest. Note that this plot
            likely includes many samples that were sequenced for non-surveilance purposes.
            '''),
        ])),

        html.Br(),

        dbc.Card(dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    html.Div([
                        dcc.Graph(id="fig_root_to_tip"),
                    ]),
                ]),
                dbc.Col([
                    dcc.RadioItems(
                        id='selector_rtt_colorby',
                        options=[
                            {'label':'sequencing batch', 'value': 'flowcell_id'},
                            {'label':'sending lab', 'value': 'collected_by'},
                            {'label':'state', 'value': 'geo_state'},
                            {'label':'purpose of sequencing', 'value': 'purpose_of_sequencing'},
                            {'label':'nextclade clade', 'value': 'nextclade_clade'},
                        ],
                        value='flowcell_id'
                    ),
                ], width=3, align="center"),
            ]),
            dbc.Row([
                dbc.Col([
                    html.P('''A "root-to-tip plot" plots the genetic distance of each sample from Wuhan Hu-1
                    against the date it was collected. It is generally somewhat linear. Outliers on this plot
                    may be indicative of laboratory or metadata errors, or of evolutionarily unusual lineages
                    (such as B.1.1.7).'''),
                ], width=12, align="center"),
            ]),
        ])),

        html.Br(),

        dbc.Card(dbc.CardBody([
            dcc.Graph(id="fig_nextclade_over_time"),
            html.P('''This shows the breakdown of major phylogenetic clades over time, using the Nextclade
            naming system. Variants of Concern (VoCs) are highlighted as specially named Nextclade clades.
            Nextclade clade 20I/501Y.V1 corresponds to PANGO lineage B.1.1.7, 20H/501Y.V2 corresponds to B.1.351,
            and 20J/501Y.V3 corresponds to P.1.'''),
        ])),

        html.Br(),

        dbc.Card(dbc.CardBody([
            html.P(children='Sequencing activity by CDC epiweek of sequencing run date.'),
            html.Div(children=dash_table.DataTable(
                id='table_seq_by_week',
                columns=[
                    {'name':'sequencing epiweek', 'id':'run_epiweek'},
                    {'name':'sequencing epiweek end', 'id':'run_epiweek_end'},
                    {'name':'samples attempted', 'id':'n_attempted'},
                    {'name':'genomes assembled', 'id':'n_good'},
                    {'name':'genomes submittable', 'id':'n_submittable'},
                ],
                sort_action='native',
                export_format='xlsx',
                export_headers='names',
            ))
        ])),

        html.Br(),

        dbc.Card(dbc.CardBody([
            html.P(children='Reportable VoC counts by CDC epiweek of sample collection.'),
            html.Div(children=dash_table.DataTable(
                id='table_vocs',
                columns=[{'name':'sample collection epiweek end', 'id':'collection_epiweek_end'}
                        ] + [{'name':p,'id':p} for p in reportable_vocs],
                sort_action='native',
                export_format='xlsx',
                export_headers='ids',
            ))
        ])),

        html.Br(),

        dbc.Card(dbc.CardBody([
            html.P(children='Details for each sample.'),
            html.Div([
                dcc.RadioItems(
                    id='selector_sample_details',
                    options=[
                        {'label':'VoCs only', 'value': 'vocs'},
                        {'label':'all samples', 'value': 'all'},
                    ],
                    value='vocs'
                ),
            ]),
            html.Div(children=dash_table.DataTable(
                id='table_vocs_by_sample',
                columns=[
                    {'name':'sample', 'id':'sample'},
                    {'name':'collaborator_id', 'id':'collaborator_id'},
                    {'name':'biosample_accession', 'id':'biosample_accession'},
                    {'name':'pango_lineage', 'id':'pango_lineage'},
                    {'name':'nextclade_clade', 'id':'nextclade_clade'},
                    {'name':'collection_date', 'id':'collection_date'},
                    {'name':'geo_loc_name', 'id':'geo_loc_name'},
                    {'name':'run_date', 'id':'run_date'},
                    {'name':'assembly_length_unambiguous', 'id':'assembly_length_unambiguous'},
                    {'name':'amplicon_set', 'id':'amplicon_set'},
                    {'name':'vadr_num_alerts', 'id':'vadr_num_alerts'},
                    {'name':'collected_by', 'id':'collected_by'},
                    {'name':'purpose_of_sequencing', 'id':'purpose_of_sequencing'},
                    {'name':'bioproject_accession', 'id':'bioproject_accession'},
                ],
                sort_action='native',
                filter_action='native',
                page_action='native',
                page_size=20,
                export_format='xlsx',
                export_headers='ids',
                style_table={'overflowX': 'auto'},
            ))
        ])),
    ]),

    html.Br(),

    html.Footer([
        html.P(
            id='bottom',
            children='SARS-CoV-2 Sequencing by the Broad Institute. All data was generated at the Broad Institute in partnership with the listed collaborating institutions that provided patient samples and metadata.'),
    ]),

])


# -------------------------- MAIN ---------------------------- #


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8080, debug=True, use_reloader=False)