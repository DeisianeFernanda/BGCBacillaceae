"""This code takes the Antismash summary data and compares to Gecco. Is still wonder of construction :)
So problably this will not run in your computer yet :("""

import pandas as pd
from pathlib import Path
from IPython.display import display, Markdown, HTML
import json
import altair as alt
import numpy as np

# Function to find overlaps between query regions and lookup targets
def find_overlap(query, lookup):
    overlap_key = None
    overlap_length = 0
    query_length = query['end'] - query['start']
    target_length = 0
    overlap_proportion = 0
    is_overlapping = False

    for key, value in lookup.items():
        if value['start_pos'] <= query['end'] and query['start'] <= value['end_pos']:
            overlap_key = key
            overlap_length = min(query['end'], value['end_pos']) - max(query['start'], value['start_pos'])
            target_length = value['end_pos'] - value['start_pos']
            overlap_proportion = overlap_length / target_length
            is_overlapping = True
            break  # Stop at the first valid overlap

    return {
        "antismash_target": overlap_key,
        "overlap_length": overlap_length,
        "overlap_proportion": overlap_proportion,
        "query_length": query_length,
        "target_length": target_length,
        "is_overlapping": is_overlapping
    }

# Function to create BGC visualizations
def bgc_charts(
        source, 
        var_selection="type", 
        var_y_axis="probability_score", 
        var_x_axis="length",
        var_tooltip=None,
        var_scatter_x_axis_title="BGC Region Length (bp)",
        var_scatter_y_axis_title="GECCO BGC Type Probability (%)",
        var_bar_x_title="BGC Count",
        var_legend_title="GECCO BGC Type",
        var_dropdown_title="GECCO BGC Type"):
    if var_tooltip is None:
        var_tooltip = ['cluster_id', 'sequence_id', 'start', 'end', 'average_p', 'max_p',
                       'type', 'alkaloid_probability', 'nrp_probability',
                       'polyketide_probability', 'ripp_probability', 'saccharide_probability',
                       'terpene_probability', 'proteins', 'domains', 'genome_id', 'length',
                       'probability_score']

    options = [i for i in source[var_selection].unique()]
    resize = alt.selection_interval(bind='scales')
    base = alt.Chart(source)
    input_dropdown = alt.binding_select(options=options + [None],
                                        name=f'{var_dropdown_title} ')
    selection = alt.selection_point(fields=[var_selection], bind=input_dropdown)

    color = alt.condition(
        selection,
        alt.Color(f'{var_selection}:N').legend(None),
        alt.value('lightgray')
    )

    scatter = base.mark_circle(size=75).encode(
        x=alt.X(f'{var_x_axis}:Q', title=var_scatter_x_axis_title),
        y=alt.Y(f'{var_y_axis}:Q', title=var_scatter_y_axis_title).axis(format='%'),
        color=color,
        tooltip=var_tooltip
    ).add_params(selection).properties(height=400, width=600).add_params(resize)

    legend = base.mark_circle(size=75).encode(
        alt.Y(f'{var_selection}:N', title=var_legend_title).axis(orient='right'),
        color=color
    )

    chart2 = base.mark_bar().encode(
        x=alt.X('count()', title=var_bar_x_title),
        y=alt.Y(f'{var_y_axis}:Q', title="").axis(format='%').bin(maxbins=30),
        color=color
    ).add_params(selection).properties(height=400, width=100).add_params(resize)

    return scatter | chart2 | legend

# Load metadata and data files
report_dir = Path("../")
dependency_version_file = report_dir / "metadata/dependency_versions.json"
with open(dependency_version_file, "r") as file:
    dependency_version = json.load(file)
    antismash_version = dependency_version.get("antismash")
    gecco_version = dependency_version.get("gecco")

gecco_table = report_dir / f"gecco/{gecco_version}/gecco_clusters.csv"
source = pd.read_csv(gecco_table)

# Calculate BGC region lengths
source["length"] = source["end"] - source["start"]

# Process probability scores
probability_category = [c for c in source.columns if "probability" in c]
for i in source.index:
    bgc_type = source.loc[i, "type"]
    probability_score = {}
    for t in bgc_type.split(";"):
        for category in [c for c in probability_category if t.lower() in c]:
            probability_score[category] = source.loc[i, category]

    max_value = max(probability_score.values(), default=np.NaN)
    source.loc[i, "probability_score"] = max_value

# Ensure missing probability columns exist
for c in ["alkaloid_probability", "nrp_probability", "polyketide_probability", 
          "ripp_probability", "saccharide_probability", "terpene_probability"]:
    if c not in source.columns:
        source[c] = np.NaN

# Display initial data table
display(HTML(DT(source.loc[:, ["cluster_id", "genome_id", "average_p", "max_p", "type",
                                "length", "probability_score"]], 
                columnDefs=[{"className": "dt-center", "targets": "_all"}], scrollX=True)))

# Process antiSMASH overlap if file exists
antismash_table = report_dir / f"tables/df_regions_antismash_{antismash_version}.csv"
if antismash_table.is_file():
    display(Markdown("## Overlap with antiSMASH BGCs"))
    df_antismash = pd.read_csv(antismash_table).set_index("bgc_id")

    data = []
    for i in source.index:
        accession_id = source.loc[i, "sequence_id"]
        cluster_id = source.loc[i, "cluster_id"]
        subset_antismash = df_antismash[df_antismash["accession"] == accession_id]
        query = source.loc[i, ["start", "end"]].to_dict()
        lookup = subset_antismash.loc[:, ["start_pos", "end_pos"]].T.to_dict()
        overlap = find_overlap(query, lookup)
        overlap["cluster_id"] = cluster_id
        data.append(overlap)

    antismash_overlap = pd.DataFrame(data)
    source2 = source.merge(antismash_overlap, on="cluster_id").drop_duplicates().reset_index(drop=True)

    # Create overlap charts
    chart2 = bgc_charts(source2, 
                        var_selection="is_overlapping", 
                        var_y_axis="probability_score", 
                        var_x_axis="length",
                        var_tooltip=['cluster_id', 'antismash_target', 'overlap_length', 'overlap_proportion', 
                                     'query_length', 'target_length', 'sequence_id', 'start', 'end', 'average_p', 
                                     'max_p', 'type', 'alkaloid_probability', 'nrp_probability',
                                     'polyketide_probability', 'ripp_probability', 'saccharide_probability',
                                     'terpene_probability', 'proteins', 'domains', 'genome_id', 'length',
                                     'probability_score'],
                        var_scatter_x_axis_title="BGC Region Length (bp)",
                        var_scatter_y_axis_title="GECCO BGC Type Probability (%)",
                        var_bar_x_title="BGC Count",
                        var_legend_title="Overlap",
                        var_dropdown_title="Overlap with antiSMASH BGCs?")
    display(chart2)

    # Display the detailed overlap table
    display(HTML(DT(source2.loc[:, ["cluster_id", "genome_id", "average_p", "max_p", "type",
                                    "length", "probability_score", "antismash_target", "overlap_length", 
                                    "overlap_proportion", "is_overlapping"]], 
                    columnDefs=[{"className": "dt-center", "targets": "_all"}], scrollX=True)))

