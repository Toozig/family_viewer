from Bio import  Entrez, SeqIO
# from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd
import seaborn as sns
from dna_features_viewer import GraphicFeature, GraphicRecord
import argparse
from pygenometracks.tracksClass import PlotTracks
from pygenometracks.plotTracks import parse_arguments
#from pygenometracks._version import __version__
from pygenometracks.utilities import InputError, get_region
import matplotlib.pyplot as plt
import numpy as np


MAX_SCORE_TFBS = 'maximal score TFBS per site'
SHOW_SEQ = 'Show sequence'
DEFAULT_FIGURE_WIDTH = 40  # in centimeters



def get_seq(enh_dict, email):
    Entrez.email = email
    chrom_dict = {'X':'23', 'Y':'24'}
    cur_chrom = enh_dict['CHROM'].replace('chr','').upper()
    cur_chrom = chrom_dict[cur_chrom] if cur_chrom in chrom_dict else  cur_chrom
    chrom_id  = "NC_000001"[:-len(cur_chrom)] + cur_chrom
    # return
    handle = Entrez.efetch(db="nucleotide", 
                           id=chrom_id, 
                           rettype="fasta", 
                           strand=1, 
                           seq_start=enh_dict['from'], 
                           seq_stop=enh_dict['to'])
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record.seq


#### plot function

def get_variant_features(variants: pd.DataFrame ):
    yellow_pastel = "#FFFFCC"

    variant_features = variants.apply(lambda x: GraphicFeature(start=x.POS, end=x.POS,
                                                                strand=+1,
                                                                color=yellow_pastel,
                                                                label=x.label, axis=1),
                                      axis=1)
    return variant_features


def plot_get_bs_feature(seq_df, bs_col, cmap='pastel'):
    """
    mooma
    """
    colors = sns.color_palette(cmap)


    bs_list = seq_df[bs_col].unique()
    bs_list = bs_list[bs_list != '']

    return seq_df.apply(lambda x: GraphicFeature(start=x.start, end=x.end,
                                       strand=x.strand,
                                        color=colors[0],  #colors_dict[x[bs_col]],
                                         label=x[bs_col],
                                        ), axis=1)


def make_track_plot(config_file: str,
              enh_dict: dict):
    peak_str = f'{enh_dict["CHROM"]}:{enh_dict["from"]}-{enh_dict["to"]}'
    ARGS = [
    "--tracks", config_file,
    "--region", peak_str,
    "--trackLabelFraction", "0.2",
    "--dpi", "130",
    "--trackLabelHAlign", "center",
    "-o", "tmp.png"
    ]

    args = parse_arguments().parse_args(ARGS)
        # Identify the regions to plot:
    if args.BED:
        regions = []
        for line in args.BED.readlines():
            try:
                chrom, start, end = line.strip().split('\t')[0:3]
            except ValueError:
                continue
            try:
                start, end = map(int, [start, end])
            except ValueError as detail:
                    raise ValueError(f"Invalid value found at line\t{line}\t. {detail}\n")
            regions.append((chrom, start, end))
    else:
        regions = [get_region(args.region)]

    if len(regions) == 0:
        raise InputError("There is no valid regions to plot.")

    # Create all the tracks
    trp = PlotTracks(args.tracks.name, args.width, fig_height=args.height,
                     fontsize=args.fontSize, dpi=args.dpi,
                     track_label_width=args.trackLabelFraction,
                     plot_regions=regions, plot_width=args.plotWidth)

    current_fig = trp.plot(args.outFileName, *regions[0], title=args.title,
                               h_align_titles=args.trackLabelHAlign,
                               decreasing_x_axis=args.decreasingXAxis)
    
    return current_fig





def make_TFBS_plot(tfbs_df: pd.DataFrame ,
              var_df: pd.DataFrame,
              enh_dict: dict,
              conservation: list,
             n_lines: int,
             show_seq: bool):
    

    tfbs_df = tfbs_df.copy()
    var_df = var_df.copy()

    features = plot_get_bs_feature(tfbs_df,'name').tolist()

    labels = '>>>>' + (pd.Series(range(len(var_df))) + 1).astype('str') + '<<<<'
    var_df.insert(0,'label',labels.tolist())

    var_df = var_df.astype({'label':str})
    features += get_variant_features(var_df).tolist()
    seq_len = enh_dict['to'] - enh_dict['from']
    

    enh13 = GraphicRecord(sequence_length=seq_len,
                              sequence = None if not show_seq else get_seq(enh_dict, 'manoneh322@konican.com') ,
                          features=features,
                          first_index=enh_dict['from'] ) 
    # plot the enhancer
    all_enh = enh13.plot_on_multiple_lines_with_density(n_lines=round(n_lines),
                                                        nucl_per_line=2000,
                                                        plot_sequence=show_seq,
                                                        figure_width=15,
                                                        density_list = conservation,
                                                        density_label='conservation')
    return all_enh

def get_TFBS_plot(app_stats,
                n_lines,
                checked_box,
                peak_id = None,
                source_name= None,
                score_threshold = None,
                family_id = None):
    # print(f'setting plot with {peak_id['INTERVAL_ID']}, {source_name}, {score_threshold}, {family_id}, {n_lines}')
    one_per_site = MAX_SCORE_TFBS in checked_box
    tfbs_df = app_stats.get_tfbs_filtered_df(one_per_site)
    var_df = app_stats.get_variant_df()
    peak_dict = app_stats.get_peak_data().to_dict('records')[0]
    conservation_list = app_stats.get_peak_conservation(peak_dict)
    show_seq = SHOW_SEQ in checked_box
    return make_TFBS_plot(tfbs_df, var_df, peak_dict, conservation_list, n_lines, show_seq)


def get_track_plot(app_stats,
                n_lines,
                checked_box,
                peak_id = None,
                source_name= None,
                score_threshold = None,
                family_id = None):
    print(f'setting plot with {peak_id}, {source_name}, {score_threshold}, {family_id}, {n_lines}')
    one_per_site = MAX_SCORE_TFBS in checked_box
    tfbs_df = app_stats.get_tfbs_filtered_df(one_per_site)
    var_df = app_stats.get_variant_df()
    peak_dict = app_stats.get_peak_data().to_dict('records')[0]
    conservation_list = app_stats.get_peak_conservation(peak_dict)
    show_seq = SHOW_SEQ in checked_box
    return make_track_plot(tfbs_df, var_df, peak_dict, conservation_list, n_lines, show_seq)



def create_legend():
    fig = plt.figure(figsize=(3, 2))  # Adjust the size as needed
 
    # Create an axis (not displayed, only for the legend)
    ax = fig.add_subplot(111)
 
    # Define the lines with different colors and line styles
    from_mom = plt.Line2D([], [], color='magenta', linestyle='--', label='XX Inheritance')
    from_dad = plt.Line2D([], [], color='blue', linestyle='--', label='XY Inheritance')
    from_both = plt.Line2D([], [], color='purple', linestyle='--', label='Both Inheritance')
    no_data = plt.Line2D([], [], color='gray', linestyle='-', label='No data')
 
    # Add the lines to the legend
    ax.legend(handles=[from_mom, from_dad,from_both, no_data], loc='center')
 
    # Remove the axis to display only the legend
    ax.axis('off')
    return fig