from Bio import  Entrez, SeqIO
# from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd
import seaborn as sns

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


def make_plot(config_file: str,
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



def get_peak_plot(app_stats,
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
    return make_plot(tfbs_df, var_df, peak_dict, conservation_list, n_lines, show_seq)
