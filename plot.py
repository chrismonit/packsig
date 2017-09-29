import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
import numpy as np
import util
import pandas as pd
import tdg.utils.readInFastaSequences as reader
import argparse

def plot_sl(axis, df):
#    parser = argparse.ArgumentParser(description="")
#    parser.add_argument("sl_tsv", help="")
#    args = parser.parse_args()
#   
    #df = pd.read_csv(args.sl_tsv, sep="\t", header=0)
    #n_sites = int(args.length)
    
    #f, ax = plt.subplots(1, figsize=(24,12))
    
    seq_ids = df["seq_id"].unique()
    n_seqs = seq_ids.shape[0]
    height = n_seqs

    for seq_id in seq_ids:
        for row in df[df["seq_id"] == seq_id][["sl_start", "sl_end", "shift"]].itertuples():
            sl_index, sl_start, sl_end, shift = row
            
            sl = [ sl_start, sl_end ]
            y = [height, height]
            
            line_colour = "k" if shift == 0.0 else "r"
            axis.plot(sl, y, color=line_colour)

        height -= 1
#    plt.show() 
    return axis

def plot():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("count_features_tsv", help="")
    parser.add_argument("fasta", help="")
    parser.add_argument("sl_tsv", help="")
    args = parser.parse_args()
    
    df = pd.read_csv(args.count_features_tsv, sep="\t", header=0)
    sequences = reader.ReadInFasta(args.fasta)
    
#    sequences = make_aln_name_seq_id(input_sequences)
    n_sites = len(sequences[0][1])
    n_seqs = len(sequences)

    df = pd.read_csv(args.count_features_tsv, sep="\t")
    
    if n_sites != df.shape[0]:
        raise ValueError("Number of rows in df (%d) does not match number of aln sites (%d)" % (df.shape[0], n_sites))

    f, axarr = plt.subplots(4)
    
    xticks = [ x*1000 for x in range(n_sites/1000) ]
    
    for ax in axarr:
        ax.grid()
        ax.set_xticks(xticks)

    axarr[0].set_title('Non-gap counts')
    axarr[0].set_ylabel("Sequences without alignment gap (total %d)" % n_seqs)
    axarr[0].set_ylim([0, n_seqs])
    axarr[0].bar(range(1, n_sites+1), df["non_gap"])

    axarr[1].set_title('Proportion of non-gap sequences containing SL (non_gap_sl/non_gap)')
    axarr[2].set_ylabel("Proportion")
    axarr[1].set_ylim([0, 1.0])
    axarr[1].bar(range(1, n_sites+1), df["(non_gap_sl/non_gap)"])
   
    axarr[2].set_title('Proportion of SL-containing sequences (sl_total/n_seq)')
    axarr[1].set_ylim([0, 1.0])
    axarr[2].bar(range(1, n_sites+1), df['sl_total/n_seq'])
   
    
    sl_df = pd.read_csv(args.sl_tsv, sep="\t", header=0)
    axarr[3] = plot_sl(axarr[3], sl_df)

    plt.show()


if __name__ == "__main__":
    plot()
