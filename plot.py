import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
import numpy as np
import util
import annotations
import pandas as pd
import tdg.utils.readInFastaSequences as reader
import argparse
import time


class Colour:
    """ This class takes a series of values X and determines appropriate bins, as one would for a histogram
        It assigns colours to each of these bins, so for a single value x, can provide the appropriate colour label """
    def __init__(self, X, bins=5, tolerance=1e-12):
        self.colours = [ 'blue', 'green', 'yellow', 'orange', 'red'  ]
        self.nan_colour = "black"
        if bins != len(self.colours):
            raise ValueException("Number of colour bins not implemented! Maximum of %d (%d given)"%(len(colours), bins))
        self.bins = bins

        self.tolerance = tolerance
        self.X = X
        self.t = max(X)
        self.b = min(X)
        r = float(self.t -self.b) # range
        fbins = float(bins)
        self.bounds = [ self.b + i * r/fbins for i in range(bins+1)  ]
        #print "min: %f, max: %f, range: %f, bins: %d: range/bins: %f" % (self.b, self.t, r, self.bins, r/fbins)
        #print "bounds:", self.bounds
    
    def get(self, x):
        if np.isnan(x):
            return self.nan_colour
        if x < self.b or x > self.t:
            raise ValueError("Value %f is not within range of min (%f) or max (%f)"%(x, self.b, self.t))
        for i in range(self.bins):
            if self.bounds[i]-self.tolerance < x < self.bounds[i+1]+self.tolerance:
                return self.colours[i]
        else:
            raise ValueError("ERROR in Colour.get")
def test():
    data = np.arange(-10, 11, 1.)
    print data
    col = Colour(data, 5)
    
    print col.get(np.NaN)


def plot_sl(axis, df, annotations={}):
#    parser = argparse.ArgumentParser(description="")
#    parser.add_argument("sl_tsv", help="")
#    args = parser.parse_args()
#   
    #df = pd.read_csv(args.sl_tsv, sep="\t", header=0)
    #n_sites = int(args.length)
    
    #f, ax = plt.subplots(1, figsize=(24,12))
    
    #df = df[ df["group"] == "C" ] # for testing with smaller groups
    
    #df = df[ ["seq_id", "sl_start", "sl_end", "shift", "s_energy"] ]
    df = df[ ["seq_id", "group", "start_a", "end_a", "s_energy"] ]
    
    seq_ids = df["seq_id"].unique()
    n_seqs = seq_ids.shape[0]
    height = n_seqs + len(annotations)

    feature_colour = "black"
    for k in sorted(annotations.keys()): # TODO this is a bad quick fix. Sorted order may not be correct order!
        for feature in annotations[k]: # feature is probably a gene
            axis.plot( feature, [height, height], feature_colour )
        height -= 1

    # This is a quick fix, need to address data properly
    # have 5 values with s_energy >> other values, which I assume are junk
    # ie ~10^6, while next highest is 7.56
    #df.ix[df.s_energy > 20., 's_energy'] = np.NaN
    mask = df.s_energy > 20
    df.loc[mask, "s_energy"] = np.NaN # https://stackoverflow.com/questions/21608228/conditional-replace-pandas
    
    col = Colour(df["s_energy"])
    
    i_seq_id = 0
    for seq_id in seq_ids:
        if i_seq_id % 100 == 0: print i_seq_id
        i_seq_id += 1
        
        y = [height, height]
        #drop = []
        for row in df[df["seq_id"] == seq_id].itertuples(): # TODO inefficient
            sl_index, seq_id, group, sl_start, sl_end, energy = row
            #drop.append(sl_index) 
            sl = [ sl_start, sl_end ]
            
            #line_colour = "k" if shift == 0.0 else "r"
            line_colour = col.get(energy)
            axis.plot(sl, y, color=line_colour)
        
        #df.drop(df.index[ drop ], inplace=True) # lets the df get smaller as we move through it
        height -= 1
#    plt.show() 
    return axis

def plot_sl_alone():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("sl_tsv", help="")
    parser.add_argument("struct_tsv", help="")
    parser.add_argument("seq_tsv", help="")
    parser.add_argument("aligned_ref", help="fasta file containing reference sequence which has been aligned to the dataset you are using")
    args = parser.parse_args()
    
    aligned_ref = reader.ReadInFasta(args.aligned_ref)[0] # use first sequence, though there may only be one
    print "Using sequence %s as aligned reference sequence" % aligned_ref[0]
    
    aligned_indices = util.unalign_indices(aligned_ref[1])
    ben = annotations.make_ben() 
    # changin values so they are aligned indices
    for k in ben.keys():
        ben[k] = [  [ aligned_indices.index(pair[0]), aligned_indices.index(pair[1]) ]  for pair in ben[k]     ]


    sl_df = pd.read_csv(args.sl_tsv, sep="\t", header=0)
    struct_df = pd.read_csv(args.struct_tsv, sep="\t", header=0)
    struct_df = struct_df[ [ "run_id", "struct_id", "s_energy"] ]
    seq_df = pd.read_csv(args.seq_tsv, sep="\t", header=0)
    seq_df = seq_df[ ["run_id", "seq_id", "group"] ]
    
    # original that doesnt work:
    #join = pd.merge(sl_df, struct_df, how="inner", on=['struct_id'])
    
    join_sl_struct = pd.merge(sl_df, struct_df, how="inner", on=["run_id", "struct_id"]).sort_values(["run_id", "struct_id"])
    
    join_sl_struct_seq = pd.merge(join_sl_struct, seq_df, how="inner", on=["run_id", "seq_id"]).sort_values(["run_id", "struct_id"])
    
    f, axis = plt.subplots(1)
    
    t0 = time.time()
    plot_sl(axis, join_sl_struct_seq, ben)
    t1 = time.time()
    print "time: ", t1-t0
    plt.show()

def plot_all():
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

def energies_hist():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("energies_tsv", help="a tsv containing a column with name 's_energy' ")
    args = parser.parse_args()
    
    df = pd.read_csv(args.energies_tsv, sep="\t", header=0)

    f, ax = plt.subplots(1)
#    mu, sigma = 0, 1
#    x = mu + sigma*np.random.randn(500)
#    data=x
    data_orig = df["s_energy"].dropna().values

    C = 1000 # some of the energy values are huge (and therefore, presumably, nonsense)
    data = [ x for x in data_orig if x < C]
   
    ax.hist(data, bins=5)
    plt.show()


if __name__ == "__main__":
    pd.set_option('display.width', 140) # set width before introducing line breaks when printing df
    #plot_all()
    plot_sl_alone()
    #energies()
    #test()
