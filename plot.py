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
        self.X = np.sort(list(X)) # NB must use numpy sort, or nan is not handled correctly
        n = len(self.X)
        self.t = max(X)
        self.b = min(X)
        r = float(self.t -self.b) # range
        fbins = float(bins)
        #self.bounds = [ self.b + i * r/fbins for i in range(bins+1)  ] # creates bounds where each bin is of equal width
        self.indices = [ i*int(round(n/fbins)) for i in range(bins) ]
        self.bounds = [ self.X[index] for index in self.indices  ] # creats boundaries where roughly equal number of elements in each bin
        self.bounds.append(self.t)
        
    def get(self, x):
        if np.isnan(x):
            return self.nan_colour
        if x < self.b or x > self.t:
            raise ValueError("Value %f is not within range of min (%f) or max (%f). Bounds: %s"%(x, self.b, self.t, self.bounds))
        for i in range(self.bins):
            if self.bounds[i]-self.tolerance < x < self.bounds[i+1]+self.tolerance:
                return self.colours[i]
        else:
            raise ValueError("ERROR in Colour.get. input value: %f" % x)
def test():
    data = sorted(np.random.rand(2000))
    print "min",min(data)
    print "max",max(data)
    print "data", data
    col = Colour(data, 5)
    print "bounds", col.bounds
    
    freq = dict(zip(col.colours, [0]*len(col.colours)))

    for d in data:
        freq[col.get(d)] += 1
    print freq


def plot_sl_diff_struct():
    """ Want to visualise possible overlap of stem loops"""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("all_sl_tsv", type=str, help="table containing data for all sl, used for determining colours")
    parser.add_argument("struct_specific_sl_tsv", type=str, help="table containing data for only the given structure number (1, 2, 3, or 4, etc)")
    parser.add_argument("struct_number", type=int, help="")
    args = parser.parse_args()
    figsize = (18, 12)
    
    all_data = pd.read_csv(args.all_sl_tsv, sep="\t")
    to_plot = pd.read_csv(args.struct_specific_sl_tsv, sep="\t")
    to_plot.columns = ["index", "run_id", "seq_id", "group", "length", "loop_id", "win_id", "struct_id", "loop_number", "start_a", "end_a", "start_u", "end_u", "struct_number", "s_energy", "dotbracket"]
    struct_number_str = str(args.struct_number)
    
    mask1 = all_data.s_energy > 20
    all_data.loc[mask1, "s_energy"] = np.NaN # https://stackoverflow.com/questions/21608228/conditional-replace-pandas
    mask2 = to_plot.s_energy > 20
    to_plot.loc[mask2, "s_energy"] = np.NaN

    col = Colour(all_data["s_energy"]) # want the colours to be determined using the whole set of energies, not jsut the set associated with this structure
    
    fig, axis = plt.subplots(1, figsize=figsize)
   
    plot_sl(axis, to_plot, col, title="struct_number=%s"%struct_number_str) 

    n_sites = 13890
    xticks = [ x*500 for x in range(n_sites/500 + 1) ] + [n_sites] 
    axis.set_xticks(xticks)
    axis.set_xticklabels(xticks)

    plt.savefig( "%s.struct_number.sl.pdf"%struct_number_str )
   


def plot_sl(axis, df, energy_colouring, title="", shadings=[]): # NB any NaN values must be dealt with in advance
    axis.patch.set_facecolor('lightgray') 
    df = df[ ["seq_id", "group", "start_a", "end_a", "s_energy"] ]
    
    
    df = df.sort_values(["group", "seq_id"]).reset_index(drop=True) # TODO probably want to define order manually

    #df.to_csv("tmp.tsv", sep="\t" )
    groups_list = list( df["group"].unique() )
    seq_list = list( df["seq_id"].unique() )[::-1]

    n_seqs = len(seq_list)

    axis.set_yticks(range(0, n_seqs))
    y_tick_positions = [ n_seqs ] 
    for group in groups_list:
        
        group_df = df[ df["group"] == group ]
        i = 0
        for row in group_df.itertuples():
            if i % 100 == 0:
                print group, i
            i += 1
            
            sl_index, seq_id, group, sl_start, sl_end, energy = row
            
            seq_pos = seq_list.index(seq_id)
            axis.plot([sl_start, sl_end], [seq_pos+0.5, seq_pos+0.5], color=energy_colouring.get(energy))
        y_tick_positions.append(seq_pos)
    
    
    axis.set_yticks(y_tick_positions)
    axis.set_yticklabels(groups_list)
    axis.grid(which="major", axis='y')
    axis.set_ylabel("sequences (by group)")
    axis.set_xlabel("Nucleotide position in alignment")
    axis.set_title(title)
    #axis.set_aspect(0.9)
    return axis

def plot_genome_features(axis, annotations={}, title=""):
    axis.patch.set_facecolor('lightgray') 
    feature_colour = "blue"
    sorted_keys = sorted(annotations.keys())[::-1] # TODO this is a bad quick fix. Sorted order may not be correct order!
    y_pos = 0
    for k in sorted_keys: 
        for feature in annotations[k]:
            axis.plot( feature, [y_pos, y_pos], feature_colour )
        y_pos += 1
    axis.set_yticks(range(len(annotations)))
    axis.set_yticklabels(sorted_keys)
    axis.set_ylim([-0.5, len(sorted_keys)+0.5])
    #axis.set_xticklabels([])
    axis.set_title(title)
    return axis

def draw_shading(axis, shadings, colour="black", transparency=0.75):
    """ Draw rectangles on a genome plot to show, for example, regions bound by NC """
    height = axis.get_ylim()[1]
    for start, end in shadings:
        lower_left = (start, 0)
        width = end - start
        axis.add_patch(Rectangle( lower_left, width, height, alpha=transparency, color=colour ))

def plot_energies(axis, energies, energy_colouring): # TODO could pass in Colour instance and use same as in plot_sl
    axis.patch.set_facecolor('lightgray') 
    col = energy_colouring

    for i in range(len(col.indices)-1): # plot each set of energies separately, with their own colour
        energies = col.X[ col.indices[i] :  col.indices[i+1]] 
        axis.scatter(energies, [1]*len(energies), color=col.colours[i])
    energies = col.X[col.indices[-1]:] 
    axis.scatter(energies, [1]*len(energies), color=col.colours[-1])


    axis.set_ylabel("")
    axis.set_xlabel("Structure energy for window (dG)")
    axis.set_yticks([1])
    axis.set_ylim([1.-0.01, 1.+0.01])
    # x tick labels 1?
    axis.set_yticklabels([])


def plot_sl_alone():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("sl_tsv", help="")
    parser.add_argument("struct_tsv", help="")
    parser.add_argument("seq_tsv", help="")
    parser.add_argument("aligned_ref", help="fasta file containing reference sequence which has been aligned to the dataset you are using")
    parser.add_argument("ref_id", help=" 'hxb2' or 'ben' ")
    parser.add_argument("-t", help="Title for SL figure")
    parser.add_argument("-g", nargs="+", help="Groups to be included in figure")
    parser.add_argument("-n", help="path to fasta containing clip seq nl43 sequence (as first or only sequence)")
    args = parser.parse_args()

    sl_title = "" if args.t == None else args.t
    figsize = (18, 12)

    #plt.figure(figsize=(3,4))
    fig, axes = plt.subplots(2, gridspec_kw = {'height_ratios':[1, 10]}, figsize=figsize, sharex=True)
    
    # making genome feature plot
    aligned_ref = reader.ReadInFasta(args.aligned_ref)[0] # use first sequence, though there may only be one
    ref_indices = util.unalign_indices(aligned_ref[1])
    genome_features = annotations.get_genome_features(args.ref_id) 
    # changing values so they are aligned indices
    for k in genome_features.keys():
        genome_features[k] = [  [ ref_indices.index(pair[0]), ref_indices.index(pair[1]) ]  for pair in genome_features[k]     ]

    plot_genome_features(axes[0], genome_features, title="Genome features (aligned sequence %s)"%aligned_ref[0])
    
    # making SL plot
    sl_df = pd.read_csv(args.sl_tsv, sep="\t", header=0)
    struct_df = pd.read_csv(args.struct_tsv, sep="\t", header=0)
    struct_df = struct_df[ [ "run_id", "struct_id", "s_energy"] ]
    seq_df = pd.read_csv(args.seq_tsv, sep="\t", header=0)
    seq_df = seq_df[ ["run_id", "seq_id", "group"] ]
    
    join_sl_struct = pd.merge(sl_df, struct_df, how="inner", on=["run_id", "struct_id"]).sort_values(["run_id", "struct_id"])
    join_sl_struct_seq = pd.merge(join_sl_struct, seq_df, how="inner", on=["run_id", "seq_id"]).sort_values(["run_id", "struct_id"])
   
    if "all" in args.g:
        to_plot = join_sl_struct_seq
    else:
        to_plot = join_sl_struct_seq[ join_sl_struct_seq["group"].isin(args.g)  ]

   # This is a quick fix, need to address data properly
    # have values with s_energy >> other values, which I assume are junk
    mask = to_plot.s_energy > 20
    to_plot.loc[mask, "s_energy"] = np.NaN # https://stackoverflow.com/questions/21608228/conditional-replace-pandas


    shadings = []
    if args.n != None:
        # need alignment indices for the nl43 sequence used to determine NC binding regions
        clipseq_indices = util.unalign_indices( reader.ReadInFasta(args.n)[0][1] )
        nc_zones = annotations.find_peak_regions(annotations.sites, annotations.mature)
        for s, e in nc_zones: 
            shadings.append( (clipseq_indices.index(s), clipseq_indices.index(e-1) )) # NB making the end inclusive
    
    col = Colour(to_plot["s_energy"])

    t0 = time.time()
    plot_sl(axes[1], to_plot, col, title=sl_title, shadings=shadings) 
    t1 = time.time()
    print "SL plot time: ", t1-t0
    

    draw_shading(axes[0], shadings, transparency=0.25)
    draw_shading(axes[1], shadings, transparency=0.75)

    n_sites = len(ref_indices)
    xticks = [ x*1000 for x in range(n_sites/1000 + 1) ] + [n_sites] 
    for axis in axes:
        axis.set_xticks(xticks)
        axis.set_xticklabels(xticks)
    
    fig2, ax = plt.subplots(1, figsize=(figsize[0], 1))
    plot_energies(ax, to_plot["s_energy"], col)
    

    #plt.savefig( ".".join(args.g)+".pdf" )
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

def energies():
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
    
    print "max",max(data)
    print "med",np.median(data)
    print "min",min(data)

    #print np.percentile(data, [0.0, 20., 40., 60., 80., 100. ])
    #ax.hist(data, bins=5)
    
    # for testing new version for assigning bins in Colour
    col = Colour(data)
    print col.bounds

    freq = dict(zip(col.colours, [0]*len(col.colours)))

    for d in data:
        freq[col.get(d)] += 1
    print freq

    ax.boxplot(data, 0, "gD")
    #ax.scatter([1.0]*len(data), data)
    plt.show()


def multiple_energy_distributions():
    d = "/big_disk/capsid/packaging_signal/analyses/large/merge_orig_restart/run_again/"
    files = [ "1.tmp", "2.tmp", "3.tmp", "4.tmp" ]
    
    f, axis = plt.subplots(1)
    
    for i in range(len(files)):
        df = pd.read_csv(d+files[i], sep="\t")
        df.columns = ["index", "run_id", "seq_id", "group", "length", "loop_id", "win_id", "struct_id", "loop_number", "start_a", "end_a", "start_u", "end_u", "struct_number", "s_energy", "dotbracket"]
        
        data = df["s_energy"]
        axis.scatter(data, [i]*len(data), marker='+')
    plt.show()

if __name__ == "__main__":
    pd.set_option('display.width', 140) # set width before introducing line breaks when printing df
    #plot_all()
    #plot_sl_alone()
    #energies()
    #test()
    #multiple_energy_distributions()
    plot_sl_diff_struct()
