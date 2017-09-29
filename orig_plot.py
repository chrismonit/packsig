import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
import numpy as np
import util
import annotations


class Genome_Plot:
    def __init__(self, length):
        self.length = length
        self.y_pos = 0
        self.x = np.arange(0.0, length, 0.1)
        self.y_values = [ ]
        self.y_labels = []
        
        plt.xlim([0.0, length])
    
    def plot(self, indices, colour=util.DEFAULT_COLOUR):
        assert (isinstance(indices, tuple) or isinstance(indices, list))
        assert len(indices) == 2
        assert -1 < indices[0] < indices[1]
        line = np.arange(indices[0], indices[1], 1.0)
        y = [self.y_pos]*len(line)
        plt.plot(line, y, colour)
    
    def new_line(self, increment, identifier):
        self.y_pos += increment
        self.y_values.append(self.y_pos)
        self.y_labels.append(identifier)

    def finalise_plot(self):
        #plt.ylim([0, max(self.y_values)])
        matplotlib.rc('ytick', labelsize=30) 
        #matplotlib.rcParams.update({'font.size': 5})
        #print self.y_values
        #print self.y_labels
        


        plt.yticks(self.y_values, self.y_labels)
        plt.xlabel("Alignment nucleotide site")
        #plt.grid()
        plt.tight_layout()


    def show(self):
        self.finalise_plot()
        plt.show()

    def save(self, save_path): # needs file extension for format
        self.finalise_plot()
        plt.savefig(save_path)
    
    def plot_hxb2(self, unalign_indices):
        f1 = [ [1, 634], # 5' LTR
            [790, 2292], # gag
            [5041, 5619], # vif
            [8379, 8469], # tat2
            [8797,9417] # nef
        ]# f1
        f2 = [ [5831,6045], # tat1
                [6062,6310], # vpu 
                [8379,8653], #rev2
                [9086,9719] #3' LTR
                ] #f2
        f3 = [ [2085,5096], #pol
                [5559, 5850],# vpr
                [5970,6045], #rev1
                [6225,8795]# env
                ] #f3
        f1 = [ [l[0]-1, l[1]-1] for l in f1 ] # make zero based
        f2 = [ [l[0]-1, l[1]-1] for l in f2 ] # make zero based
        f3 = [ [l[0]-1, l[1]-1] for l in f3 ] # make zero based
        frames = [f1, f2, f3]
        ids = [ "HXB2 F1", "HXB2 F2", "HXB2 F3" ]
        self.new_line(50, "")# make some space
        
        for i in range(len(frames))[::-1]: # go backwards so that the order, downwards on plot, is F1, F2, F3
            #print ids[i]
            self.new_line(70.0, ids[i] )
            for gene in frames[i]:
                align_gene = (unalign_indices.index(gene[0]), unalign_indices.index(gene[1]))
                self.plot(align_gene)
    
    def shading(self, shade_ranges):
        axis = plt.gca()
        for shade in shade_ranges:
            lower_left = (shade[0], 0)
            w = shade[1] - shade[0]
            h = self.y_pos # the (current) highest point
            transparency = 0.5
            axis.add_patch(Rectangle( lower_left, w, h, alpha=transparency ))

def test():

    gp = Genome_Plot(100)

    gp.plot( (2, 7) )
    gp.show()

def main():
    import argparse
    import tdg.utils.readInFastaSequences as reader
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("idxa", help="")
    parser.add_argument("alignment", help="NB HXB2 must be first sequence!!!")
    parser.add_argument("-s", help="file path for saving figure. default is to display immediately")
    args = parser.parse_args()
    
    sequences = reader.ReadInFasta(args.alignment)
    print "assuming sequence %s is HXB2" % sequences[0][0]
    hxb2_unalign_indices = util.unalign_indices(sequences[0][1])
    align_len = len(sequences[0][1])

    f = open(args.idxa)
    data = [line.rstrip().split(util.DELIM) for line in f.readlines()]
    f.close()
    
    for row in data:
        row[1] = row[1].split(".")
        row[2] = float(row[2]) # p frac
        row[3] = int(row[3]) # loop pattern start
        row[4] = int(row[4]) # loop pattern end


    p_fractions = sorted(list(set([ row[2] for row in data ]))) # make a list of the unique p_fracs present

    colours = dict(zip(p_fractions, util.COLOURS)) # colours used for p fractions in plot
    #print "colours", colours

    fig = Genome_Plot(align_len)
    
    WIN_GAP = 2.5
    SEQ_GAP = 5.0

    current_p_frac = None
    current_seq = None # name of sequence
    
    current_group = None
    label = "-------"

    line_num = len(data)
    for row in data[::-1]:
        line_num -= 1
        #print row
        
        increment = 0.0
        new_line = False
        if row[1][0] != current_seq:
            current_seq = row[1][0]
            increment += SEQ_GAP
            new_line = True
        if row[2] != current_p_frac:
            current_p_frac = row[2]
            increment += WIN_GAP
            new_line = True
        
#        if current_p_frac != 0.5: 
#            print "\nWARNING: only showing one set of window results, as a temporary thing\n"
#            continue

        if new_line:
            # silly cludge, to print only names for 0.0 p frac for clarity
            #if current_p_frac == 0.0: identifier = current_seq + "_" + str(current_p_frac)
            #else: identifier = ""
            
            line_group = current_seq.split(".")[0]
            if current_group != line_group:
                identifier = current_seq
                print line_num, row[1]
                current_group = line_group
            else:
                identifier = ""
            fig.new_line(increment, identifier)

        seq_id = row[1][0]
        identifier = seq_id + "_" + str(current_p_frac)
        fig.plot( row[3:5], colours[current_p_frac])
   
    fig.plot_hxb2(hxb2_unalign_indices)
   
    # get rectangle indices, expressed as in unaligned sequence
#    rectangles_unalign = annotations.find_peak_regions(annotations.sites, annotations.mature)
#    
#    # map to aligned sequence positions
#    clip_seq = sequences[-1] 
#    print "Using seqeunce %s as reference for CLIP Seq data" % clip_seq[0]
#    unalign_indices = util.unalign_indices(clip_seq[1])
#    rectangles_align = [  (unalign_indices.index(r[0]), unalign_indices.index(r[1]) ) for r in rectangles_unalign  ]
#    
#    fig.shading( rectangles_align )

    
    if args.s != None:
        fig.save(args.s)
    else:
        fig.show()

if __name__ == "__main__":
    main()
