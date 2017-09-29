NEWLINE = "\n"
DELIM = "\t"
MOTIF = "[Aa]([AaCcGgUuTt]{3,4}){4}[Aa]"
MIN_LOOP_LENGTH = 17
# TODO this is hideous. Want minimum two '(' and ')' within space of 3 characters (allowing for '.' either side of '(...)'
STEM_LOOP_PATTERN = "\(\(\(\.{%d,}\)\)\)|\(\(\.\(\.{%d,}\)\.\)\)|\(\.\(\(\.{%d,}\)\)\.\)" % ((MIN_LOOP_LENGTH,)*3)
PWD="/path/to/working/dir/"
PREF = "XX:" # prefix for printed output lines
MSG = PREF+"MESG" # prefix for printed output lines
RST = PREF+"rnas"
INIT= PREF+"INIT"
INDX = PREF+"IDX"
WARN = PREF+"WARN"
MAX_P_FRACS = 3
COLOURS = [ "black", "red", "green", "blue" ]
DEFAULT_COLOUR = "black"
GAP = "-"
# technically the fasta is an input file I produce for mfold. Missing .ct, .ps and .ss because this is caught with glob, since there can be several
MFOLD_EXTENSIONS = [ ".pnt", ".det", ".out", ".log", ".ann", ".plot", ".sav", ".count", ".h-num", ".ss-count", ".fasta" ]


# extensions you need a wildcard for, since there can be more than one
MFOLD_MULTIPLE_EXTENSIONS = [ ".ct", ".ps", ".ss" ]

class SeqStruct:
    """ RNA sequence and secondary structure """
    def __init__(self, name, seq, ss, dG):
        assert isinstance(name, str)
        assert isinstance(seq, str)
        assert isinstance(ss, str)
        assert isinstance(dG, float)
        self.name = name
        self.seq = seq
        self.ss = ss  # secondary structure
        self.dG = dG
    def __str__(self):
        return NEWLINE.join( [ 
            self.name,
            self.seq,
            self.ss+DELIM+"(%.2f)"%self.dG
            ] )


def readSeqSS(path):
    """ Read in fasta-style file with names, sequences and secondary structures """
    assert isinstance(path, str)
    f = open(path)
    records = f.read().split(">")
    f.close()

    seqStructs = []
    for record in records:
        if record == "": continue
        record = record.split(NEWLINE)  # change record into a list
        record = filter(lambda line: line != "", record)  # filter out any empty lines
        sss = SeqStruct(record[0], record[1], record[2])
        seqStructs.append(sss)
    return seqStructs


def writeSeqSS(data, path):
    """ Write SeqSS data to a file """
    f = open(path, 'w')
    for seqStruct in data:
        f.write(seqStruct + NEWLINE)
    f.close()

def writeSingleSequenceFasta(name, seq, path):
    f = open(path, 'w')
    f.write(">"+name+NEWLINE)
    f.write(seq+NEWLINE)
    f.close()


def msg(run_id, message, prefix=MSG, delim=DELIM):
    print delim.join([ prefix, run_id, message ])

def init_msg(message, prefix=INIT, delim=DELIM):
    print delim.join([ prefix, message ])

def result_msg(seqStruct, prefix=RST, delim=DELIM):
    print prefix+delim+str(seqStruct).replace(NEWLINE, NEWLINE+prefix+delim)


def indx_msg(index_type, loop_id, p_frac, aligned_index_start, aligned_index_end, subseq, prefix=INDX, delim=DELIM):
    #index_type can be "A" or "U" for aligned or unaligned. struct_num is the number of struct given my mfold
    print delim.join( [prefix+index_type, loop_id, str(p_frac), str(aligned_index_start), str(aligned_index_end), subseq] )


def warning_msg(message, prefix=WARN, delim=DELIM):
    print delim.join([prefix, message])


def unalign_indices(aln_seq, gap="-", gap_value=None):
    """ For a sequence of length N interspersed with gap characters, 
        return a list of length N where the ith element is the index
        of the ith site if the sequence contained no gap characters.
        Indices for sites containing gaps in the sequence are
        represented by <gap_value>"""
    indices = []
    gap_count = 0
    for i in range(len(aln_seq)):
        if aln_seq[i] == gap:
            gap_count += 1
            indices.append(gap_value)
        else:
            indices.append(i-gap_count)
    return indices




def test():
    print "test 1"
    aln = "-aaa---t-tt"
    for i in range(len(aln)):
        print unalign_index(aln, i)
    
    print "test 2"
    unalign = unalign_indices(aln)
    print unalign
    
    print unalign.index(5)

if __name__ == "__main__":
    test()
