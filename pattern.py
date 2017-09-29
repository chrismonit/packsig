import util
import re

def has_motif(seq, motifPattern=util.MOTIF):
    """ basic check to see if a string sequence contains the pentamer motif"""
    assert isinstance(seq, str)
    return ( re.search(motifPattern, seq) != None)

def find_motif_loops(seq_struct, loop_pattern=util.STEM_LOOP_PATTERN, motifPattern=util.MOTIF):
    """ Returns a list of tuples containing indices for subsequences that contain a motif-containing loop.
        Returns an empty list if there are none. """ 
    
    loop_motif_spans = []
    for loop_match in re.finditer(loop_pattern, seq_struct.ss):
        loop_seq = seq_struct.seq[loop_match.start():loop_match.end()] # subsequence corresponding to loop
        if has_motif(loop_seq):
            loop_motif_spans.append( loop_match.span() )

    return loop_motif_spans # will be empty if there are no matches
       
def test():
    print "TESTING ONLY"
    rna = "...(((.....................)))...(((...................)))...."
    seq = "A"*len(rna)
    result = find_motif_loops(util.SeqStruct("name", seq, rna, 0.3))
    print result

if __name__ == "__main__":
    test()


