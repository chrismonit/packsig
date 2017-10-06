import tdg.utils.readInFastaSequences as reader
import pattern, ss_prediction, util
import os, argparse

def check_names_length(names, alignment_length):
    """ mfold seems to allow file lengths no greater than 30 characters, including the extensions it adds, so checking in advance"""
    nameTooLong = False
    alignment_str_len = len(str(alignment_length)) # the length of the string representation of the alignment length (how many digits)
    max_name_len = 30 - ((2*alignment_str_len) + 5) # +2 for dots between values in name and +3 for mfold file extensions
    for name in names:
        if len(name) > max_name_len:
            print "WARNING sequence %s has a name which is too long for mfold. Max name length is %d. Please change" % (name, max_name_len)
            nameTooLong = True
    if nameTooLong: exit()




def check_window_parameters(W, P, min_W=10, min_P=0.0, max_P=0.99):
    """ check the settings used to control window sizes and movements make sense"""
    if W < min_W:
        raise ValueError("The set window width (%d) is less than minimum value allowed (%d)" % (W,min_W))
    if not (min_P <= P < max_P):
        raise ValueError("The set window shift proportion %f does not respect bounds (min %f, max %f)"%(P, min_P, max_P))


def get_windows(L, W, P):
    windows = []

    position = int(round(W * P)) # starting position for first window
    if P > 0.0: 
        windows.append( (0, position) ) # the first window does not start at the first base in the sequence, so we have a special window which is less than W in length
    while position+W < L: # move along windows
       windows.append((position, min(position+W, L))) # lower and upper bounds of the window. using min function prevents us having a window which goes beyond sequence length
       position += W
    if position < L: # if we haven't yet reached the end of the sequence...
        windows.append((position, L)) # ...then add a smaller window including what's left
    return windows 
    
def windows_sanity_check(windows, sequence):
    """ sanity check. Checking the windows add up to the whole sequence """
    testSeq = ""
    for window in windows:
        testSeq += sequence[window[0]:window[1]]
    assert testSeq == sequence, "Sanity check for sliding windows has failed"




def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("alignment", help="Fasta format alignment of sequences to be analysed")
    parser.add_argument("width", help="Width of nucleotide window to analyse")
    parser.add_argument("-p", help="Overlap proportions for windows", nargs="+", default=[])
    parser.add_argument("-wd", help="Working directory", default=None)
    args = parser.parse_args()

    if args.wd == None:
        WD = os.getcwd() + "/"
    else:
        raise ValueError("User specified working directory not implemented!")
        #if os.path.exists(args.wd): WD = args.wd
        #else: raise ValueError("Working directory path does not exist")
    W = int(args.width)
    
    try: user_p_fractions = [ float(x) for x in args.p]
    except ValueError:
        print "Error: Overlap proportion arguments must be numerical values, separated by space character"
        exit()
     
    if len(user_p_fractions) > util.MAX_P_FRACS: raise ValueError("Error: Maximum of %d user-defined window shifts" % util.MAX_P_FRACS)
    

    p_fractions = [ 0.0 ] + user_p_fractions
    
    colours = dict(zip(p_fractions, util.COLOURS)) # colours used for p fractions in plot

    sequences = reader.ReadInFasta(args.alignment)
    
    util.init_msg("motif_pattern='%s'"%util.MOTIF)
    util.init_msg("stem_loop_pattern='%s'"%util.STEM_LOOP_PATTERN)
   
    # NB for current setup to work, the ref seq needs to be the longest seq.
    # sites not found in the ref seq can't be described for others
    reference_seq = sequences[0][1] # expecting this will include gaps
    
    ref_unalign_indices = util.unalign_indices(reference_seq) # indices for reference sequence unaligned
    ref_unalign_len = max(ref_unalign_indices) # assuming first seq in file is the reference
    
    util.init_msg("width=%d" % W)

    align_len = len(reference_seq) # assuming aligned properly and all the same length
    check_names_length([pair[0] for pair in sequences], align_len) # most unaligned sequences will have length less than align_len, but this is the max length possible
    util.init_msg("aln_length=%d"%align_len)
    
    SEQ_NAME_DELIM = "."
    
    HEAD = "HEADER_"
    TAG_SEQ = "SEQ:"
    TAG_WIN = "WIN:"
    TAG_STRUCT = "STRUCT:"
    TAG_LOOP = "SL:"
    print util.tabmsg([HEAD+TAG_SEQ, "seq_id", "group", "length"])
    print util.tabmsg([HEAD+TAG_WIN, "win_id", "seq_id", "shift", "start_a", "end_a", "start_u", "end_u", "has_motif"])
    print util.tabmsg([HEAD+TAG_STRUCT, "struct_id", "seq_id", "win_id", "struct_number", "energy", "dotbracket"])
    print util.tabmsg([HEAD+TAG_LOOP, "loop_id", "seq_id", "win_id", "struct_id", "loop_number", "start_a", "end_a", "start_u", "end_u"])
    
    win_id, struct_id, loop_id = (0, 0, 0)

    # do analysis
    for nameSeqPair in sequences:
    
        genome_with_gaps = nameSeqPair[1].upper().replace("T", "U")
        unalign_indices = util.unalign_indices(genome_with_gaps)
        
        genome = genome_with_gaps.replace(util.GAP, "")

        group, seq_id = nameSeqPair[0].split(SEQ_NAME_DELIM) # assuming name is in format GROUP.SEQ_ID
        
        L = len(genome)

        #util.init_msg("sequence=%s, unalign_length=%d" % (seq_id, L))
        print util.tabmsg([TAG_SEQ, seq_id, group, L])

        for P in p_fractions:
            
            util.init_msg("shift=%.2f" % P)
            
            check_window_parameters(W, P)
            windows = get_windows(L, W, P)

            for window in windows: 
                util.init_msg("windows:%d.%d" % (window[0]+1, window[1])) # correct for zero based. Don't +1 for end so numbers are inclusive

            windows_sanity_check(windows, genome)

            for window in windows:
                win_start_u, win_end_u = window # unaligned numbers. NB indices are like slice numbers (upper is exclusive)
                win_start_a, win_end_a = unalign_indices.index(win_start_u), unalign_indices.index(win_end_u-1)+1  # aligned numbers
                
                run_id = ".".join( [ seq_id, str(win_start_u+1), str(win_end_u) ] )  # +1 to correct for zero based, and don't -1 for win_end_u for the same reason
                win_seq = genome[win_start_u:win_end_u]
                
                has_motif = pattern.has_motif(win_seq)
                print util.tabmsg([TAG_WIN, win_id, seq_id, P, win_start_a, win_end_a, win_start_u, win_end_u, has_motif ])
                
                #util.msg(run_id, "seq:"+win_seq)
                if has_motif:
                    util.msg(run_id, "contains motif, estimating SS")
                    
                    seq_structs = ss_prediction.predict_ss(win_seq, run_id, WD) # will return empty list if no .ct file available
                    
                    util.msg(run_id, "%d SS produced"%len(seq_structs))
                    for iSeqStruct in range(len(seq_structs)):
                        print util.tabmsg([TAG_STRUCT, struct_id, seq_id, win_id, iSeqStruct+1, seq_structs[iSeqStruct].dG, seq_structs[iSeqStruct].ss ])

                        motif_loops = pattern.find_motif_loops(seq_structs[iSeqStruct]) # list of subsequence indices for stem loops containing pattern in this ss
                        if len(motif_loops) == 0:
                            util.msg(run_id, "SS_%d does NOT have loop_motif"%(iSeqStruct+1))
                        else:
                            util.msg(run_id, "SS_%d DOES have %d motif_loop(s):"%(iSeqStruct+1, len(motif_loops)))
                            util.result_msg(seq_structs[iSeqStruct]) # print the whole SeqStruct
                            for iLoop in range(len(motif_loops)): # NB motif_loop indices are like slices: start is inclusive, end is exclusive
                                loop_start_u, loop_end_u = motif_loops[iLoop][0]+win_start_u, motif_loops[iLoop][1]+win_start_u # _u meaning unaligned
                                loop_start_a, loop_end_a  = unalign_indices.index(loop_start_u), unalign_indices.index(loop_end_u)
                                
                                print util.tabmsg([TAG_LOOP, loop_id, seq_id, win_id, struct_id, iLoop+1, loop_start_a, loop_end_a, loop_start_u, loop_end_u])

                                #loop_id = ".".join([run_id, str(iSeqStruct+1), str(iLoop+1)])
                                #stem_loop_seq = win_seq[motif_loops[iLoop][0]:motif_loops[iLoop][1]]
                                # for loop_start_u/loop_end_u, +1 correct for zero based for start. Not +1 for end because index value is exclusive
                                #util.indx_msg("U", loop_id, P, (loop_start_u+1), loop_end_u, stem_loop_seq) # this seq's indices
                                #util.indx_msg("A", loop_id, P, (loop_start_a+1), loop_end_a, stem_loop_seq) # the ref seq's indices
                                loop_id += 1
                        struct_id += 1
                else: # if has motif
                    util.msg(run_id, "no motif in this window")
                win_id += 1
            
    print "End of analysis"

if __name__ == "__main__":
    main()
