"""Implementation for converting RNA ss mfold output (.ct files) to dot-bracket format"""

import util, argparse

def get_dot_brackets(ct_file_path, run_id=""):
    all_lines = []
    with open(ct_file_path) as f:
        for line in f:
            line = line.rstrip()
            if line != "": all_lines.append( line.split("\t") )
    
    seq_length = int(all_lines[0][0]) # ct files give the sequence length as first entry in header for each structure
    assert len(all_lines) % (seq_length+1) == 0, ".ct file apparently does not contain discrete number of estimates secondary structures" # +1 to allow for header    

    dot_bracket_results = []
    for i in range( len(all_lines) / (seq_length+1) ):
        subset_lines = all_lines[i*(seq_length+1): i*(seq_length+1)+seq_length+1] # lines from .ct file corresponding to a single structure. NB each structure's representation (including header) will be seq_length+1 lines long
        name = run_id + "." + str(i+1) # not using zero based, to match mfold output file names
        sss = make_dot_bracket(name, subset_lines)
        dot_bracket_results.append( sss )
    return dot_bracket_results

def make_dot_bracket(name, ss_lines):
    """Given the ct data for a single predicted RNA structure, produce a SeqStruct instance"""

    header = ss_lines[0] # eg: "100     dG = -19.24   [initially -21.60] nl43.1.100"
    dG_entry = header[1] 
    
    try:
        dG = float(dG_entry.split(" ")[3])
    except ValueError:
        print "WARNING no dG value found for %s" % name
        dG = float("NaN")
    
    dotbracket = []
    sequence = []
    i = 0
    for row in ss_lines[1:]: # first line is header
        i += 1
        value = int(row[4])
        if value == 0:
            symbol = "." # no base pairing
        elif value > i:
            symbol = "(" # pairing opens
        elif value < i:
            symbol = ")" # pairing closes
        dotbracket.append(symbol)
        sequence.append(row[1])
    return util.SeqStruct(name, "".join(sequence), "".join(dotbracket), dG )



def main():
    parser = argparse.ArgumentParser(description="Inspired by Mike Zucker's Ct2B.pl script, ct2b.py converts RNA secondary structure *.ct files output by mFold to dot-bracket notation")
    parser.add_argument("ct_file", help="File containing one or more RNA secondary structures in *.ct format")
    args = parser.parse_args()
    results = get_dot_brackets(args.ct_file)
    print results[0].seq
    for i in range(len(results)):
        print results[i].ss + "\t(" + str(results[i].dG) + ")"


if __name__ == "__main__":
    main()
