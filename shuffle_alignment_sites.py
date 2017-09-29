import random, argparse
import tdg.utils.readInFastaSequences as reader

DELIM = "\t"

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("alignment", help="Fasta format alignment of sequences to be shuffled")
    args = parser.parse_args()

    sequences = reader.ReadInFasta(args.alignment)
    
    #check all sequences are the same length
    for i, j in zip(range(len(sequences)), range(len(sequences))):
        if i >= j: continue
        seq1 = sequences[i][1]
        seq2 = sequences[j][1]
        if len(seq1) != len(seq2):
            raise ValueError("All sequences in alignment must be the same length")

    alignment_len = len(sequences[0][1]) # length of first sequence

    orig_indices = range(alignment_len)
    shuffled_indices = range(alignment_len) # not yet shuffled
    random.shuffle(shuffled_indices)

    header = DELIM.join(["i", "f(i)"])
    print header
    for i in range(len(orig_indices)):
        print DELIM.join([ str(orig_indices[i]), str(shuffled_indices[i]) ])
    
    shuffled = []
    for pair in sequences:
        name = pair[0]
        seq = pair[1]

        shuffled_seq = "".join([ seq[index] for index in shuffled_indices ])
        shuffled.append([name, shuffled_seq])

    reader.WriteFasta(shuffled, "shuffled.fasta")

if __name__ == "__main__":
    main()
