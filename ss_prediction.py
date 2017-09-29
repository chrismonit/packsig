"""given an RNA sequence, provide a predicted secondary structure"""
import subprocess, os
import util
import ct2b
from glob import glob # "The glob module finds all the pathnames matching a specified pattern according to the rules used by the Unix shell, although results are returned in arbitrary order."

def predict_ss(seq, run_id, working_directory):

    run_mfold(seq, run_id, working_directory)
    ct_file_path = working_directory+run_id+".ct"
    
    if os.path.exists(ct_file_path):
        seqStructs = ct2b.get_dot_brackets(ct_file_path, run_id=run_id)
    else:
        print util.warning_msg(run_id+util.DELIM+".ct file not found (mfold may not have run analysis): "+ct_file_path)
        seqStructs = []

    rm_mfold_files(run_id, working_directory)
    return seqStructs

def run_mfold(seq, run_id, working_directory):

    fastaPath = working_directory + run_id + ".fasta"
    util.writeSingleSequenceFasta(run_id, seq, fastaPath)
    #print ">>>>> mFold start: %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" % run_id
    FNULL = open(os.devnull, 'w')
    mFoldProcess = subprocess.Popen(["mfold", "SEQ="+fastaPath], cwd=working_directory, stdout=FNULL, stderr=subprocess.STDOUT
) # NB using default inputs here
    mFoldProcess.wait() # wait until this child process has finished
    FNULL.close()
    #print "<<<<< mFold end: %s <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" % run_id



def rm_mfold_files(run_id, working_directory):
    # choosing not to use glob to delete all files, in case a file I haven't yet seen is created and associated with an error

    for extension in util.MFOLD_EXTENSIONS:
        path = working_directory+run_id+extension
        try:
            os.remove(path)
        except OSError:
            print "A Cannot remove file %s" % path
            # do nothing else

    # remove *.ct, *.ps and *.ss files separately, because there may be several with <run_id>_N.ct
    for extension in util.MFOLD_MULTIPLE_EXTENSIONS:
        g = glob(working_directory + run_id + "*" + extension)
        for path in g:
            try:
                os.remove(path)
            except OSError:
                    print "B Cannot remove ct file(s) for run_id %s" % run_id 
                    # do nothing else

def test():
    seq = "ACGAUCGAUUUAUCGGAAAAAUUUUUCCGAUAUAGC"
    run_id = "test_sequence"
    wd = "/big_disk/capsid/packaging_signal/ca-ps_pipleine/pipeline_test/"
    run_mfold(seq, run_id, wd, wd+run_id+".fasta")

    ss = ct2b.get_dot_bracket(wd+run_id+".ct")
    print ss[0]
    print ss[1]

if __name__ == '__main__':
    test()
