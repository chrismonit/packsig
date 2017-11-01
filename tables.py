""" Methods for manipulating tables with pandas, mostly one-off functions"""


import pandas as pd
import numpy as np


def prep():
    # merge orig results and restart jobs
    # assign categories
    # sort them by sequence name, or category
    # give unique indexers

    wd = "/big_disk/capsid/packaging_signal/ca-ps_pipleine/pipeline_test/large/merge_orig_restart/"
    f1 = "orig.idxa.tsv"
    f2 = "corrected_restart.idxa.tsv"

    df1 = pd.read_csv(wd+f1,sep='\t', header=0)
    df2 = pd.read_csv(wd+f2,sep='\t', header=0)


    df3 = pd.concat([df1, df2], ignore_index=True) # ignore index means new df has a new set of row numbers


    df3.drop('id', axis=1, inplace=True) # drop column

    # want to add a "groups" column 
    # add a second groups column, with values tbc
    # eg subtype, and useful classification
    groups = []
    for i in range(df3.shape[0]):
        name = df3.loc[i, "name"].replace("-", "_")
        elements = name.split("_")
        groups.append(elements[0])
        new_name = "_".join(elements[1:])
    #    if len(elements) > 2:
    #        print name, elements, new_name
        # reassign name
        df3.set_value(i, "name", new_name)

    df3.insert(0, "group", groups)
    df3.rename(index=str, columns={"name":"seq_id"}, inplace=True)
    df3['loop'] = df3['loop'].apply(lambda x: int(x)) # original values are floats, dont know why

    # found there are some duplicate entries for "700010607_tfIMC", which are a subset of entries for "700010607"
    df3.drop( df3[df3["seq_id"] == "700010607_tfIMC"].index, inplace=True)
    
    # currently have two HXB2 sequences, "HXB2" and "ref.K03445"
    df3.drop( df3[df3["seq_id"] == "HXB2"].index, inplace=True)
    
    df3 = df3.sort_values(["group", "seq_id", "w_start", "struct", "loop", "shift"]).reset_index(drop=True)

    df3.to_csv("large_dataset_out.tsv", sep="\t", index_label="index")

def get_names():
    """ Takes the data and produces the sequence names in form <group>.<seq_id>
        For comparing with sequence names in alignment file
    """
    wd = "/big_disk/capsid/packaging_signal/ca-ps_pipleine/pipeline_test/large/merge_orig_restart/"
    f = "large_dataset_out.tsv"
    df = pd.read_csv(wd+f, sep="\t", header=0)
    
    df.drop_duplicates(["seq_id"], inplace=True)

    to_drop = set(df.columns.values) - set(["group", "seq_id"])
    df.drop( to_drop, axis=1, inplace=True )

    for i in range(df.shape[0]):
        print str(df.iloc[i].loc["group"]) + "." + df.iloc[i].loc["seq_id"]
    
def prep_dotbracket():
    wd = "/big_disk/capsid/packaging_signal/ca-ps_pipleine/pipeline_test/large/merge_orig_restart/"
    orig = "orig_structures.tsv"
    restart = "restart_structures.tsv"
    orig_df = pd.read_csv(wd+orig, sep="\t", header=0)
    restart_df = pd.read_csv(wd+restart, sep="\t", header=0)
    
    cat_df = pd.concat([orig_df, restart_df], ignore_index=True) # ignore index means new df has a new set of row numbers
    
    cat_df.drop_duplicates(inplace=True)
    
    # sort by group
    cat_df = cat_df.sort_values(["group", "seq_id", "w_start", "struct"]).reset_index(drop=True)
    
    new_order = ["group", "seq_id", "w_start", "w_end", "struct", "energy", "w_seq", "dotbracket"]
    cat_df = cat_df[ new_order ]

    # NB I later renamed "energy" as "w_energy" manually in the file
    cat_df.to_csv("structures.tsv", sep="\t", index_label="index", na_rep="NA")

def join_large_structures():
    """ So far have made two tables, one with summaries of SLs and one with RNA structures and energies
    Now want to join these together """

    wd = "/big_disk/capsid/packaging_signal/ca-ps_pipleine/pipeline_test/large/merge_orig_restart/"
    large = "large_dataset_out.tsv"
    structures = "structures.tsv"
    df1 = pd.read_csv(wd+large, sep="\t", header=0)
    df2 = pd.read_csv(wd+structures, sep="\t", header=0)

    df2.drop( ["w_seq", "dotbracket"], axis=1, inplace=True )

    join = pd.merge(df1, df2, how="left", on=["group", "seq_id", "w_start", "w_end", "struct"])

    join.drop( ["index_x", "index_y"], axis=1, inplace=True )
    join.to_csv("join.tsv", sep="\t", index_label="index")

def test_join1():
    left = pd.DataFrame(
            {'key': ['K0', 'K1', 'K2', 'K3'],
               'A': ['A0', 'A1', 'A2', 'A3'],
               'B': ['B0', 'B1', 'B2', 'B3']}
            )

    right = pd.DataFrame(
            {'key': ['K1', 'K3'],
               'C': ['C1', 'C3'],
               'D': ['D1', 'D3']}
            )
    result = pd.merge(left, right, how="inner", on='key') 
    print "L"
    print left
    print "R"
    print right
    print "join"
    print result


def test_join2():
    struct = pd.DataFrame(
            {'struct_id': ['K0', 'K1', 'K2', 'K3', 'K4', 'K5'],
                'run_id':['0', '0', '1', '1', '2', '2'],
               'w_energy': ['A0', 'A1', 'A2', 'A3', 'A4', 'A5']}
            )

    sl = pd.DataFrame(
            {'struct_id': ['K1', 'K3', 'K1'],
               'loop_id': ['C1', 'C3', 'C4'],
                'run_id':['0', '1', '2'],
               'win_id': ['D1', 'D3', 'D1']}
            )
    result = pd.merge(struct, sl, how="inner", on=['struct_id', 'run_id']) 
    print "L"
    print struct
    print "R"
    print sl
    print "join"
    print result


def check_overlapping_sl():
    d = "/big_disk/capsid/packaging_signal/ca-ps_pipleine/pipeline_test/large/merge_orig_restart/run_again/"
    sl = pd.read_csv(d+"sl.tsv", sep="\t", header=0)
    seq = pd.read_csv(d+"seq.tsv", sep="\t", header=0)
    struct = pd.read_csv(d+"struct.tsv", sep="\t", header=0)
    win = pd.read_csv(d+"win.tsv", sep="\t", header=0)

    join_seq_sl = pd.merge(seq, sl, how="inner", on=["run_id", "seq_id"])
    join_seq_sl_struct = pd.merge(join_seq_sl, struct, how="inner", on=["run_id", "seq_id", "struct_id"])

    #print sl.shape, join_seq_sl.shape, join_seq_sl_struct.shape
    b = join_seq_sl_struct[ join_seq_sl_struct["group"] == "B"  ].sort_values(["start_a", "end_a", "loop_number"])
    
    print b.head()
    print b.tail()

    b.to_csv("b.tmp.tsv", sep="\t", index_label="index")
    print "done"

if __name__ == "__main__":
    pd.set_option('display.width', 140) # set width before introducing line breaks when printing df
    check_overlapping_sl()
