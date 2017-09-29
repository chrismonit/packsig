""" Methods for manipulating tables with pandas, mostly one-off functions"""


import pandas as pd
import numpy as np

pd.set_option('display.width', 140) # set width before introducing line breaks when printing df

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
    
    

if __name__ == "__main__":
    get_names()
