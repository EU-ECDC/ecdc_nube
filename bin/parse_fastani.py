#!/usr/bin/env python

import pandas as pd
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("-in", "--inp", type=str, 
                    help="fasta file with the sequences to add to the alignment", required=True)
parser.add_argument("-d", "--data_summary", type=str,
                    help="path to the curated data summary csv", required=True)
args = parser.parse_args()

# Load the results from fastANI
tab_r = pd.read_csv("fastani_out.txt", sep="\t", header = None)
tab_r.columns = ["in","ref","ani","orthologous_matches","total_fragments"]
isolate_key = tab_r["in"][0].split("/")[-1].replace(".fasta","")

# Reshape the ref entries
tab_r["ref"] = [i.split("/")[-1].replace(".fna","") for i in tab_r["ref"].tolist()]

# Load the data on the species
tab_s = pd.read_csv(args.data_summary, sep="\t")
d_sci_names = dict(zip(tab_s["Assembly Accession"].tolist(), tab_s["Organism Scientific Name"].tolist()))

# I define the column species in tab_r
tab_r["species"] = [d_sci_names[i] for i in tab_r["ref"].tolist()]
tab_r["query"] = isolate_key

# Number of matches above the threshold
#threshold = 95
#tab_matches = tab_r.query('ani > @threshold')

tab_r[["query","ref","species","ani"]].to_csv(f"{isolate_key}_species.tsv", sep ="\t", index = False)