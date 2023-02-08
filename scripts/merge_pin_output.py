
import pandas as pd
import argparse
from matplotlib import pyplot as plt

def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("-t", "--target",
						help="target search results file (.pin) format",
						type=str,required=True)
	parser.add_argument("-d", "--decoy",
						help="decoy search results file (.pin) format",
						type=str,required=True)
	parser.add_argument("-o", "--output",
						help="output file (.pin) format",
						type=str, required=True)
	parser.add_argument("--score",
						help="discriminant score for target decoy competition",
						type=str, default="sim2")
	parser.add_argument("--update_delta_scores",
						help="Update delta scores when target/decoy are 1st and 2nd ranked hit",
						action='store_true')
	parser.add_argument("--main_features",
						help="Track main features only. Otherwise all features will be tracked, which can reduce separation performance.",
						action='store_true')
 


	args = parser.parse_args()
	return args


def update_delta_scores(df1, idx1, df2, idx2):
    if df1.at[idx1, "delta_similarity"] > df1.at[idx1, "similarity"] - df2.at[idx2, "similarity"]:
        df1.at[idx1, "delta_similarity"] = df1.at[idx1, "similarity"] - df2.at[idx2, "similarity"]
    if df1.at[idx1, "delta_dot"] > df1.at[idx1, "dot_product"] - df2.at[idx2, "dot_product"]:
        df1.at[idx1, "delta_dot"] = df1.at[idx1, "dot_product"] - df2.at[idx2, "dot_product"]
    if df1.at[idx1, "delta_annotation_similarity"] > df1.at[idx1, "annotation_similarity"] - df2.at[idx2, "annotation_similarity"]:
        df1.at[idx1, "delta_annotation_similarity"] = df1.at[idx1, "annotation_similarity"] - df2.at[idx2, "annotation_similarity"]
    if df1.at[idx1, "delta_sim2"] > df1.at[idx1, "sim2"] - df2.at[idx2, "sim2"]:
        df1.at[idx1, "delta_sim2"] = df1.at[idx1, "sim2"] - df2.at[idx2, "sim2"]
    return

def merge_files(args):
    df = pd.read_csv(args.target, sep='\t', comment='#', low_memory=False)
    df_decoy = pd.read_csv(args.decoy, sep='\t', comment='#', low_memory=False)

    scans = df["ScanNr"].unique()
    
    for num in scans:
        decoy_match = df_decoy[df_decoy["ScanNr"] == num]
        if len(decoy_match) == 0:
            continue
        elif len(decoy_match) > 1:
            print("Error: multiple occurance of a ScanNr")
            exit(1)
        else:
            decoy_idx = decoy_match.index[0]
            decoy_match = decoy_match.iloc[0]
        
        target_match = df[df["ScanNr"] == num]
        target_idx = target_match.index[0]
        target_match = target_match.iloc[0]
        
        # Equal peptide -> Drop decoy
        if target_match["Peptide"].replace("L", "I") == decoy_match["Peptide"].replace("L", "I"):
            df_decoy.drop(decoy_idx, inplace=True)
            continue
        
        #print(target_idx, decoy_idx)
        # Compare score -> Keep higher scoring match
        if target_match[args.score] > decoy_match[args.score]:
            if args.update_delta_scores:
                update_delta_scores(df, target_idx, df_decoy, decoy_idx)
            df_decoy.drop(decoy_idx, inplace=True)
        else:
            if args.update_delta_scores:
                update_delta_scores(df_decoy, decoy_idx, df, target_idx)
            df.drop(target_idx, inplace=True)
        

    
    df = pd.concat([df, df_decoy], ignore_index=True)
    #df["sim2_half"] = df["sim2"] / 2.0
    #df["sim2_double"] = 2.0 * df["sim2"]
    
    #cols = ["PSMId", "Label", "ScanNr", "sim2", "sim2_half", "sim2_double", "Peptide", "Proteins"]
    #df = df[cols]
    
    if args.main_features:
        features = ["charge", "similarity", "bias", "delta_similarity", "annotation_similarity", "annotation_bias", "annotation_sim2", "delta_annotation_similarity", "peak_count_ref", "abs_mass_difference", "peptide_length"]
        col = ["PSMId", "Label", "ScanNr"] + features + ["Peptide", "Proteins"]
        df = df[col]
    df.to_csv(args.output, sep="\t", index=False)
	
    num_targets = sum(df["Label"] == 1)
    num_decoys = sum(df["Label"] == -1)
    
    
    print(f"Files merged successfully! {num_targets} targets and {num_decoys} decoys remaining after competition")



def main():
	args = parse_args()
	merge_files(args)
 
if __name__ == "__main__":
    main()


# bug fix: Remove -inf values
# sed -i 's/-inf/-9999/g' yeast_td.pin