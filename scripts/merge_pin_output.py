
import pandas as pd
import argparse


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

	args = parser.parse_args()
	return args


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
            df_decoy.drop(decoy_idx, inplace=True)
        else:
            df.drop(target_idx, inplace=True)
        

    
    df = pd.concat([df, df_decoy], ignore_index=True)
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
# sed -i 's/-inf/-9999/g' yeast_td_test.pin