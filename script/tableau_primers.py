import argparse
import pandas as pd

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Connect primer information to Tableau monitor dashboard",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-primer', type=str, required=True, help="input primer info")
    parser.add_argument('--output-primer', type=str, required=True, help="output primer connection")
    args = parser.parse_args()

    print("Generating primer connection list:")
    pi = pd.read_table(args.input_primer,sep='\t',header=None,names=["Reference","Start","End","Info","Score","Strand"])
    pi['Length'] = pi['End'] - pi['Start']
    pi['Set'] = pi['Info'].str[:-2]
    pi['Primer'] = pi['Info'].str[-1:]
    pi = pi[["Set","Primer","Start","End","Strand","Length","Info"]]
    pi.to_csv(args.output_primer,sep='\t',index=False)
