import argparse
import re
import pandas as pd

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Generate datasource for Tableau primer monitor dashboard",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-primer', type=str, required=True, help="input primer info")
    parser.add_argument('--input-nextclade', type=str, required=True, help="input Nextclade results")
    parser.add_argument('--input-metadata', type=str, required=True, help="input raw metadata")
    parser.add_argument('--param-clade', type=str, required=True, help="define clade build")
    parser.add_argument('--output-metadata', type=str, required=True, help="output metadata")
    args = parser.parse_args()

    print("Extracting primer regions:")
    pi = open(args.input_primer).readlines()
    primer = {}
    for i in range(0,len(pi)):
        line = pi[i].split("\t")
        item = line[3]
        primer[item] = list(range(int(line[1])+1,int(line[2])+1))

    print("Wrangling Nextclade results and metadata:")
    nextclade = pd.read_csv(args.input_nextclade,sep='\t',usecols=["seqName","clade","lineage","substitutions","deletions","insertions"])
    nextclade.rename(columns={"seqName":"accession"},inplace=True)
    nextclade.fillna({"substitutions":"0","deletions":"0","insertions":"0"},inplace=True)
    metadata = pd.read_csv(args.input_metadata,sep='\t',usecols=["accession","date","date_submitted","region","country","division"])
    metadata.rename(columns={"date_submitted":"submission"},inplace=True)
    metadata.fillna("Unknown",inplace=True)
    df = pd.merge(metadata,nextclade,on="accession")
    if args.param_clade == "Clade-I":
        df = df.loc[(df["clade"]).isin(["I","Ia","Ib"])].reset_index(drop=True)
    elif args.param_clade == "Clade-II":
        df = df.loc[(df["clade"]).isin(["II","IIa","IIb"])].reset_index(drop=True)

    print("Extracting variant positions:")
    var_pos = {}
    for i in range(0,len(df)):
        deletion = []
        del_range = df.loc[i,"deletions"].split(",")
        for j in range(0,len(del_range)):
            deletion = deletion + list(range(int(re.sub("-.*","",del_range[j])),int(re.sub(".*-","",del_range[j]))+1))
        snp = re.findall(r"\d+",df.loc[i,"substitutions"])
        insertion = re.findall(r"\d+",df.loc[i,"insertions"])
        mut = list(map(int,(snp+deletion+insertion)))
        mut = sorted(list(set(mut)))
        var_pos[df.loc[i,"accession"]] = mut

    print("Generating final variant summary file:")
    data = df[["accession","date","submission","region","country","division","clade","lineage"]].assign(set="",primer="",position="")
    for i in range(0,len(primer)):
        overlap = [list(filter(lambda x: x in list(primer.values())[i], sublist)) for sublist in list(var_pos.values())]
        for j in range(0,len(overlap)):
            for k in overlap[j]:
                data.loc[len(data.index)] = list(data.iloc[j,0:8])+[list(primer.keys())[i][:-2],list(primer.keys())[i][-1:],k]
    data.to_csv(args.output_metadata,sep='\t',index=False)
