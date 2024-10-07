import argparse
import re
import pandas as pd

def variants(cluster):
    if re.fullmatch("all",cluster,re.IGNORECASE):
        file_name = "ALL"
        nc_sub = nc
    elif re.search("I",cluster):
        clade = re.search("I.*",cluster).group(0)
        file_name = cluster
        if clade == "I":
            nc_sub = nc[nc["clade"].isin(["I","Ia","Ib"])]
        elif clade == "II":
            nc_sub = nc[nc["clade"].isin(["II","IIa","IIb"])]
        else:
            nc_sub = nc[nc["clade"]==clade]
    elif re.search(r"\*",cluster):
        lineage = re.sub(r"\*","",cluster)
        file_name = cluster
        if lineage == "A":
            nc_sub = nc[nc["lineage"].str.contains("A",na=False)]
        elif lineage == "B":
            nc_sub = nc[nc["lineage"].str.contains("B|C",na=False)]
    else:
        file_name = cluster
        nc_sub = nc[nc["lineage"]==cluster]
    variant = ""
    snp = nc_sub["substitutions"].dropna().reset_index(drop=True)
    pool = []
    for i in range(0,len(snp.index)):
        pool = pool + str(snp[i]).split(",")
    freq = pd.Series(pool).value_counts()
    freq_show = freq[freq>=len(nc_sub.index)*0.05]
    if len(freq_show) > 0:
        for i in range(0,len(freq_show)):
            ref = "REFERENCE"
            end = int(re.search(r"\d+",freq_show.index[i]).group(0))
            start = end - 1
            mut = re.sub(r".*\d+","",freq_show.index[i])
            pct = freq_show.iloc[i] / len(nc_sub.index) * 100
            strand = "."
            color = "0,0,0\n"
            variant = variant + "\t".join(list(map(str,[ref,start,end,mut,pct,strand,start,start,color])))
    deletion = nc_sub["deletions"].dropna().reset_index(drop=True)
    pool = []
    for i in range(0,len(deletion.index)):
        pool = pool + str(deletion[i]).split(",")
    freq = pd.Series(pool).value_counts()
    freq_show = freq[freq>=2]
    if len(freq_show) > 0:
        for i in range(0,len(freq_show)):
            ref = "REFERENCE"
            start = int(re.sub("-.*","",freq_show.index[i])) - 1
            end = int(re.sub(".*-","",freq_show.index[i]))
            mut = "".join([str(end-start),"-"])
            pct = freq_show.iloc[i] / len(nc_sub.index) * 100
            strand = "."
            color = "0,204,0\n"
            variant = variant + "\t".join(list(map(str,[ref,start,end,mut,pct,strand,start,start,color])))
    insertion = nc_sub["insertions"].dropna().reset_index(drop=True)
    pool = []
    for i in range(0,len(insertion.index)):
        pool = pool + str(insertion[i]).split(",")
    freq = pd.Series(pool).value_counts()
    freq_show = freq[freq>=2]
    if len(freq_show) > 0:
        for i in range(0,len(freq_show)):
            ref = "REFERENCE"
            start = int(re.sub(":.*","",freq_show.index[i]))
            end = start + 1
            mut = "".join([str(len(re.sub(".*:","",freq_show.index[i]))),"+"])
            pct = freq_show.iloc[i] / len(nc_sub.index) * 100
            strand = "."
            color = "127,0,255\n"
            variant = variant + "\t".join(list(map(str,[ref,start,end,mut,pct,strand,start,start,color])))
    return file_name, variant

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Construct BED files of variants for Mpox clusters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-nextclade', type=str, required=True, help="input Nextclade results")
    parser.add_argument('--param-main', type=str, required=True, help="define main clade")
    parser.add_argument('--param-sub1', type=str, required=True, help="define sub-clade-1")
    parser.add_argument('--param-sub2', type=str, required=True, help="define sub-clade-2")
    parser.add_argument('--param-chrom', type=str, required=True, help="define variant reference")
    args = parser.parse_args()

    print("Generating variants .bed files:")
    nc = pd.read_csv(args.input_nextclade,sep='\t',usecols=["seqName","clade","lineage","substitutions","deletions","insertions"])
    clusters = [args.param_main,args.param_sub1,args.param_sub2]
    for i in clusters:
        print(i + " variants")
        file_name,variant = variants(cluster=i)
        variant = variant.replace("REFERENCE",args.param_chrom)
        open("".join(["cluster/",args.param_main,"/",file_name,".bed"]),'w').writelines(variant)
    if "Clade-IIb" in clusters:
        for i in nc["lineage"].value_counts().index:
            if i in ["Ia","Ib"]:
                continue
            print(i + " variants")
            file_name,variant = variants(cluster=i)
            variant = variant.replace("REFERENCE",args.param_chrom)
            open("".join(["cluster/Clade-IIb/",file_name,".bed"]),'w').writelines(variant)
