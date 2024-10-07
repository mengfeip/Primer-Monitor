import argparse
import re
import pandas as pd

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Generate generic tab delimited results summary",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-summary', type=str, required=True, help="input raw summary")
    parser.add_argument('--input-sub1', type=str, required=True, help="sub-clade-1 variants profile")
    parser.add_argument('--input-sub2', type=str, required=True, help="sub-clade-2 variants profile")
    parser.add_argument('--param-primer', type=str, required=True, help="query primer set")
    parser.add_argument('--param-reference', type=str, required=True, help="define reference")
    parser.add_argument('--param-main', type=str, required=True, help="define main clade")
    parser.add_argument('--param-sub1', type=str, required=True, help="define sub-clade-1")
    parser.add_argument('--param-sub2', type=str, required=True, help="define sub-clade-2")
    parser.add_argument('--output-summary', type=str, required=True, help="output final summary")
    args = parser.parse_args()

    PRIMER = "Primer Set: " + args.param_primer
    Pv_MAIN = "Prevalence in " + args.param_main
    Pv_SUB1 = "Prevalence in " + args.param_sub1
    Pv_SUB2 = "Prevalence in " + args.param_sub2

    sub1 = pd.read_table(args.input_sub1,sep='\t',header=None,names=["Reference","Start","End","Site","Score","Strand","ThickStart","ThickEnd","RGB"])
    sub2 = pd.read_table(args.input_sub2,sep='\t',header=None,names=["Reference","Start","End","Site","Score","Strand","ThickStart","ThickEnd","RGB"])

    test = open(args.input_summary,'r').read().splitlines()
    if test == ["All Perfect Match"]:
        site = {PRIMER:['All Perfect Match'], 'Variant Description':['Not Found Any'], 'Impact Score':['0.000'], Pv_MAIN:['N/A'], Pv_SUB1:['N/A'], Pv_SUB2:['N/A'], 'Start':['N/A'], 'End':['N/A'], 'Reference':[args.param_reference]}
        summary = pd.DataFrame(site)
    else:
        site = pd.read_table(args.input_summary,sep='\t',header=None,names=["Reference","PrimerStart","PrimerEnd",PRIMER,"Raw Score","Strand","VarRef","Start","End","Mutation","Prevalence","VarStrand","ThickStart","ThickEnd","RGB","Blank"])
        site[PRIMER] = site[PRIMER].replace({".*-F$":" Forward",".*-P$":" Probe",".*-R$":" Reverse"},regex=True)
        for i in range(0,len(site)):
            site.loc[i,'Impact Score'] = f"{site.loc[i,'Raw Score']:.3f}"
            site.loc[i,Pv_MAIN] = f"{site.loc[i,'Prevalence']:.3f}%"
            if re.search(r"\-",site.loc[i,'Mutation']):
                site.loc[i,'Variant Description'] = "Deletion: " + re.sub("-","",site.loc[i,'Mutation']) + "bp"
            elif re.search(r"\+",site.loc[i,'Mutation']):
                site.loc[i,'Variant Description'] = "Insertion: " + re.sub("+","",site.loc[i,'Mutation']) + "bp"
            else:
                site.loc[i,'Variant Description'] = "SNP: " + str(site.loc[i,'End']) + site.loc[i,'Mutation']
            for j in range(0,len(sub1)):
                if ((site.loc[i,'End']==sub1.loc[j,'End']) and (site.loc[i,'Mutation']==sub1.loc[j,'Site'])):
                    site.loc[i,Pv_SUB1] = f"{sub1.loc[j,'Score']:.3f}%"
                    break
                else:
                    site.loc[i,Pv_SUB1] = "Not Present"
            for k in range(0,len(sub2)):
                if ((site.loc[i,'End']==sub2.loc[k,'End']) and (site.loc[i,'Mutation']==sub2.loc[k,'Site'])):
                    site.loc[i,Pv_SUB2] = f"{sub2.loc[k,'Score']:.3f}%"
                    break
                else:
                    site.loc[i,Pv_SUB2] = "Not Present"
        columns = [PRIMER,'Variant Description','Impact Score',Pv_MAIN,Pv_SUB1,Pv_SUB2,'Start','End','Reference']
        summary = site[columns]
        summary = summary.sort_values(by=['Variant Description',PRIMER])
    summary.to_csv(args.output_summary,sep='\t',index=False)
