import argparse
import pandas as pd

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Generate generic tab delimited sites for IGV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-variant', type=str, required=True, help="input variants info")
    parser.add_argument('--output-site', type=str, required=True, help="output site columns")
    args = parser.parse_args()

    site = pd.read_table(args.input_variant,sep='\t',header=None,names=["Reference","Start","End","Site","Raw Score","Strand","ThickStart","ThickEnd","RGB"])
    
    site['Tracked Site'] = site['Site'].replace({"-F$":" Forward","-P$":" Probe","-R$":" Reverse"},regex=True)
    site['Impact Rank'] = ["High" if x>=1 else "Medium" if x>0 else "Low" for x in site['Raw Score']]
    for i in range(0,len(site)):
        site.loc[i,'Total Score'] = f"{site.loc[i,'Raw Score']:.3f}"
    site.loc[site['Site']=="MPXV Genome",'Strand'] = '.'
    site.loc[site['Site']=="OPG057/F13L",'Strand'] = '.'
    site.loc[site['Site']=="OPG065/E3L",'Strand'] = '.'
    site.loc[site['Site']=="OPG071/E9L",'Strand'] = '.'
    site.loc[site['Site']=="tag048/F13L",'Strand'] = '.'
    site.loc[site['Site']=="tag055/E3L",'Strand'] = '.'
    site.loc[site['Site']=="tag062/E9L",'Strand'] = '.'
    site.loc[site['Site']=="MPXV Genome",['Total Score','Impact Rank']] = None
    site.loc[site['Site']=="OPG057/F13L",['Total Score','Impact Rank']] = None
    site.loc[site['Site']=="OPG065/E3L",['Total Score','Impact Rank']] = None
    site.loc[site['Site']=="OPG071/E9L",['Total Score','Impact Rank']] = None
    site.loc[site['Site']=="tag048/F13L",['Total Score','Impact Rank']] = None
    site.loc[site['Site']=="tag055/E3L",['Total Score','Impact Rank']] = None
    site.loc[site['Site']=="tag062/E9L",['Total Score','Impact Rank']] = None
    columns = ['Tracked Site','Reference','Start','End','Strand','Total Score','Impact Rank']
    site_display = site[columns]
    site_display.to_csv(args.output_site,sep='\t',index=False)
