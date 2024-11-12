import argparse
import re
import math

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Calculate impact score of variants on primer regions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-intersect', type=str, required=True, help="input intersection between primer and variants")
    parser.add_argument('--input-nonintersect', type=str, required=True, help="input non-intersection between primer and variants")
    parser.add_argument('--output-monitor', type=str, required=True, help="output monitor .bed")
    parser.add_argument('--output-summary', type=str, required=True, help="output summary results")
    args = parser.parse_args()

    print("Assigning impact score to primer regions:")
    ist = open(args.input_intersect,'r').read().splitlines()
    nonist = open(args.input_nonintersect,'r').read().splitlines()
    intersect = {}
    summary = {}
    for i in ist:
        line = i.split("\t")
        primer_info = line[3]
        primer_start = int(line[1])
        primer_end = int(line[2])
        primer_strand = line[5]
        for j in range(1,len(line)):
            if line[j] == line[0]:
                break
        if re.search(r"\-",line[j+3]):
            variant_type = "del"
            variant_length = int(re.sub(r"\-","",line[j+3]))
        elif re.search(r"\+",line[j+3]):
            variant_type = "ins"
            variant_length = int(re.sub(r"\+","",line[j+3]))
        else:
            variant_type = "snp"
            variant_length = int(len(line[j+3]))
        variant_start = int(line[j+1])
        variant_end = int(line[j+2])
        overlap_start = variant_start if variant_start>=primer_start else primer_start
        overlap_end = variant_end if variant_end<=primer_end else primer_end
        overlap_length = overlap_end - overlap_start
        length_5p = primer_end-overlap_end if primer_strand=="-" else overlap_start-primer_start
        length_3p = overlap_end-primer_start-1 if primer_strand=="-" else primer_end-overlap_start-1
        total_score = 0
        for k in range(0,overlap_length):
            score = 0
            if length_5p <= 3:
                score += 1
            elif length_3p <= 2:
                score += 2
            elif 2<length_3p<=5:
                score += 4
            else:
                score += 3
            if variant_type == "ins":
                score *= 10 * variant_length
            elif variant_type == "del":
                score *= 10
            else:
                score *= 1
            total_score += score
            length_5p += 1
            length_3p -= 1
        summary[i] = line
        summary[i][4] = math.log(2 ** total_score)
        if primer_info not in intersect.keys():
            intersect[primer_info] = line[0:6]
            intersect[primer_info][4] = 2 ** total_score
        else:
            intersect[primer_info][4] += 2 ** total_score
    for i in list(intersect.keys()):
        intersect[i][4] = math.log(intersect[i][4])
    for i in nonist:
        line = i.split("\t")
        intersect[line[3]] = line[0:6]
        intersect[line[3]][4] = 0
    for i in list(intersect.keys()):
        if float(intersect[i][4]) > 1:
            line = "\t".join(list(map(str,intersect[i]))+[intersect[i][1],intersect[i][1],"255,0,0\n"])
            open(args.output_monitor,'a').writelines(line)
        elif 0<float(intersect[i][4])<=1:
            line = "\t".join(list(map(str,intersect[i]))+[intersect[i][1],intersect[i][1],"255,128,0\n"])
            open(args.output_monitor,'a').writelines(line)
        else:
            line = "\t".join(list(map(str,intersect[i]))+[intersect[i][1],intersect[i][1],"0,128,255\n"])
            open(args.output_monitor,'a').writelines(line)
    if len(summary)==0:
        open(args.output_summary,'a').writelines("All Perfect Match")
    else:
        for i in list(summary.keys()):
            summary_line = "\t".join(list(map(str,summary[i]))+["\n"])
            open(args.output_summary,'a').writelines(summary_line)
