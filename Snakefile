if not config:
    configfile: "config/config.yaml"

BUILD = config["tableau_clades"]
PRIMERS = config["igv_primers"]
MAIN = config["igv_main"]
SUB1 = config["igv_sub1"]
SUB2 = config["igv_sub2"]

rule all:
    input:
        primer=expand("tableau/{datasource}/monitor_primers.tsv", datasource=BUILD),
        tableau=expand("tableau/{datasource}/monitor_{datasource}.tsv", datasource=BUILD),
        igv=expand("igv/{clade}/{primer}.html", clade=MAIN, primer=PRIMERS),

rule import_data:
    "Import GenBank Mpox Nextclade metadata"
    input:
        primers=expand("primer/{datasource}/{primer}.bed", datasource=BUILD, primer=PRIMERS),
        nextclade=config["nextclade"],
        metadata=config["metadata"],
    output:
        primer=f"data/{BUILD}/primers.bed",
        nextclade=f"data/{BUILD}/nextclade.tsv",
        metadata=f"data/{BUILD}/metadata.tsv",
    shell:
        """
        cat {input.primers} > {output.primer}
        cp {input.nextclade} {output.nextclade}
        cp {input.metadata} {output.metadata}
        """

rule connect_tableau_primers:
    "Connect primer information to Tableau monitor dashboard"
    input:
        primer=rules.import_data.output.primer,
    output:
        primer=f"tableau/{BUILD}/monitor_primers.tsv",
    shell:
        """
        python3 script/tableau_primers.py \
            --input-primer {input.primer} \
            --output-primer {output.primer}
        """

rule generate_tableau_datasource:
    "Generate datasource for Tableau primer monitor dashboard"
    input:
        primer=rules.import_data.output.primer,
        nextclade=rules.import_data.output.nextclade,
        metadata=rules.import_data.output.metadata,
    params:
        clade=config["tableau_clades"],
    output:
        metadata=f"tableau/{BUILD}/monitor_{BUILD}.tsv",
    shell:
        """
        python3 script/tableau_datasource.py \
            --input-primer {input.primer} \
            --input-nextclade {input.nextclade} \
            --input-metadata {input.metadata} \
            --param-clade {params.clade} \
            --output-metadata {output.metadata}
        """

rule construct_cluster_variants:
    "Construct BED files of Mpox variants for clusters"
    input:
        nextclade=rules.import_data.output.nextclade,
    params:
        main=config["igv_main"],
        sub1=config["igv_sub1"],
        sub2=config["igv_sub2"],
        chrom=config["reference"],
    output:
        main=f"cluster/{MAIN}/{MAIN}.bed",
        sub1=f"cluster/{MAIN}/{SUB1}.bed",
        sub2=f"cluster/{MAIN}/{SUB2}.bed",
    shell:
        """
        python3 script/igv_variants.py \
            --input-nextclade {input.nextclade} \
            --param-main {params.main} \
            --param-sub1 {params.sub1} \
            --param-sub2 {params.sub2} \
            --param-chrom {params.chrom}
        """

rule search_variant_intersection:
    "Search genomic intersection between primer and cluster variants"
    input:
        primer=f"primer/{BUILD}/{{PRIMERS}}.bed",
        clade=rules.construct_cluster_variants.output.main,
    output:
        intersect=f"tmp/{MAIN}/{{PRIMERS}}_{MAIN}_intersects.bed",
        nonintersect=f"tmp/{MAIN}/{{PRIMERS}}_{MAIN}_nonintersect.bed",
    shell:
        """
        module load BEDTools
        bedtools intersect -wa -wb \
            -a {input.primer} -b {input.clade} \
            > {output.intersect}
        bedtools subtract -A \
            -a {input.primer} -b {input.clade} \
            > {output.nonintersect}
        """

rule assign_impact_score:
    "Assign impact score on primer regions"
    input:
        intersect=rules.search_variant_intersection.output.intersect,
        nonintersect=rules.search_variant_intersection.output.nonintersect,
    output:
        monitor_raw=f"monitor/{MAIN}/{{PRIMERS}}_raw.bed",
        monitor=f"monitor/{MAIN}/{{PRIMERS}}.bed",
        summary=f"monitor/{MAIN}/{{PRIMERS}}_overlap.bed",
    shell:
        """
        python3 script/igv_scores.py \
            --input-intersect {input.intersect} \
            --input-nonintersect {input.nonintersect} \
            --output-monitor {output.monitor_raw} \
            --output-summary {output.summary}
        sort -k4 {output.monitor_raw} -o {output.monitor}
        """

rule summarize_monitor_regions:
    "Summarize primer regions for IGV variant monitor"
    input:
        genome=config["genome"],
        monitor=rules.assign_impact_score.output.monitor,
    output:
        variant=f"variant/{MAIN}/{{PRIMERS}}.bed",
        site=f"variant/{MAIN}/{{PRIMERS}}.tsv",
    shell:
        """
        cat {input.genome} {input.monitor} > {output.variant}
        python3 script/igv_sites.py \
            --input-variant {output.variant} \
            --output-site {output.site}
        """

rule summarize_overlap_variants:
    "Summarize intersected variants as results summary"
    input:
        summary=rules.assign_impact_score.output.summary,
        sub1=rules.construct_cluster_variants.output.sub1,
        sub2=rules.construct_cluster_variants.output.sub2,
    params:
        primer=f"{{PRIMERS}}",
        reference=config["reference"],
        main=config["igv_main"],
        sub1=config["igv_sub1"],
        sub2=config["igv_sub2"],
    output:
        summary=f"summary/{MAIN}/{{PRIMERS}}.tsv",
    shell:
        """
        python3 script/igv_summary.py \
            --input-summary {input.summary} \
            --input-sub1 {input.sub1:q} \
            --input-sub2 {input.sub2:q} \
            --param-primer {params.primer} \
            --param-reference {params.reference} \
            --param-main {params.main} \
            --param-sub1 {params.sub1} \
            --param-sub2 {params.sub2} \
            --output-summary {output.summary}
        """

rule build_igv:
    "Build IGV primer monitor visualization browser"
    input:
        variant=rules.summarize_monitor_regions.output.site,
        summary=rules.summarize_overlap_variants.output.summary,
        reference=config["ref_fasta"],
        gene=config["ref_gff"],
        main=rules.construct_cluster_variants.output.main,
        sub1=rules.construct_cluster_variants.output.sub1,
        sub2=rules.construct_cluster_variants.output.sub2,
        monitor=rules.assign_impact_score.output.monitor,
    params:
        sequence=config["chrom"],
        begin=config["start"],
        end=config["end"],
        flank=config["flank"],
        template=config["template"],
    output:
        igv=f"igv/{MAIN}/{{PRIMERS}}.html",
    shell:
        """
        module load Python/3.8.3
        python script/igv_reports/report.py {input.variant} \
            --summary {input.summary} \
            --sequence {params.sequence} \
            --begin {params.begin} \
            --end {params.end} \
            --fasta {input.reference} \
            --tracks {input.gene} {input.main:q} {input.sub1:q} {input.sub2:q} {input.monitor} \
            --flanking {params.flank} \
            --template {params.template} \
            --output {output.igv}
        """

rule clean_build:
    "Clean up build directories"
    params:
        tableau=config["tableau_clades"],
        igv=config["igv_main"],
    shell:
        """
        rm -vrf "tmp/{params.igv}" "data/{params.tableau}" "cluster/{params.igv}" "monitor/{params.igv}" "variant/{params.igv}" "summary/{params.igv}" "igv/{params.igv}" "tableau/{params.tableau}"
        """

rule clean_all:
    "Clean up all directories"
    params:
        "tmp",
        "cluster",
        "monitor",
        "variant",
        "summary",
        "igv",
        "tableau",
    shell:
        """
        rm -vrf {params}
        """
