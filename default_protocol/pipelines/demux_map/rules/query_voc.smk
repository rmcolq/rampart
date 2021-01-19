from Bio import SeqIO

rule all:
    input:
        expand(config["output_path"]+ "/{filename_stem}_all.csv", filename_stem=config["filename_stem"])

rule first_ref:
    """
    Extract FIRST reference in the panel for coordinates.
    """
    input:
        ref=config["references_file"]
    output:
        first_ref=temp(config["output_path"] + "/temp/first_ref.fasta")
    run:
        records = list(SeqIO.parse(input.ref, "fasta"))
        reference_seq = records[0]
        with open(output.first_ref, "w") as output_handle:
            SeqIO.write(reference_seq, output_handle, "fasta")

rule read_to_ref_coordinates:
    """
    This rule takes the FASTQ and maps it to the FIRST reference in the panel.
    It uses the FASTQ independent of demuxing -- i.e. it doesn't matter whether the file has barcode
    information in the header.
    It translates the resulting match/mismatch cigar to a string in the SAM file to something with a 1-1
    relationship with the reference (ignores insertions).
    """
    input:
        ref= config["output_path"] + "/temp/first_ref.fasta",
        fastq= config["reads_fastq"]
    output:
        sam= temp(config["output_path"] + "/temp/{filename_stem}.ref_aligned.sam"),
        fasta= temp(config["output_path"] + "/temp/{filename_stem}.ref_aligned.fasta")
    shell:
        """
        minimap2 -a -x map-ont \
            --cs \
            --sam-hit-only \
            {input.ref:q} \
            {input.fastq:q} \
            > {output.sam:q}

        datafunk sam_2_fasta \
            -s {output.sam:q} \
            -r {input.ref:q} \
            -o {output.fasta:q} \
            --pad
        """

rule query_voc:
    """
    This rule queries a list of mutations to see if there is evidence in the reads.
    """
    input:
        fasta= config["output_path"] + "/temp/{filename_stem}.ref_aligned.fasta",
        ref= config["output_path"] + "/temp/first_ref.fasta",
        variants= config["variants_file"]
    output:
        temp(config["output_path"] + "/temp/{filename_stem}_vars.temp.csv")
    shell:
        """
        python3 type_variants.py \
            --fasta-in {input.fasta:q} \
            --variants-config {input.variants:q} \
            --reference {input.ref:q} \
            --variants-out {output:q} \
            --append-genotypes
        """

rule make_csv:
    """
    This makes a variants csv in correct format.
    """
    input:
        config["output_path"] + "/temp/{filename_stem}_vars.temp.csv",
    params:
        path_to_script = workflow.current_basedir
    output:
        temp(config["output_path"] + "/temp/{filename_stem}_vars.csv")
    shell:
        """
        python {params.path_to_script}/parse_voc.py \
            --in-file {input:q} \
            --out-file {output:q} \
            --random
        """

rule combine_csv:
    """
    This rule takes the output CSV from map and variant calling and combines.
    """
    input:
        vars_csv=config["output_path"] + "/temp/{filename_stem}_vars.csv",
        map_csv=config["mapped_csv"]
    output:
        config["output_path"] + "/{filename_stem}_all.csv"
    shell:
        """
        fastafunk add_columns \
            --in-metadata {input.map_csv} \
            --in-data {input.vars_csv:q} \
            --index-column read_name \
            --join-on read_name \
            --new-columns variants \
            --out-metadata {output:q}
        """