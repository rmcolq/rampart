rule first_ref:
    """
    Extract FIRST reference in the panel for coordinates.
    """
    input:
        ref= config["references_file"]
    output:
        temp(config["output_path"] + "/temp/first_ref.fasta")
    run:
        from Bio import SeqIO
        records = list(SeqIO.parse(input[0], "fasta"))
        reference_seq = records[0]
        with open(output[0], "w") as output_handle:
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
        fastq= get_unzipped_fastq,
        ref= config["output_path"] + "/temp/first_ref.fasta"
    output:
        sam= temp(config["output_path"] + "/temp/{filename_stem}.ref_aligned.sam"),
        fasta= temp(config["output_path"] + "/temp/{filename_stem}.ref_aligned.fasta")
    threads: config["threads"]
    shell:
        """
        minimap2 -a -x map-ont \
        --cs \
        --sam-hit-only \
        {input.ref:q} \
        {input.fastq:q} \
        > {output.sam:q}

        apps/gofasta/gofasta sam toMultiAlign \
        --reference {input.ref:q} \
        --samfile {output.sam:q} \
        > {output.fasta:q}
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
    threads: config["threads"]
    shell:
        """
        python3 apps/type_variants/type_variants.py \
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
        config["output_path"] + "/temp/{filename_stem}_vars.temp.csv"
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
#produces a csv report
