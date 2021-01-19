rule do_variant_calling:
    input:
        map=config["output_path"]+ "/temp/{filename_stem}_map.csv",
        fastq=get_unzipped_fastq
    params:
        path=workflow.current_basedir,
        config=workflow.current_basedir+"/../config.yaml",
        filename_stem=config["filename_stem"],
        output_path=config["output_path"],
        references_file=config["references_file"],
        variants_file=config["variants_file"],
        inbetween=config["output_path"] + "/{filename_stem}_all.csv"
    output:
        metadata=config["output_path"]+ "/{filename_stem}.csv"
    log:
        config["output_path"]+ "/temp/{filename_stem}.variant_call.log"
    run:
        if params.variants_file is not None and params.variants_file!= "":
            shell("snakemake --nolock --snakefile {params.path}/query_voc.smk "
                        "--configfile {params.config} "
                        "--config "
                        "filename_stem={params.filename_stem} "
                        "output_path={params.output_path} "
                        "references_file={params.references_file} "
                        "variants_file={params.variants_file} "
                        "reads_fastq={input.fastq} "
                        "mapped_csv={input.map} &> {log}")
            shell("mv {params.inbetween} {output.metadata}")
        else:
            with open(input.map,'r') as in_file, open(output.metadata,'w') as out_file:
                line=in_file.readline().strip()
                out_file.write("%s,variants\n" %line)
                line=in_file.readline().strip()
                while line:
                    out_file.write("%s,\n" %line)
                    line=in_file.readline().strip()


