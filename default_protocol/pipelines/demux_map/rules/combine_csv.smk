rule combine_csv:
    """
    This rule takes the output CSV from map and variant calling and combines.
    """
    input:
        map_csv=config["output_path"] + "/temp/{filename_stem}_map.csv",
        vars_csv=config["output_path"] + "/temp/{filename_stem}_vars.csv"
    output:
        report = config["output_path"] + "/{filename_stem}.csv"
    shell:
        """
        fastafunk add_columns \
        --in-metadata {input.map_csv} \
        --in-data {input.vars_csv} \
        --index-column read_name \
        --join-on read_name \
        --new-columns variants \
        --out-metadata {output}
        """
#produces a csv report
