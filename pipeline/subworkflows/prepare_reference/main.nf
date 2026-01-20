process BSBOLT_INDEX{

    tag "${reference_file}"

    input:
        path reference_file
        path panel_bed
    output:
        path "bisulfite_reference", type: "dir"
    script:
    """
        mkdir -p bisulfite_reference
        python -m bsbolt Index \
        -G ${reference_file} \
        -DB bisulfite_reference \
        -MR ${panel_bed}
    """
    stub:
    """
    mkdir -p bisulfite_reference
    """
}

process GET_CONTEXTS{
    tag "${reference_file}"

    input:
        path reference_file
        path panel_bed
    output:
        path "cpg_context.bed", emit: cpg
        path "conversion_context.bed", emit: conversion
    script:
    """
        python ${workflow.projectDir}/scripts/get_contexts.py \
        --input ${reference_file} \
        --panel ${panel_bed} \
        --out_cpg cpg_context.bed \
        --out_conversion conversion_context.bed
    """
    stub:
    """
    touch cpg_context.bed
    touch conversion_context.bed
    """
}

workflow PREPARE_REFERENCE{
    take:
        bisulfite_reference
    main:
        prepared_bisulfite_reference = bisulfite_reference.ifEmpty{
            BSBOLT_INDEX(params.usual_reference, params.panel_bed)
        }
        
        contexts = GET_CONTEXTS(params.usual_reference, params.panel_bed)
    emit:
        bisulfite_reference = prepared_bisulfite_reference
        cpg_context = contexts.cpg
        conversion_context = contexts.conversion
}