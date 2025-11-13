/*
 * Reusable module for concatenating files
 *
 * This module handles concatenation of large numbers of files (15k+)
 * with optional comment stripping and deduplication
 */

process CONCATENATE_FILES {
    tag "${tag}"
    cpus 1
    memory '4 GB'
    time '1 h'

    input:
    tuple path(input_files), val(pattern), val(output_name), val(strip_comments), val(sort_dedupe), val(tag)

    output:
    path "${output_name}"

    script:
    def strip_cmd = strip_comments ? "grep -h -v '^#'" : "cat"
    def sort_cmd = sort_dedupe ? "| sort -u" : ""

    """
    # Concatenate files matching pattern (exclude output file to avoid circular dependency)
    find . -name '${pattern}' -type l -not -name '${output_name}' -print0 | xargs -0 ${strip_cmd} ${sort_cmd} > ${output_name}

    # Ensure output file exists (even if empty)
    touch ${output_name}
    """
}
