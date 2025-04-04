#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def parseTsvFileToMap(filePath) {
    def result = [:]
    def lines = new File(filePath).readLines()
    def header = lines[0].split("\t")

    lines.drop(1).eachWithIndex { line, index ->
        def values = line.split("\t")
        def lineMap = [:]
        header.eachWithIndex { colName, colIndex ->
            lineMap[colName] = values[colIndex]
        }

        result[lineMap['sample']] = lineMap
    }

    return result
}

library_info = parseTsvFileToMap(params.library_info)

process cellbender {

    cpus 1
    memory '64 GB'
    publishDir "${params.results}/cellbender_modified"
    container 'docker://porchard/cellbender:0.3.0'
    time '8h'
    errorStrategy 'ignore'

    input:
    tuple val(library), path(solo_out)

    output:
    path("${library}*")
    path("${library}*.h5"), emit: h5_files

    script:
    learning_rate = library_info[library].learning_rate
    expected_cells = library_info[library].expected_cells
    total_droplets_included = library_info[library].total_droplets_included

    """
    echo "library: ${library}"
    echo "solo_output Path: ${solo_out}"

    cp ${solo_out}/GeneFull_ExonOverIntron/raw/matrix.mtx matrix.mtx
    cp ${solo_out}/GeneFull_ExonOverIntron/raw/features.tsv genes.tsv
    cp ${solo_out}/GeneFull_ExonOverIntron/raw/barcodes.tsv barcodes.tsv

    cellbender remove-background --cuda --epochs 150 --fpr 0.05 0.1 0.15 --input . --output ./${library}.cellbender.h5 --learning-rate ${learning_rate} --expected-cells ${expected_cells} --total-droplets-included ${total_droplets_included}
    cp .command.log ${library}.log
    """

}


workflow {
    libraries = library_info.keySet()

    solo_out = Channel.from(libraries.collect({it -> [it, file(library_info[it].solo_out)]})) // library, solo_out

    cellbender(solo_out)

}
