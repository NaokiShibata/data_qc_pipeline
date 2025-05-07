nextflow.enable.dsl=2

params.input_csv = null

process FASTP {
    tag "$sample_id"
    publishDir "Results/fastp/fastq", pattern: "*_R{1,2}.fastq.gz", mode: "copy" 
    publishDir "Results/fastp/summary", pattern: "*.{html,json}", mode: "copy"

    input:
        tuple val(sample_id), path(read1), path(read2)
    
    output:
        tuple(path("${sample_id}_R{1,2}.fastq.gz")), emit: fastp_fastq
        path("${sample_id}_report.json"), emit: fastp_json
        path("${sample_id}_report.html"), emit: fastp_html

    script:
        """
        fastp \
            --in1 ${read1} \
            --in2 ${read2} \
            --out1 ${sample_id}_R1.fastq.gz \
            --out2 ${sample_id}_R2.fastq.gz \
            --html ${sample_id}_report.html \
            --json ${sample_id}_report.json \
            --qualified_quality_phred 30 \
            --length_required 50 \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --cut_front \
            --thread ${task.cpus}
        """
}

process MULTIQC {
    publishDir "Results/multiqc", mode: "copy"

    input:
        path fastp_json

    output:
        path "*"

    script:
        """
        multiqc .
        """
}

workflow {
    if (! params.input_csv) {
        error "Error: Please provide a CSV file with --input_csv"
    }

    reads = Channel.fromPath( params.input_csv )
        .splitCsv(header: true)
        .map { row -> 
            tuple(row.sample_id, file(row.read1), file(row.read2))
        }
    
    FASTP(reads)
    MULTIQC(FASTP.out.fastp_json.collect())
}