version 1.0

workflow DMRNAseq {
	input {
        Array[File] input_bams
        Array[String] sample_ids

        String project_id

        File annotation_gff
		File ref_fasta

        Boolean unaligned = true

        Int htseq_minMAPQ = 10
        String htseq_mode = "union"
        String htseq_nonunique = "none"
        String htseq_feature_type = "exon"
        String htseq_feature_ID = "gene_id"
        String htseq_additional_flags = ""
        
        Int memSizeGB = 512
        Int threadCount = 32
	}

    scatter (i in range(length(input_bams))) {
        if (unaligned) {
            call samtools_fastq {
                input:
                    ubam = input_bams[i],
                    sample_id = sample_ids[i],
                    memSizeGB = memSizeGB,
                    threadCount = threadCount
            }

            call minimap2_rna_alignment {
                input:
                    fastq = samtools_fastq.fastq,
                    ref_fasta = ref_fasta,
                    sample_id = sample_ids[i],
                    memSizeGB = memSizeGB,
                    threadCount = threadCount
            }
        } 

        call samtools_sort_index {
            input:
                sam = select_first([minimap2_rna_alignment.aligned_sam, input_bams[i]]),
                ref_fasta = ref_fasta,
                sample_id = sample_ids[i],
                memSizeGB = memSizeGB,
                threadCount = threadCount
        }

        call nanoplot_qc {
            input:
                bam = samtools_sort_index.sorted_bam, 
                sample_id = sample_ids[i],
                memSizeGB = memSizeGB,
                threadCount = threadCount
        }
    }

    call HTSeq_count {
        input:
            sorted_bams = samtools_sort_index.sorted_bam,
            annotation_gff = annotation_gff,
            project_id = project_id,
            htseq_minMAPQ = htseq_minMAPQ,
            htseq_mode = htseq_mode,
            htseq_nonunique = htseq_nonunique,
            htseq_feature_type = htseq_feature_type,
            htseq_feature_ID = htseq_feature_ID,
            htseq_additional_flags = htseq_additional_flags,
            memSizeGB = memSizeGB,
            threadCount = threadCount
    }

	output {
        File htseq_count_matrix = HTSeq_count.htseq_count_matrix
        Array[File] htseq_tagged_bams = HTSeq_count.htseq_tagged_bams

        Array[File] nanoplot_gzips = nanoplot_qc.nanoplot_gzip
	}

	meta {
		author: "Julian Menendez"
		email: "jmmenend@ucsc.edu"
		description: "Direct RNA differential expression & methylation pipeline."
	}
}

task samtools_fastq {
	input {
		File ubam
		String sample_id

        Int memSizeGB
		Int threadCount 
	}

	# Estimate disk size required
	Int input_bam_size = ceil(size(ubam, "GB"))       
	Int final_disk_dize = input_bam_size * 10

	# set outputs as wdl variables
	String fastq_output = "~{sample_id}.fastq"

	command <<<
		set -eux -o pipefail

		samtools fastq \
            -T * \
            -@ ~{threadCount} \
            ~{ubam} > ~{fastq_output}
	>>>

	output {
		File fastq = fastq_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/samtools:1.20" 
		preemptible: 1
	}
}

task minimap2_rna_alignment {
	input {
		File fastq
		File ref_fasta
		String sample_id

        Int memSizeGB
        Int threadCount
	}

    # Estimate disk size required
	Int input_fastq_size = ceil(size(fastq, "GB"))       
	Int final_disk_dize = input_fastq_size * 10

	# set outputs as wdl variables
	String aligned_sam_output = "~{sample_id}_minimap2.sam"

	command <<<
        set -eux -o pipefail

        minimap2 -ax splice \
            -u f \
            -k 13 \
            --eqx \
            -y \
            -Y \
            -t 32 \
            -N 20 \
            --junc-bed ${anno_bed} \
            -I 16G \
            ~{ref_fasta} ~{fastq} > ~{aligned_sam_output}
	>>>

	output {
		File aligned_sam = aligned_sam_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/minimap2:2.28"
		preemptible: 1
	}
}

task samtools_sort_index {
	input {
		File sam
        File ref_fasta
		String sample_id

        Int memSizeGB
        Int threadCount
	}

    # Estimate disk size required
	Int input_sam_size = ceil(size(sam, "GB"))       
	Int final_disk_dize = input_sam_size * 10

	# set outputs as wdl variables
	String sorted_bam_output = "~{sample_id}_minimap2.bam"
    String sorted_bai_output = "~{sample_id}_minimap2.bam.bai"

	command <<<
        set -eux -o pipefail

        samtools view -@ ~{threadCount} -bh -T ~{ref_fasta} ~{sam} | \
			samtools sort -@ ~{threadCount} - > ~{sorted_bam_output}

        samtools index -@ ~{threadCount} ~{sorted_bam_output}
	>>>

	output {
		File sorted_bam = sorted_bam_output
        File sorted_bai = sorted_bai_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/samtools:1.20"
		preemptible: 1
	}
}

task nanoplot_qc {
	input {
		File bam
		String sample_id

        Int memSizeGB
        Int threadCount
	}

    # Estimate disk size required
	Int input_bam_size = ceil(size(bam, "GB"))       
	Int final_disk_dize = input_bam_size * 10

	# set outputs as wdl variables
	String nanoplot_folder_output = "~{sample_id}_nanoplot"
    String nanoplot_gzip_output = "~{sample_id}_nanoplot.gz"

	command <<<
        set -eux -o pipefail

        NanoPlot \
            --threads ~{threadCount} \
            --outdir ~{nanoplot_folder_output} \
            --prefix ~{sample_id} \
            --raw \
            --loglength \
            --huge \
            --color "olivedrab" \
            --format png \
            --plots dot \
            --N50 \
            --dpi 1200 \
            --bam ~{bam}

        gzip -c ~{nanoplot_folder_output} > ~{nanoplot_gzip_output}
	>>>

	output {
		File nanoplot_gzip = nanoplot_gzip_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/nanoplot:1.42.0"
		preemptible: 1
	}
}

task HTSeq_count {
	input {
		Array[File] sorted_bams
        File annotation_gff
        String project_id

        Int htseq_minMAPQ
        String htseq_mode
        String htseq_nonunique 
        String htseq_feature_type
        String htseq_feature_ID
        String htseq_additional_flags

        Int memSizeGB
        Int threadCount
	}

    # Estimate disk size required
	Int input_bam_size = ceil(size(sorted_bams, "GB"))       
	Int final_disk_dize = input_bam_size * 10

    # set inputs as wdl variables
    String sorted_bams_string_input = sep(" ", sorted_bams)

	# set outputs as wdl variables
    String htseq_tagged_bams_output = sub(sorted_bams_string_input, "\\.bam$","_htseq.sam")
    String htseq_count_matrix_output = "~{project_id}_htseq_count_matrix.tsv"

	command <<<
        set -eux -o pipefail

        htseq-count \
            --format=bam \
            --order=pos \
            --nprocesses=~{threadCount} \
            --a=~{htseq_minMAPQ} \
            --mode=~{htseq_mode}
            --nonunique=~{htseq_nonunique} \
            --type=~{htseq_feature_type} \
            --idattr=~{htseq_feature_ID} \
            ~{htseq_additional_flags} \
            --samout=~{htseq_tagged_bams_output} \
            --counts_output=~{htseq_count_matrix_output} \
            ~{sorted_bams_string_input} ~{annotation_gff} # FIX SORTED BAMS
	>>>

	output {
		Array[File] htseq_tagged_bams = glob("*_htseq.sam")
        File htseq_count_matrix = htseq_count_matrix_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/htseq:2.0.5"
		preemptible: 1
	}
}