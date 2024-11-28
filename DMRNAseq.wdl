version 1.0

workflow DMRNAseq {
    input {
        Array[File] input_bams

        Array[String] sample_ids
		Array[String] treatments

        String project_id

        File annotation_gff
        File ref_fasta

        Boolean unaligned = true

        # HTSEQ FLAGS
        Int htseq_minMAPQ = 10
        String htseq_mode = "union"
        String htseq_nonunique = "none"
        String htseq_feature_type = "exon"
        String htseq_feature_ID = "Parent"

        # MODKIT FLAGS
        Array[String] mod_codes = ["a", "m", "17802"]
        Array[String] ml_thresholds = ["205", "205", "205"]

		# EXTRA FLAGS OPTIONS
		String extra_samtools_flags = ""
        String extra_htseq_flags = ""		
		String extra_modkit_flags = ""

        Int memSizeGB = 256
        Int threadCount = 16
    }

    call gff_to_bed {
        input:
            annotation_gff = annotation_gff,
            memSizeGB = memSizeGB,
            threadCount = threadCount
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
                    anno_bed = gff_to_bed.anno_bed,
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

        call bedtools_coverage {
            input:
                bam = samtools_sort_index.sorted_bam, 
                sample_id = sample_ids[i],
                memSizeGB = memSizeGB,
                threadCount = threadCount
        }

        call HTSeq_count {
            input:
                sorted_bam = samtools_sort_index.sorted_bam,
                annotation_gff = annotation_gff,
                sample_id = sample_ids[i],
                htseq_minMAPQ = htseq_minMAPQ,
                htseq_mode = htseq_mode,
                htseq_nonunique = htseq_nonunique,
                htseq_feature_type = htseq_feature_type,
                htseq_feature_ID = htseq_feature_ID,
                htseq_additional_flags = htseq_additional_flags,
                memSizeGB = memSizeGB,
                threadCount = threadCount
        }
    }

    call DESeq2 {
        input:
            count_matrices = HTSeq_count.htseq_count_matrix,
            project_id = project_id,
			memSizeGB = memSizeGB,
			threadCount = threadCount
    }

    # Scatter over mod_codes and ml_thresholds
    scatter (j in range(length(mod_codes))) {
        scatter (i in range(length(input_bams))) {
			call ml_plot {
				input:
					bam = samtools_sort_index.sorted_bam[i], 
					mod_code = mod_codes[j],
					ml_threshold = ml_thresholds[j],
					sample_id = sample_ids[i],
					memSizeGB = memSizeGB,
					threadCount = threadCount
			}

            call modkit_pileup_by_htseq_tag {
                input:
                    bam = samtools_sort_index.sorted_bam[i],
                    bai = samtools_sort_index.sorted_bai[i],
                    ml_threshold = ml_thresholds[j],
                    mod_code = mod_codes[j],
                    extra_modkit_flags = extra_modkit_flags,
                    ref_fasta = ref_fasta,
                    sample_id = sample_ids[i],
                    memSizeGB = memSizeGB,
                    threadCount = threadCount
            }
        }

		Array[File] flattened_pileups = flatten(modkit_pileup_by_htseq_tag.modkit_pileups)

        call diff_methylation {
            input:
                pileups = flattened_pileups,
                mod_code = mod_codes[j],
                sample_ids = sample_ids,
                project_id = project_id,
				memSizeGB = memSizeGB,
				threadCount = threadCount
        }
    }

	Array[File] flattened_mlplots = flatten(flatten(ml_plot.ml_plots))

    output {
        Array[File] htseq_count_matrices = HTSeq_count.htseq_count_matrix
        Array[File] htseq_tagged_bams = HTSeq_count.htseq_tagged_bam

        Array[File] nanoplot_gzips = nanoplot_qc.nanoplot_gzip
		Array[File] ml_plots = flattened_mlplots
        Array[File] coverage_bedgraph = bedtools_coverage.coverage_bedgraph

        Array[File] diff_exp = DESeq2.differential_expression_matrices
		Array[File] diff_methyl_site = diff_methylation.diff_methyl_site
        Array[File] diff_methyl_region = diff_methylation.diff_methyl_region
    }

    meta {
        author: "Julian Menendez"
        email: "jmmenend@ucsc.edu"
        description: "Direct RNA differential expression & methylation pipeline."
    }
}

task gff_to_bed {
	input {
		File annotation_gff

        Int memSizeGB
		Int threadCount 
	}

	# Estimate disk size required
	Int input_gff_size = ceil(size(annotation_gff, "GB"))       
	Int final_disk_dize = input_gff_size * 6

	# set outputs as wdl variables
    String annotation_gff_basename = basename(annotation_gff)
	String anno_bed_output = "~{annotation_gff_basename}.bed"

	command <<<
		set -eux -o pipefail

		gff2bed < ~{annotation_gff} > ~{anno_bed_output}
	>>>

	output {
		File anno_bed = anno_bed_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "biocontainers/bedops:v2.4.35dfsg-1-deb_cv1" 
		preemptible: 1
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
            -T '*' \
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
		docker: "staphb/samtools:1.21" 
		preemptible: 1
	}
}

task minimap2_rna_alignment {
	input {
		File fastq
		File ref_fasta
        File anno_bed
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
            --junc-bed ~{anno_bed} \
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
		String extra_samtools_flags
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

        samtools view -@ ~{threadCount} ~{extra_samtools_flags} -bh -T ~{ref_fasta} ~{sam} | \
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
    String nanoplot_gzip_output = "~{sample_id}_nanoplot.tar.gz"

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

        tar -czf ~{nanoplot_gzip_output} ~{nanoplot_folder_output}
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

# make docker env for this
task ml_plot {
	input {
		File bam
		String mod_code
		Int ml_threshold
		String sample_id

        Int memSizeGB
        Int threadCount
	}

    # Estimate disk size required
	Int input_bam_size = ceil(size(bam, "GB"))       
	Int final_disk_dize = input_bam_size * 3

	command <<<
        set -eux -o pipefail

        # plot the global ML distribution of bam file
		MLplot.py ~{bam} ~{sample_id} ~{mod_code} 
	>>>

	output {
		Array[File] ml_plots = glob("*.png")
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/ML_plot: .... "
		preemptible: 1
	}
}

task bedtools_coverage {
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
	String coverage_bedgraph_output = "~{sample_id}_coverage.bedgraph"

	command <<<
        set -eux -o pipefail

        bedtools genomecov -bga -split -ibam ~{bam} > ~{coverage_bedgraph_output}
	>>>

	output {
		File coverage_bedgraph = coverage_bedgraph_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "staphb/bedtools:2.31.1"
		preemptible: 1
	}
}

task HTSeq_count {
	input {
		File sorted_bam
        File annotation_gff
        String sample_id

        Int htseq_minMAPQ
        String htseq_mode
        String htseq_nonunique 
        String htseq_feature_type
        String htseq_feature_ID
        String extra_htseq_flags

        Int memSizeGB
        Int threadCount
	}

    # Estimate disk size required
	Int input_bam_size = ceil(size(sorted_bam, "GB"))       
	Int final_disk_dize = input_bam_size * 10

	# set outputs as wdl variables
    String htseq_tagged_bam_output = "~{sample_id}_htseq.bam"
    String htseq_count_matrix_output = "~{sample_id}_htseq_matrix.tsv"

	command <<<
        set -eux -o pipefail

        htseq-count \
            --order pos \
            --stranded no \
            --minaqual ~{htseq_minMAPQ} \
            --type ~{htseq_feature_type} \
            --idattr ~{htseq_feature_ID} \
            --mode ~{htseq_mode} \
            --nonunique ~{htseq_nonunique} \
            --samout ~{htseq_tagged_bam_output} \
            --samout-format bam \
            --nprocesses ~{threadCount} \
            ~{extra_htseq_flags} \
            ~{sorted_bam} ~{annotation_gff} \
            --counts_output ~{htseq_count_matrix_output}
	>>>

	output {
		File htseq_tagged_bam = htseq_tagged_bam_output
        File htseq_count_matrix = htseq_count_matrix_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/htseq:2.0.9"
		preemptible: 1
	}
}


task DESeq2 {
	input {
		Array[File] count_matrices
		Array[String] sample_ids
		Array[String] treatments

		String project_id

		Int memSizeGB
        Int threadCount
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(count_matrices, "GB"))   
	Int final_disk_dize = input_file_size * 5

	String differential_expression_output_dir = "deseq2_results"
    String differential_expression_matrix = "~{project_id}.deseq2_results"

	command <<<
		set -eux -o pipefail

		# Create a coldata file with three columns: sample, condition, type
        printf "sample,condition,type\n" > coldata.csv
        
        # Iterate through sample_ids and treatments to populate coldata
        for i in $(seq 0 $((~{length(sample_ids)} - 1))); do
            sample_id=~{sample_ids[i]}
            treatment=~{treatments[i]}
            echo -e "${sample_id},${treatment},single-read" >> coldata.tsv
        done

		# Create output directory
        mkdir -p ~{differential_expression_output_dir}

		# create coldata from sample_ids, treatments, third column aka type is always single-read

        # Merge count matrices in bash
        # Assuming first column is gene names and subsequent columns are counts
        # Use awk to merge files based on the first column
        awk 'NR==FNR {a[$1]=$0; next} 
             {if(FNR==1) {h=$0} 
              else {if($1 in a) {print a[$1]"\t"$2} 
                    else {print $0"\t0"}}}' \
            ~{sep=' ' count_matrices} > merged_counts_matrix.tsv

        # Run R script for DESeq2 analysis
        Rscript /opt/run_DESeq2.R \
            merged_counts_matrix.tsv \
            ~{coldata_file} \
            ~{differential_expression_output_dir}/~{differential_expression_matrix}
	>>>

	output {
		Array[File] differential_expression_matrices = glob("~{differential_expression_output_dir}/*")
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/deseq2: ... "
		preemptible: 1
	}
}

task modkit_pileup_by_htseq_tag {
	input {
		File bam
		File bai
		File ref_fasta
		String sample_id

		String extra_modkit_flags
		String mod_code
		Int ml_threshold

		Int memSizeGB
        Int threadCount
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(bam, "GB"))   
	Int final_disk_dize = input_file_size * 10

	command <<<
		set -eux -o pipefail

		# turn ML score into 0-1 probability
		ml_prob=$(echo "scale=4; ~{ml_threshold}/256" | bc)

		# figure out how to fix motif
		modkit pileup \
			-t ~{threadCount} \
			--force-allow-implicit \
			--mod-thresholds a:${ml_prob} \
			--motif DRACH 2 \
			--ref ~{ref_fasta} \
			--prefix ~{sample_id} \
			--partition-tag XF \
			~{extra_modkit_flags} \
            ~{bam} ~{sample_id}_~{mod_code}_modkit_pileups

		for pileup in ~{sample_id}_~{mod_code}_modkit_pileups/*; do
			# Create bedgraph from filtered pileup files
			awk -F'\t' '{print $1, $2, $3, $11}' OFS='\t' ${pileup} > "${pileup}graph"
		done
	>>>

	output {
		Array[File] modkit_pileups = glob("~{sample_id}_~{mod_code}_modkit_pileups/*.bed")
		Array[File] modkit_pileup_bedgraphs = glob("~{sample_id}_~{mod_code}_modkit_pileups/*.bedgraph")
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "ontresearch/modkit:latest"
		preemptible: 1
	}
}

# write entire python script for this...
task diff_methylation {
	input {
		Array[File] pileups
		String mod_code
		Array[String] sample_ids
		String project_id

		Int memSizeGB
        Int threadCount
	}

	# Estimate disk size required
	Int input_file_size = ceil(size(pileups, "GB"))   
	Int final_disk_dize = input_file_size * 3

	String diff_methyl_site_output = "~{project_id}.dmr_site.tsv"
	String diff_methyl_region_output = "~{project_id}.dmr_regions.tsv"

	command <<<
		set -eux -o pipefail

	>>>

	output {
		File diff_methyl_site = diff_methyl_site_output
		File diff_methyl_region = diff_methyl_region_output
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "jmmenend/dmr: ...."
		preemptible: 1
	}
}