#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.first_fastq = ''
params.second_fastq = ''
params.output_dir = ''
params.nanopore_fastq = ''
params.adapters = ''
params.kraken2_db = ''
params.gtdbtk_db = ''
params.utax_reference_db = ''
params.virsorter2_db = ''


process fastQC_hybrid {

    container 'staphb/fastqc:latest'
    publishDir "${params.output_dir}/fastqc_output", mode: 'copy'

    input:
    path first_fastq
    path second_fastq

    script:
    """
    mkdir -p ${params.output_dir}/fastqc_output
    fastqc -o ${params.output_dir}/fastqc_output ${first_fastq} ${second_fastq}
    """
}

process nanoplot_hybrid {

    container 'staphb/nanoplot:latest'
    publishDir "${params.output_dir}/nanoplot_output", mode: 'copy'

    input:
    path nanopore_fastq

    script:
    """
    mkdir -p ${params.output_dir}/nanoplot_output
    NanoPlot --fastq ${nanopore_fastq} -o ${params.output_dir}/nanoplot_output --plots hex dot
    """
}

process trimmomatic_hybrid {

    container 'staphb/trimmomatic:latest'

    input:
    path first_fastq
    path second_fastq
    path adapters

    output:
    path("trimmomatic_output/paired_R1.fastq.gz"), emit: forward_paired
    path("trimmomatic_output/paired_R2.fastq.gz"), emit: reverse_paired

    script:
    """
    mkdir -p trimmomatic_output
    trimmomatic PE ${first_fastq} ${second_fastq} \\
    trimmomatic_output/paired_R1.fastq.gz trimmomatic_output/forward_unpaired.fq.gz \\
    trimmomatic_output/paired_R2.fastq.gz trimmomatic_output/reverse_unpaired.fq.gz \\
    ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    cp -r trimmomatic_output ${params.output_dir}/
    """
}

process chopper_hybrid {

    container 'ufuomababatunde/chopper:v0.7.0'
    publishDir "${params.output_dir}/chopper_output", mode: 'copy'

    input:
    path nanopore_fastq
    
    output:
    path "chopper_output/nanopore_trimmed.fastq.gz", emit: trimmed_fastq_nanopore

    script:
    """
    mkdir -p chopper_output
    gunzip -c ${nanopore_fastq} | chopper -q 8 -l 500 | gzip > chopper_output/nanopore_trimmed.fastq.gz 
    """
}

process metaspades_hybrid {

    container 'staphb/spades:latest'
    publishDir "${params.output_dir}/metaspades_hybrid_output", mode: 'copy'

    input:
    path forward_paired
    path reverse_paired
    path trimmed_fastq_nanopore
    
    output:
    path "metaspades_hybrid_output/contigs.fasta", emit: metaspades_hybrid_contigs

    script:
    """
    mkdir -p metaspades_hybrid_output
    spades.py --meta -o metaspades_hybrid_output -1 ${forward_paired} -2 ${reverse_paired} --nanopore ${trimmed_fastq_nanopore} 
    """
}

process medaka_hybrid {

    container 'nanozoo/medaka:1.9.1--b3671e5'
    publishDir "${params.output_dir}/medaka_hybrid_output", mode: 'copy'

    input:
    path metaspades_hybrid_contigs
    path nanopore_fastq
    
    output:
    path "medaka_hybrid_output/consensus.fasta", emit: medaka_contigs

    script:
    """
    mkdir -p medaka_hybrid_output
    gunzip -c ${nanopore_fastq} > nanopore.fastq
    medaka_consensus -i nanopore.fastq -d ${metaspades_hybrid_contigs} -o medaka_hybrid_output -t 16 -m r941_min_high_g360 
    """
}

process polca_hybrid {

    container 'pnatarajan84/pypolca_revised:latest'
    publishDir "${params.output_dir}/polca_hybrid_output", mode: 'copy'
    
    input:
    path medaka_contigs
    path first_fastq
    path second_fastq
    
    output:
    path "polca_hybrid_output/polca_results/polca_corrected.fasta", emit: polca_contigs
    
    script: 
    """
    mkdir -p polca_hybrid_output
    pypolca run -a ${medaka_contigs} -1 ${first_fastq} -2 ${second_fastq} -t 4 -m 14G -o polca_hybrid_output/polca_results -p polca -f
    """
}

process metawrapBinning_hybrid {

    container 'pnatarajan84/metawrap_1.3.2_binning_fastq_gz_modified:latest'
    publishDir "${params.output_dir}/metawrap_output", mode: 'copy'

    input:
    path polca_contigs
    path forward_paired
    path reverse_paired

    output:
    path "metawrap_output/*", emit: metawrap_output_files
    path "metawrap_output/concoct_bins/*_bins_contigs.tsv", emit: concoct_bins_tsv
    path "metawrap_output/maxbin2_bins/*_bins_contigs.tsv", emit: maxbin2_bins_tsv
    path "metawrap_output/metabat2_bins/*_bins_contigs.tsv", emit: metabat2_bins_tsv

    script:
    """
    metawrap binning -o metawrap_output -t 40 -a ${polca_contigs} --metabat2 --maxbin2 --concoct ${forward_paired} ${reverse_paired}
    for dir in metawrap_output/concoct_bins metawrap_output/maxbin2_bins metawrap_output/metabat2_bins; do
        if [ -d "\$dir" ]; then
           bin_tool=\$(basename \$dir)
           tsv_file="\$dir/\${bin_tool}_bins_contigs.tsv"
           for bin_file in "\$dir"/*.fa; do
               bin_name=\$(basename \$bin_file .fa);
               grep "^>" \$bin_file | sed "s/>//" | awk -v bin="\$bin_name" '{print \$1 "\t" bin}' >> \$tsv_file;
           done;
       fi;
    done
    """
}	


process Semibin2_hybrid {
	container 'bladerunner2945/semibin2_minimap2_samtools:latest'

	input:
	path polca_contigs
	path nanopore_fastq

	output:
	path "semibin2_hybrid_output/semibin2_bins/output_bins", emit: semibin2_output_bins
	path "semibin2_hybrid_output/semibin2_bins/contig_bins.tsv", emit: semibin2_bins_tsv

	script:
	"""
	mkdir -p semibin2_hybrid_output
	minimap2 -ax map-ont ${polca_contigs} ${nanopore_fastq} > alignment_consensus.sam
	samtools view -bS alignment_consensus.sam | samtools sort -o alignment_consensus.sort.bam
	samtools index alignment_consensus.sort.bam
	SemiBin2 single_easy_bin -i ${polca_contigs} -b alignment_consensus.sort.bam -o semibin2_hybrid_output/semibin2_bins --environment wastewater

	for bin in semibin2_hybrid_output/semibin2_bins/output_bins/*.fa.gz; do
	gunzip \$bin
	done
	
	sed -i '1d' semibin2_hybrid_output/semibin2_bins/contig_bins.tsv
	chmod -R 777 semibin2_hybrid_output/semibin2_bins

	cp -r semibin2_hybrid_output/semibin2_bins ${params.output_dir}
	"""
	}


process DAStool_hybrid {
    
    container 'bladerunner2945/dastools1.1.6:latest'
    publishDir "${params.output_dir}/dastool_output", mode: 'copy'

    input:
    path concoct_bins_tsv
    path maxbin2_bins_tsv
    path semibin2_bins_tsv
    path metabat2_bins_tsv
    path polca_contigs

    output:
    path "dastool_output_DASTool_bins", emit: dastool_output_bins

    script:
    """
    # Ensure the output directory exists
    mkdir -p dastool_output

    # Run DAS_Tool
    DAS_Tool -i ${concoct_bins_tsv},${maxbin2_bins_tsv},${semibin2_bins_tsv},${metabat2_bins_tsv} -l concoct,maxbin2,semibin2,metabat2 -c ${polca_contigs} -o dastool_output --write_bin_evals --threads 20 --write_bins --write_unbinned

    """
}



process QUAST_hybrid {

    container 'staphb/quast:latest'
    publishDir "${params.output_dir}/quast_hybrid_output", mode: 'copy'
        
    input:
    path polca_contigs

    output:
    path "quast_hybrid_output", emit: quast_hybrid_report

    script:
    """
    mkdir -p quast_hybrid_output
    metaquast.py ${polca_contigs} -o quast_hybrid_output
    """
}	

process Kraken2_hybrid {

    container 'bladerunner2945/kraken2_krona:latest'
    publishDir "${params.output_dir}/kraken2_hybrid_output", mode: 'copy'

    input:
    path polca_contigs
    path kraken2_db
    
    output:
    path "kraken2_hybrid_output", emit: kraken2_hybrid_report

    script:
    """
    mkdir -p kraken2_hybrid_output
    kraken2 --db ${kraken2_db} --output kraken2_hybrid_output/kraken2_hybrid_output.txt --report kraken2_hybrid_output/kraken2_hybrid_report.txt ${polca_contigs}
    ktImportText -o kraken2_hybrid_output/kraken2_hybrid_krona_chart.html kraken2_hybrid_output/kraken2_hybrid_report.txt
    """
}

/*
process Metaphlan_hybrid {

	container 'stang/metaphlan4:v1'
	publishDir "${params.output_dir}/metaphlan_hybrid_output", mode: 'copy'

	input:
	path polca_contigs

	output:
	path "metaphlan_hybrid_output/metaphlan_hybrid_profile_contigs.txt", emit: metaphlan_profile_hybrid_contigs

	script:
	"""
	mkdir -p metaphlan_hybrid_output
	metaphlan ${polca_contigs} --input_type fasta > metaphlan_hybrid_output/metaphlan_hybrid_profile_contigs.txt
	"""
	}
*/



process CheckM_hybrid {

	container 'nanozoo/checkm:1.1.3--c79a047'
	publishDir "${params.output_dir}/checkm_hybrid_output", mode: 'copy'

	input:
	path dastool_output_bins

	output:
	path "checkm_output", emit: checkm_output_files

	script:
	"""
	mkdir -p checkm_hybrid_output
	checkm lineage_wf -x fa ${dastool_output_bins} checkm_output
	"""
	}

process GTDBTK_hybrid {

	container 'nanozoo/gtdbtk:2.3.2--641ec99'
	publishDir "${params.output_dir}/gtdbtk_hybrid_output", mode: 'copy'

	input:
	path dastool_output_bins
	path gtdbtk_db

	output:
	path "gtdbtk_hybrid_output", emit: gtdbtk_hybrid_output_files

	script:
	"""
	export GTDBTK_DATA_PATH=${gtdbtk_db}
	gtdbtk classify_wf --genome_dir ${dastool_output_bins} --out_dir gtdbtk_hybrid_output --skip_ani_screen --cpus 10 --extension fa
	"""
	}

process PrependBinNames_hybrid {
    
	    container 'bladerunner2945/prepend_bin_names:latest'
	    publishDir "${params.output_dir}/prepend_bin_names_output", mode: 'copy'

	    input: 
	    path dastool_output_bins

	    output:
	    path "merged_bins.fa", emit: merged_bins_fa

	    script:
	    """
	    prepend_bin_names.py ${dastool_output_bins} merged_bins.fa
	    """
	}

process Vsearch_hybrid  {

	container 'manutamminen/vsearch:v1.28'
        publishDir "${params.output_dir}/vsearch_fungi_output", mode: 'copy'
	
	input:
	path merged_bins_fa
	path utax_reference_db
	
	output:
	path "vsearch_fungi_output/output_hybrid_fungi_taxa_bins.txt", emit: vsearch_hybrid_fungi_output
	
	script:
	"""
	mkdir -p vsearch_fungi_output
	vsearch --usearch_global ${merged_bins_fa} --db ${utax_reference_db} --id 0.97 --blast6out vsearch_fungi_output/output_hybrid_fungi_taxa_bins.txt --top_hits_only
        """}
        
process AMRfinder_hybrid {

	container 'staphb/ncbi-amrfinderplus:latest'
	publishDir "${params.output_dir}/amrfinder_output", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "amrfinder_hybrid_output/amr_finder_hybrid_report.txt", emit: amr_finder_hybrid_report

	script:
	"""
	mkdir -p amrfinder_hybrid_output
	amrfinder -n ${merged_bins_fa} -o amrfinder_hybrid_output/amr_finder_hybrid_report.txt
	"""
	}  
	
	
	
process Deep_Arg_hybrid {

	container 'bladerunner2945/deep_arg:latest'
	publishDir "${params.output_dir}/deep_arg_output_hybrid", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "deep_arg_output_hybrid", emit: deep_arg_report_hybrid

	script:
	"""
	mkdir -p deep_arg_output_hybrid
	deeparg predict --model LS --input ${merged_bins_fa} --output deep_arg_output_hybrid/deep_arg_results_hybrid --d /data/deeparg_data
	"""
	}	      
	
	
process PlasmidFinder_hybrid {

	container 'staphb/plasmidfinder:latest'
	publishDir "${params.output_dir}/plasmidfinder_hybrid_output", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "plasmidfinder_hybrid_output/*", emit: plasmidfinder_hybrid_output_files

	script:
	"""
	mkdir -p plasmidfinder_hybrid_output
	plasmidfinder.py -i ${merged_bins_fa} -o plasmidfinder_hybrid_output -x -q
	"""
	}      
  
  
process VirSorter2_hybrid {

	container 'staphb/virsorter2:latest'
	publishDir "${params.output_dir}/virsorter2_hybrid_output", mode: 'copy'

	input:
	path merged_bins_fa
	path virsorter2_db

	output:
	path "virsorter2_hybrid_output", emit: virsorter2_hybrid_output_files


	script:
	"""
	mkdir -p virsorter2_hybrid_output
	virsorter run -i ${merged_bins_fa} -w virsorter2_hybrid_output --provirus-off --max-orf-per-seq 10 all -j 12 --db-dir ${virsorter2_db} --min-score 0.9
	"""
	}      
        

workflow  {

    fastqc_results = fastQC_hybrid(params.first_fastq, params.second_fastq)
    nanoplot_results =nanoplot_hybrid(params.nanopore_fastq)
    trimmomatic_results = trimmomatic_hybrid(params.first_fastq, params.second_fastq, params.adapters)
    chopper_results = chopper_hybrid(params.nanopore_fastq)
    metaspades_results = metaspades_hybrid(trimmomatic_results.forward_paired, trimmomatic_results.reverse_paired, chopper_results.trimmed_fastq_nanopore)
    medaka_results = medaka_hybrid(metaspades_results.metaspades_hybrid_contigs, params.nanopore_fastq)
    polca_results = polca_hybrid(medaka_results.medaka_contigs, params.first_fastq, params.second_fastq)
    quast_results = QUAST_hybrid(polca_results.polca_contigs)
    kraken2_results = Kraken2_hybrid(polca_results.polca_contigs, params.kraken2_db)
    //metaphlan_results = Metaphlan_hybrid(polca_results.polca_contigs)
    semibin2_results = Semibin2_hybrid(polca_results.polca_contigs, params.nanopore_fastq)
    binning_results = metawrapBinning_hybrid(polca_results.polca_contigs, trimmomatic_results.forward_paired, trimmomatic_results.reverse_paired)
    dastool_results = DAStool_hybrid(binning_results.concoct_bins_tsv, binning_results.maxbin2_bins_tsv, semibin2_results.semibin2_bins_tsv, binning_results.metabat2_bins_tsv, polca_results.polca_contigs)
    checkm_results = CheckM_hybrid(dastool_results.dastool_output_bins)
    gtdbtk_results = GTDBTK_hybrid(dastool_results.dastool_output_bins, params.gtdbtk_db)
    prepend_results = PrependBinNames_hybrid(dastool_results.dastool_output_bins)
    vsearch_results = Vsearch_hybrid(prepend_results.merged_bins_fa,params.utax_reference_db)
    amrfinder_results = AMRfinder_hybrid(prepend_results.merged_bins_fa)
    deeparg_results = Deep_Arg_hybrid(prepend_results.merged_bins_fa)
    plasmidFinder_results = PlasmidFinder_hybrid(prepend_results.merged_bins_fa)
    virsorter2_results = VirSorter2_hybrid(prepend_results.merged_bins_fa,params.virsorter2_db)
    
    }
