#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.nanopore_fastq = ''
params.output_dir = ''
params.kraken2_db = ''
params.gtdbtk_db = ''
params.utax_reference_db = ''
params.virsorter2_db = ''

process nanoplot {

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

process chopper {

    container 'ufuomababatunde/chopper:v0.7.0'
    publishDir "${params.output_dir}/chopper_output", mode: 'copy'

    input:
    path nanopore_fastq
    
    output:
    path "chopper_output/nanopore_trimmed.fastq.gz", emit: trimmed_fastq

    script:
    """
    mkdir -p chopper_output
    gunzip -c ${nanopore_fastq} | chopper -q 8 -l 500 | gzip > chopper_output/nanopore_trimmed.fastq.gz 
    """
}

process flye {

    container 'staphb/flye:latest'
    publishDir "${params.output_dir}/flye_output", mode: 'copy'

    input:
    path nanopore_fastq
    
    output:
    path "flye_output/assembly.fasta", emit: flye_contigs

    script:
    """
    mkdir -p flye_output
    flye --nano-raw ${nanopore_fastq} --meta --out-dir flye_output --threads 10 
    """
}

process medaka {

    container 'nanozoo/medaka:1.9.1--b3671e5'
    publishDir "${params.output_dir}/medaka_output", mode: 'copy'

    input:
    path flye_contigs
    path nanopore_fastq
    
    output:
    path "medaka_output/consensus.fasta", emit: medaka_contigs

    script:
    """
    mkdir -p medaka_output
    gunzip -c ${nanopore_fastq} > nanopore.fastq
    medaka_consensus -i nanopore.fastq -d ${flye_contigs} -o medaka_output -t 16 -m r941_min_high_g360 
    """
}




process QUAST {

    container 'staphb/quast:latest'
    publishDir "${params.output_dir}/quast_output", mode: 'copy'
        
    input:
    path medaka_contigs

    output:
    path "quast_output", emit: quast_report

    script:
    """
    mkdir -p quast_output
    metaquast.py ${medaka_contigs} -o quast_output
    """
}	


process Kraken2 {

    container 'bladerunner2945/kraken2_krona:latest'
    publishDir "${params.output_dir}/kraken2_output", mode: 'copy'

    input:
    path medaka_contigs
    path kraken2_db
    
    output:
    path "kraken2_output", emit: kraken2_report

    script:
    """
    mkdir -p kraken2_output
    kraken2 --db ${kraken2_db} --output kraken2_output/kraken2_output.txt --report kraken2_output/kraken2_report.txt ${medaka_contigs}
    ktImportText -o kraken2_output/kraken2_nanopore_krona_chart.html kraken2_output/kraken2_report.txt
    """
}

process Metaphlan_nanopore {


	container 'stang/metaphlan4:v1'
	publishDir "${params.output_dir}/metaphlan_output", mode: 'copy'

	input:
	path medaka_contigs

	output:
	path "metaphlan_output/metaphlan_profile_nanopore_contigs.txt", emit: metaphlan_profile_nanopore_contigs

	script:
	"""
	mkdir -p metaphlan_output
	metaphlan ${medaka_contigs} --input_type fasta > metaphlan_output/metaphlan_profile_nanopore_contigs.txt
	"""
	}
	
process Semibin2_nanopore {
	container 'bladerunner2945/semibin2_minimap2_samtools:latest'

	input:
	path medaka_contigs
	path nanopore_fastq

	output:
	path "semibin2_output/semibin2_bins/output_bins", emit: semibin2_output_bins
	path "semibin2_output/semibin2_bins/contig_bins.tsv", emit: semibin2_bins_tsv

	script:
	"""
	mkdir -p semibin2_output
	minimap2 -ax map-ont ${medaka_contigs} ${nanopore_fastq} > alignment_consensus.sam
	samtools view -bS alignment_consensus.sam | samtools sort -o alignment_consensus.sort.bam
	samtools index alignment_consensus.sort.bam
	SemiBin2 single_easy_bin -i ${medaka_contigs} -b alignment_consensus.sort.bam -o semibin2_output/semibin2_bins --environment wastewater
	
	sed -i '1d' semibin2_output/semibin2_bins/contig_bins.tsv
	chmod -R 777 semibin2_output/semibin2_bins
	
	cp -r semibin2_output/semibin2_bins ${params.output_dir}
	"""
	}
	
process metawrapBinning_nanopore {

    container 'pnatarajan84/metawrap_1.3.2_binning_fastq_gz_modified:latest'
    publishDir "${params.output_dir}/metawrap_output", mode: 'copy'

    input:
    path medaka_contigs
    path nanopore_fastq

    output:
    path "metawrap_output/*", emit: metawrap_output_files
    path "metawrap_output/concoct_bins/*_bins_contigs.tsv", emit: concoct_bins_tsv
    path "metawrap_output/maxbin2_bins/*_bins_contigs.tsv", emit: maxbin2_bins_tsv
    path "metawrap_output/metabat2_bins/*_bins_contigs.tsv", emit: metabat2_bins_tsv

    script:
    """
    metawrap binning -o metawrap_output -t 40 -a ${medaka_contigs} --metabat2 --maxbin2 --concoct --single-end ${nanopore_fastq}
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
	
process DAStool_nanopore {
    
    container 'bladerunner2945/dastools1.1.6:latest'
    publishDir "${params.output_dir}/dastool_output", mode: 'copy'

    input:
    path concoct_bins_tsv
    path maxbin2_bins_tsv
    path semibin2_bins_tsv
    path medaka_contigs

    output:
    path "dastool_output_DASTool_bins", emit: dastool_output_bins

    script:
    """
    # Ensure the output directory exists
    mkdir -p dastool_output

    # Run DAS_Tool
    DAS_Tool -i ${concoct_bins_tsv},${maxbin2_bins_tsv},${semibin2_bins_tsv} -l concoct,maxbin2,semibin2 -c ${medaka_contigs} -o dastool_output --write_bin_evals --threads 20 --write_bins --write_unbinned

    """
}	



process CheckM_nanopore {

	container 'nanozoo/checkm:1.1.3--c79a047'
	publishDir "${params.output_dir}/checkm_output", mode: 'copy'

	input:
	path dastool_output_bins

	output:
	path "checkm_output", emit: checkm_output_files

	script:
	"""
	mkdir -p checkm_output
	checkm lineage_wf -x fa ${dastool_output_bins} checkm_output
	"""
	}

process GTDBTK_nanopore {

	container 'nanozoo/gtdbtk:2.3.2--641ec99'
	publishDir "${params.output_dir}/gtdbtk_output", mode: 'copy'

	input:
	path dastool_output_bins
	path gtdbtk_db

	output:
	path "gtdbtk_output/*", emit: gtdbtk_output_files

	script:
	"""
	export GTDBTK_DATA_PATH=${gtdbtk_db}
	gtdbtk classify_wf --genome_dir ${dastool_output_bins} --out_dir gtdbtk_output --skip_ani_screen --cpus 20 --extension fa
	"""
	}
	
process PrependBinNames_nanopore {
    
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


process Vsearch_nanopore  {

	container 'manutamminen/vsearch:v1.28'
        publishDir "${params.output_dir}/vsearch_fungi_output", mode: 'copy'
	
	input:
	path merged_bins_fa
	path utax_reference_db
	
	output:
	path "vsearch_fungi_output/output_fungi_taxa_bins.txt", emit: vsearch_fungi_output
	
	script:
	"""
	mkdir -p vsearch_fungi_output
	vsearch --usearch_global ${merged_bins_fa} --db ${utax_reference_db} --id 0.97 --blast6out vsearch_fungi_output/output_fungi_taxa_bins.txt --top_hits_only
        """}
        
        
        
process AMRfinder_nanopore {

	container 'staphb/ncbi-amrfinderplus:latest'
	publishDir "${params.output_dir}/amrfinder_output", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "amrfinder_output/amr_finder_report.txt", emit: amr_finder_report

	script:
	"""
	mkdir -p amrfinder_output
	amrfinder -n ${merged_bins_fa} -o amrfinder_output/amr_finder_report.txt
	"""
	}
	
	
process Deep_Arg_nanopore {

	container 'bladerunner2945/deep_arg:latest'
	publishDir "${params.output_dir}/deep_arg_output_nanopore", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "deep_arg_output_nanopore", emit: deep_arg_report_nanopore

	script:
	"""
	mkdir -p deep_arg_output_nanopore
	deeparg predict --model LS --input ${merged_bins_fa} --output deep_arg_output_nanopore/deep_arg_results_nanopore --d /data/deeparg_data
	"""
	}	
	
	
process PlasmidFinder_nanopore {



	container 'staphb/plasmidfinder:latest'
	publishDir "${params.output_dir}/plasmidfinder_output", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "plasmidfinder_output/*", emit: plasmidfinder_output_files

	script:
	"""
	mkdir -p plasmidfinder_output
	plasmidfinder.py -i ${merged_bins_fa} -o plasmidfinder_output -x -q
	"""
	}
	
process VirSorter2_nanopore {

	container 'staphb/virsorter2:latest'
	publishDir "${params.output_dir}/virsorter2_output", mode: 'copy'

	input:
	path merged_bins_fa
	path virsorter2_db

	output:
	path "virsorter2_output", emit: virsorter2_output_files


	script:
	"""
	mkdir -p virsorter2_output
	virsorter run -i ${merged_bins_fa} -w virsorter2_output --provirus-off --max-orf-per-seq 10 all -j 12 --db-dir ${virsorter2_db} --min-score 0.9
	"""
	}
	
	

workflow  {

    nanoplot_results =nanoplot(params.nanopore_fastq)
    chopper_results = chopper(params.nanopore_fastq)
    flye_results = flye(params.nanopore_fastq)
    medaka_results = medaka(flye_results.flye_contigs, params.nanopore_fastq)
    quast_results = QUAST(medaka_results.medaka_contigs)
    kraken2_results = Kraken2(medaka_results.medaka_contigs, params.kraken2_db)
    //metaphlan4_results = Metaphlan_nanopore(medaka_results.medaka_contigs) 
    semibin2_results = Semibin2_nanopore(medaka_results.medaka_contigs, params.nanopore_fastq)
    binning_results = metawrapBinning_nanopore(medaka_results.medaka_contigs, params.nanopore_fastq)
    dastool_results = DAStool_nanopore(binning_results.concoct_bins_tsv, binning_results.maxbin2_bins_tsv, semibin2_results.semibin2_bins_tsv, medaka_results.medaka_contigs)
    checkM_results = CheckM_nanopore(dastool_results.dastool_output_bins)
    gtdbtk_results = GTDBTK_nanopore(dastool_results.dastool_output_bins, params.gtdbtk_db)
    prepend_results = PrependBinNames_nanopore(dastool_results.dastool_output_bins)
    vsearch_results = Vsearch_nanopore(prepend_results.merged_bins_fa,params.utax_reference_db)
    AMRfinder_results = AMRfinder_nanopore(prepend_results.merged_bins_fa)
    deeparg_results = Deep_Arg_nanopore(prepend_results.merged_bins_fa)
    plasmidFinder_results = PlasmidFinder_nanopore(prepend_results.merged_bins_fa)
    virsorter2_results = VirSorter2_nanopore(prepend_results.merged_bins_fa,params.virsorter2_db)
    
    }
