/*
 * File: icShape_original.nd
 * Created: Wednesday, 2nd December 2020 1:56:21 pm
 * Author: Zhang Tong (zhangtong516@gmail.com)
 * -----
 * Last Modified: Wednesday, 2nd December 2020 1:56:21 pm
 * Modified By: Zhang Tong (zhangtong516@gmail.com)
 * -----
 * Copyright (c) 2020 GIS
 *
 */

params.help = false
def helpMessage() {
    log.info"""
    The typical command for running the pipeline is as follows:
    nextflow run ${baseDir}/icshape_RT.nf --sampleinfo <> --designinfo <> --batchName <>
    example:
        cd WORKING_DIR
        nextflow run ${baseDir}/icshape_RT.nf --batchName batch1 --sampleInfo RHN1259.sample.tsv

    By default, the pipeline will write the Output files into "results/" folder. It can be changed if user specify another directory.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

// Run ID
if (params.containsKey('batchName') && params.batchName){
    params.outDir = "results/icSHAPE_pipe/" + params.batchName
} else {
    exit 1, "[Pipeline error] Please specify your design info file using `--batchName`!\n"
}

// sample info file
assert (params.sampleInfo != true) : "[Pipeline error] Please specify your sample info file using `--sampleInfo`!\n"

ch_reads = process_sampleinfo(file(params.sampleInfo), params.batchName, "reads")

process collapse_reads {
    tag "${name}"
    publishDir "${params.outDir}" + "/01_collapsedFastq/", mode: 'copy'

    input:
    set val(name), file(fastq_gz) from ch_reads

    output:
    set val(name), file("*.collapsed.fq.gz") into collapsed_reads
    set val(name), file("*.collapsed.read_counts.txt") into collapsed_read_count

    script:
    """
    bash ${baseDir}/scripts/collapse_reads.sh ${fastq_gz} ${name}.collapsed.fq.gz ${task.cpus} ${params.trimLength} ${params.minLength}
    zcat ${fastq_gz} | wc -l | awk '{print "raw_reads\t" \$1/4}' > ${name}.collapsed.read_counts.txt
    zcat  ${name}.collapsed.fq.gz | wc -l | awk '{print "collapsed_reads\t"\$1/4}' >> ${name}.collapsed.read_counts.txt
    """
}

process trim_adaptor {
    tag "${name}"
    publishDir "${params.outDir}" + "/02_trimmedFastq/", mode: 'copy'
    input:
    set val(name), file(collapsed_reads_fastq_gz) from collapsed_reads

    output:
    set val(name), file("*.collapsed.trimmed.fastq.gz") into trimmed_reads
    set val(name), file("*.trimmed.read_counts.txt") into trimmed_read_count
    file("*.read.trim.log.sum")
    file("*.trimAdaptor.log")

    script:
    """
    cutadapt -j ${task.cpus} -n 2 -a ${params.adp1} -g ${params.adp2} -m ${params.minLength} \
        ${collapsed_reads_fastq_gz}  -o ${name}.collapsed.trimmed.fastq.gz  > ${name}.trimAdaptor.log

    perl ${baseDir}/scripts/parse_trim_log.pl ${name}.trimAdaptor.log  > ${name}.read.trim.log.sum

    zcat ${name}.collapsed.trimmed.fastq.gz | wc -l | awk '{print "trimmed_reads\t"\$1/4}' > ${name}.trimmed.read_counts.txt
    """
}

process mapping_rRNA {
    tag "${name}"
    publishDir "${params.outDir}" + "/03_rRNA/", mode: 'copy'

    input:
    set val(name), file(trimmed_reads_fastq_gz) from trimmed_reads

    output:
    set val(name), file("*.collapsed.trimmed.remove_rRNA.fastq") into removed_rRNA_fq
    set val(name), file("*.collapsed.trimmed.map_rRNA.sam") into mapped_rRNA_sam
    file("*.collapsed.trimmed.map_rRNA.log")


    script:
    """
    icSHAPE-pipe cleanFq -i ${trimmed_reads_fastq_gz} -o ${name}.collapsed.trimmed.remove_rRNA.fastq \
        -x ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA \
        -p ${task.cpus} --mode Local \
        --sam ${name}.collapsed.trimmed.map_rRNA.sam 2> \
        ${name}.collapsed.trimmed.map_rRNA.log

    """
}

process mapping_smallRNA {
    tag "${name}"
    publishDir "${params.outDir}" + "/04_smallRNA/", mode: 'copy'

    input:
    set val(name), file(rm_rRNA_fastq) from removed_rRNA_fq

    output:
    set val(name), file("*.collapsed.trimmed.remove_rRNA_smallRNA.fastq") into removed_smallRNA_fq
    set val(name), file("*.collapsed.trimmed.map_smallRNA.sam") into mapped_smallRNA_sam
    file("*.collapsed.trimmed.map_smallRNA.log")


    script:
    """
    icSHAPE-pipe cleanFq -i ${rm_rRNA_fastq} -o ${name}.collapsed.trimmed.remove_rRNA_smallRNA.fastq \
        -x ${params.genome_dir}/smallRNA/smallRNA \
        -p ${task.cpus} --mode Local \
        --sam ${name}.collapsed.trimmed.map_smallRNA.sam 2> \
        ${name}.collapsed.trimmed.map_smallRNA.log
    """
}


process star_mapping_genome {
    tag "${name}"
    publishDir "${params.outDir}" + "/05_mappedResult/", mode: 'copy'

    input:
    set val(name), file(rm_smallRNA_fastq) from removed_smallRNA_fq

    output:
    set val(name), file("*.map_genome.sorted.bam") into mapped_bam_for_fpkm
    set val(name), file("*.map_genome.sorted.bam") into mapped_bam_for_sam2tab
    set val(name), file("*.map_genome.sorted.bam") into mapped_bam_for_coverage
    set val(name), file("*.map_genome.unsorted.bam") into unsorted_bam
    set val(name), file("*.map_genome.Log.progress.out") into log_progress
    set val(name), file("*.map_genome.Log.final.out") into log_final
    set val(name), file("*.map_genome.Log.out") into log_out
    set val(name), file("*.map_genome.Log.std.out") into log_std


    script:
    """
    STAR --readFilesIn ${rm_smallRNA_fastq} \
        --outFileNamePrefix ${name}.map_genome. \
        --genomeDir ${params.genome_dir}/star \
        --runThreadN ${task.cpus} \
        --genomeLoad NoSharedMemory \
        --runMode alignReads \
        --outSAMtype BAM Unsorted \
        --outSAMmultNmax 1 \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 2 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMstrandField intronMotif \
        --outSJfilterOverhangMin 30 12 12 12 \
        --alignEndsType EndToEnd \
        --outSAMattributes All \
        --outSAMunmapped Within \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --outStd BAM_Unsorted > ${name}.map_genome.unsorted.bam

    samtools view -h ${name}.map_genome.unsorted.bam |\
        awk '\$0~/^@/{print \$0}\$0!~/^@/{for(i=12;i<NF;i++){if(substr(\$i,1,4)=="MD:Z"){if(and(16,\$2)==0){ if(\$i!~/^MD:Z:0/ ) print \$0; }else{if(\$i!~/^MD:Z:.*0\$/) print \$0; }}}}' |\
        samtools view --threads ${task.cpus} -bh - |\
        samtools sort - -m 2G --threads ${task.cpus} \
            -o ${name}.map_genome.sorted.bam
    """
}


process estimate_rpkm{
    tag "${name}"
    publishDir "${params.outDir}" + "/06_fpkm/", mode: 'copy'

    input:
    set val(name), file(mapped_bam) from mapped_bam_for_fpkm

    output:
    set val(name), file("*.gene.txt") into gene_rpkm
    set val(name), file("*.txn.txt") into txn_rpkm
    set val(name), file("*.gene.txt.summary") into gene_rpkm_summary
    set val(name), file("*.txn.txt.summary") into txn_rpkm_summary
    set val(name), file("*.isoforms.fpkm_tracking") into isoform_fpkm
    set val(name), file("*.isoforms.fpkm_tracking") into isoform_fpkm2
    set val(name), file("*.genes.fpkm_tracking") into gene_fpkm
    set val(name), file("*.skipped.gtf") into skipped_gtf
    set val(name), file("*.transcripts.gtf") into transcripts_gtf

    script:
    """
    icSHAPE-pipe calcFPKM -i ${mapped_bam} -o ${name} \
        -G ${params.genome_dir}/Homo_sapiens.GRCh38.100.gtf -p ${task.cpus}

    mv ./${name}/isoforms.fpkm_tracking ${name}.isoforms.fpkm_tracking
    mv ./${name}/genes.fpkm_tracking ${name}.genes.fpkm_tracking
    mv ./${name}/skipped.gtf ${name}.skipped.gtf
    mv ./${name}/transcripts.gtf ${name}.transcripts.gtf

    featureCounts -T ${task.cpus} \
        -a ${params.genome_dir}/Homo_sapiens.GRCh38.100.gtf \
        -g transcript_id \
        -o ${name}.txn.txt \
        ${mapped_bam}

    featureCounts -T ${task.cpus} \
        -a ${params.genome_dir}/Homo_sapiens.GRCh38.100.gtf \
        -g gene_id \
        -o ${name}.gene.txt \
        ${mapped_bam}
    """
}

mapped_bams = mapped_bam_for_sam2tab.join(mapped_rRNA_sam).join(mapped_smallRNA_sam)

process sam2tab {
    tag "${name}"
    publishDir "${params.outDir}" + "/07_sam2tab/", mode: 'copy'

    input:
    set val(name), file(mapped_bam), file(mapped_rRNA_sam), file(mapped_smallRNA_sam) from mapped_bams

    output:
    set val(name), file("*.genome.tab"), file("*.rRNA.tab"), file("*.smallRNA.tab") into genome_sam2tab1
    set val(name), file("*.genome.tab"), file("*.rRNA.tab"), file("*.smallRNA.tab") into genome_sam2tab2
    set val(name), file("*.genome.tab"), file("*.rRNA.tab"), file("*.smallRNA.tab") into genome_sam2tab3
    set val(name), file("*.genome.tab"), file("*.rRNA.tab"), file("*.smallRNA.tab") into genome_sam2tab4

    script:
    """
    icSHAPE-pipe sam2tab -in ${mapped_bam} -out ${name}.genome.tab
    icSHAPE-pipe sam2tab -in ${mapped_rRNA_sam} -out ${name}.rRNA.tab
    icSHAPE-pipe sam2tab -in ${mapped_smallRNA_sam} -out ${name}.smallRNA.tab
    """
}

ch_for_genome_SHAPE_single_DMSO = genome_sam2tab1.filter{it[1] =~/__DMSO__/}.map{
    [process_sample_name(it[0], "single"), it[1], it[2], it[3]] }
ch_for_genome_SHAPE_single_NAIN3 = genome_sam2tab2.filter{it[1] =~/__NAIN3__/}.map{
    [process_sample_name(it[0], "single"), it[1], it[2], it[3]] }

ch_for_genome_SHAPE_single =  ch_for_genome_SHAPE_single_DMSO.join(ch_for_genome_SHAPE_single_NAIN3)


process tab2gTab_single {
    tag "${name}"
    publishDir "${params.outDir}" + "/08_genomeSHAPE_single/", mode: 'copy'

    input:
    set val(name), file(dmso_genome_tab), file(dmso_rRNA_tab), file(dmso_smallRNA_tab), file(nain3_genome_tab), file(nain3_rRNA_tab), file(nain3_smallRNA_tab) from ch_for_genome_SHAPE_single

    output:
    set val(name), file("*.genome.single.gTab"), file("*.rRNA.single.gTab"), file("*.smallRNA.single.gTab") into single_gTab
    file("*.genome.single.gTab.param.log")
    file("*.rRNA.single.gTab.param.log")
    file("*.smallRNA.single.gTab.param.log")

    script:
    """
    icSHAPE-pipe calcSHAPE -D ${dmso_genome_tab} -N ${nain3_genome_tab} \
        -size ${params.genome_dir}/star/chrNameLength.txt \
        -ijf ${params.genome_dir}/star/sjdbList.fromGTF.out.tab \
        -genome ${params.genome_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -bases A,T,C,G \
        -out ${name}.genome.single.gTab

    icSHAPE-pipe calcSHAPE -D ${dmso_rRNA_tab} -N ${nain3_rRNA_tab} \
        -size ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA.len \
        -genome ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA.fa \
        -bases A,T,C,G \
        -non-sliding \
        -out ${name}.rRNA.single.gTab

    icSHAPE-pipe calcSHAPE -D ${dmso_smallRNA_tab} -N ${nain3_smallRNA_tab} \
        -size ${params.genome_dir}/smallRNA/smallRNA.len \
        -genome ${params.genome_dir}/smallRNA/smallRNA.fa \
        -bases A,T,C,G \
        -non-sliding \
        -out ${name}.smallRNA.single.gTab
    """
}

dmso_isoform_fpkm_single = isoform_fpkm2.filter{it[1] =~/__DMSO__/}.map{
    [process_sample_name(it[0], "single"), it[1]] }
single_gTab_rpkm = single_gTab.join(dmso_isoform_fpkm_single)

process gTab2RTBD_single {
    // slow process need to iterate all genome
    tag "${name}"
    publishDir "${params.outDir}" + "/08_genomeSHAPE_single/", mode: 'copy'

    input:
    set val(name), file(single_genome_gTab), file(single_rRNA_gTab), file(single_smallRNA_gTab), file(isoform_fpkms) from single_gTab_rpkm

    output:
    set val(name), file("*.genome.single.RTBD") into single_RTBD
    set val(name), file("*.trans.single.shape") into single_shape
    

    script:
    """
    icSHAPE-pipe genRTBDToTransRTBD -i ${single_genome_gTab} \
        -g ${params.genome_dir}/GTF/Anno.genomeCoor.bed \
        -p ${task.cpus}  -c 5,6,7,8 \
        -o ${name}.genome.single.RTBD

    icSHAPE-pipe genSHAPEToTransSHAPE -i ${single_genome_gTab} \
        -o ${name}.trans.single.shape \
        -g ${params.genome_dir}/GTF/Anno.genomeCoor.bed \
        -p ${task.cpus} \
        -r ${isoform_fpkms} \
        -c ${params.minCov} -T ${params.minRT} -M ${params.minRPKM}
    
    icSHAPE-pipe genSHAPEToTransSHAPE -i ${single_rRNA_gTab} \
        -o ${name}.trans.single.shape \
        -s ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA.len \
        -p ${task.cpus} \
        --app -c ${params.minCov} -T ${params.minRT} -M ${params.minRPKM}

    icSHAPE-pipe genSHAPEToTransSHAPE -i ${single_smallRNA_gTab} \
        -o ${name}.trans.single.shape \
        -s ${params.genome_dir}/smallRNA/smallRNA.len \
        -p ${task.cpus} \
        --app -c ${params.minCov} -T ${params.minRT} -M ${params.minRPKM}

    """
}

process coverage_rt_count {
    tag "${name}"
    publishDir "${params.outDir}" + "/12_coverageCount_Single/", mode: 'copy'

    input:
    set val(name), file(single_rtbd_file) from single_RTBD
    output:
    file("*.rtbd.txt") into single_coverage
    file("*.rtbd.txt.log")

    script:
    """
    python ${baseDir}/scripts/transform_RTBD_file.py  ${single_rtbd_file}  ${name}.rtbd.txt
    """
}


ch_for_merge_coverage = single_coverage.collect().map{ wrap_items(it[0], sep=":") }

process merge_coverage {
    tag "${name}"
    publishDir "${params.outDir}" + "/12_coverageCount_Single/", mode: 'copy'

    input:
    val(input_files) from ch_for_merge_coverage

    output:
    file("all.rt.merged.wide.csv.gz")
    file("all.cov.merged.wide.csv.gz")

    script:
    """
    Rscript ${baseDir}/scripts/merge_RTBD_txt.R  ${input_files}
    """
}


ch_for_genome_SHAPE_DMSO = genome_sam2tab3.filter{it[1] =~/__DMSO__/}.map{
    [process_sample_name(it[0], "sampleName"), it[1], it[2], it[3]] }. groupTuple().map{
        [ it[0], wrap_items(it[1],sep=","), wrap_items(it[2],sep=","), wrap_items(it[3],sep=",") ]
    }
ch_for_genome_SHAPE_NAIN3 = genome_sam2tab4.filter{it[1] =~/__NAIN3__/}.map{
    [process_sample_name(it[0], "sampleName"), it[1], it[2], it[3]] }. groupTuple().map{
        [ it[0], wrap_items(it[1],sep=","), wrap_items(it[2],sep=","), wrap_items(it[3],sep=",") ]
    }
ch_for_genome_SHAPE =  ch_for_genome_SHAPE_DMSO.join(ch_for_genome_SHAPE_NAIN3)


process tab2gTab_merged {
    tag "${name}"
    publishDir "${params.outDir}" + "/09_genomeSHAPE/", mode: 'copy'

    input:
    set val(name), val(dmso_genome_tab), val(dmso_rRNA_tab), val(dmso_smallRNA_tab), val(nain3_genome_tab), val(nain3_rRNA_tab), val(nain3_smallRNA_tab) from ch_for_genome_SHAPE

    output:
    set val(name), file("*.genome.gTab"), file("*.rRNA.gTab"), file("*.smallRNA.gTab") into merged_gTab

    script:
    """
    icSHAPE-pipe calcSHAPE -D ${dmso_genome_tab} -N ${nain3_genome_tab} \
        -size ${params.genome_dir}/star/chrNameLength.txt \
        -ijf ${params.genome_dir}/star/sjdbList.fromGTF.out.tab \
        -genome ${params.genome_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -bases A,T,C,G \
        -out ${name}.genome.gTab

    icSHAPE-pipe calcSHAPE -D ${dmso_rRNA_tab} -N ${nain3_rRNA_tab} \
        -size ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA.len \
        -genome ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA.fa \
        -bases A,T,C,G \
        -non-sliding \
        -out ${name}.rRNA.gTab

    icSHAPE-pipe calcSHAPE -D ${dmso_smallRNA_tab} -N ${nain3_smallRNA_tab} \
        -size ${params.genome_dir}/smallRNA/smallRNA.len \
        -genome ${params.genome_dir}/smallRNA/smallRNA.fa \
        -bases A,T,C,G \
        -non-sliding \
        -out ${name}.smallRNA.gTab
    """
}

dmso_isoform_fpkm = isoform_fpkm.filter{it[1] =~/__DMSO__/}.map{
    [process_sample_name(it[0], "sampleName"), it[1] ]}.groupTuple().map{
        [ it[0], wrap_items(it[1],sep=",") ]
    }

ch_for_SHAPE = merged_gTab.join(dmso_isoform_fpkm)

process calculate_shape {
    tag "${name}"
    publishDir "${params.outDir}" + "/10_SHAPE/", mode: 'copy'

    input:
    set val(name), file(genome_gTab), file(rRNA_gTab), file(smallRNA_gTab), val(isoform_fpkms) from ch_for_SHAPE

    output:
    set val(name), file("*.final.shape") into final_shape

    script:
    """
    icSHAPE-pipe genSHAPEToTransSHAPE -i ${genome_gTab} \
        -o ${name}.final.shape \
        -g ${params.genome_dir}/GTF/Anno.genomeCoor.bed \
        -p ${task.cpus} \
        -r ${isoform_fpkms} \
        -c ${params.minCov} -T ${params.minRT} -M ${params.minRPKM}

    icSHAPE-pipe genSHAPEToTransSHAPE -i ${rRNA_gTab} \
        -o ${name}.final.shape \
        -s ${params.genome_dir}/rRNA/human_rRNA_tRNA_mtRNA.len \
        -p ${task.cpus} \
        --app -c ${params.minCov} -T ${params.minRT} -M ${params.minRPKM}

    icSHAPE-pipe genSHAPEToTransSHAPE -i ${smallRNA_gTab} \
        -o ${name}.final.shape \
        -s ${params.genome_dir}/smallRNA/smallRNA.len \
        -p ${task.cpus} \
        --app -c ${params.minCov} -T ${params.minRT} -M ${params.minRPKM}
    """
}

process generate_bedgraph {
    tag "${name}"
    publishDir "${params.outDir}" + "/11_bedgraph/", mode: 'copy'

    input:
    set val(name), file(genome_shape) from final_shape
    output:
    set val(name), file("*.transcriptome.bedgraph.gz") into transcritome_bedgraph

    script:
    """
    python ${baseDir}/scripts/shape_to_bedGraph.py \
        -i ${genome_shape} \
        -o ${name}.transcriptome.bedgraph.gz \
        -n ${name}
    """
}

// extra functions to parse my sample.tsv input file
def wrap_items(input_files, sep=";") {
    // example input: [file1, file2]
    // output: file1+sep+file2
    def result =  input_files instanceof Path ? input_files.toString() : (input_files as List).join(sep)
    return result.toString()
}

def process_sample_name(condition, mode="single", sep="__") {
    // Process sample Names: example: D0__DMSO__rep1
    // mode == single: D0__DMSO__rep1  --> D0__rep1
    // mode == merge:  D0__DMSO__rep1  --> D0__DMSO
    def p = condition.toString().split(sep)
    if(mode == "single") { result = p[0] + sep + p[2] }
    if(mode == "merge") { result = p[0] + sep + p[1] }
    if(mode == "sampleName") { result = p[0]}
    return(result)

}

def process_sampleinfo(tsvFile, batchName, mode) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def sample          = row.sampleName // D0
            def libID           = row.libID // example: RHN1259Lib1
            def seqBatch        = row.seqBatch  //example: batch1
            def compBatch       = row.comparisonBatch //example: batch2
            def libType         = row.treatment // NAIN3 or DMSO
            def replicate       = row.rep // rep1
            def fastq           = row.fq // file with full path
            def condition       = row.sampleName + "__" + row.treatment
            def unique_name     = condition + "__" + replicate

            def genome_tab      = params.outDir + "/07_sam2tab/" + unique_name + ".genome.tab"
            def rRNA_tab        = params.outDir + "/07_sam2tab/" + unique_name + ".rRNA.tab"
            def smallRNA_tab    = params.outDir + "/07_sam2tab/" + unique_name + ".smallRNA.tab"

            def nrt             = params.outDir  + condition + ".merged.nrt"

            if (mode == "reads" && compBatch == batchName) return [ unique_name, file(fastq) ]
            if (mode == "sam2tab" && compBatch == batchName) return [ unique_name, file(genome_tab), file(rRNA_tab), file(smallRNA_tab) ]
            if (mode == "normalizedRT" && compBatch == batchName) return [ condition, file(nrt) ]
        }
        .unique()
    }
