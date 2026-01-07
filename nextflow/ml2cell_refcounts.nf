// end can be  for single end or '\-p' on the command line for paired end

params.endparameter = '-p'

ch_refdata=Channel.fromPath(params.refdata).collect()
ch_end=Channel
.value(params.endparameter)
ch_bams = Channel 
.fromPath(params.bamlocation)

ch_annot=ch_refdata


process CUSTOMFEATURECOUNTS{
conda "/home/rstudio/ml2cell_code/nextflow/environments/subread.yml"
container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0' 

    input:
    //val(meta) //commenting out all the meta references.
    val(end)
    path(bams)
    path(annotation)

    output:
    path("*featureCounts.txt")        , emit: counts
    path("*featureCounts.txt.summary"), emit: summary
    path "versions.yml"                                , emit: versions
    path ("*onlycounts.txt")                              , emit: countcol 


script:
    def args = task.ext.args ?: ''
    def prefix = 'testfcounts' // task.ext.prefix ?: "${meta.id}"
    //def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    /*
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    
    */
    """
    featureCounts \\
        $args \\
        $end \\
        -T $task.cpus \\
        -a $annotation \\
        -s $strandedness \\
        -F 'SAF' \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}
        
    

    cut -f 7 ${prefix}.featureCounts.txt > ${bams}.onlycounts.txt
    


     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """	
}

process MERGECOUNTS{
container  'ghcr.io/cicirello/gnu-on-alpine:latest'

input:
path files
path ref

output:

path "allmergedcounts.txt" ,  emit: mergedcounts
path "allmergedcounts_withreference.txt", emit:refandcounts

publishDir params.outdir, mode:'copy'

script:

"""
# count the total number of rows we are dealing with. probably best to obtain this from the root peaks file in the future. 
#num_rows=\$(wc -l < $files[0])
# Create an array quoted of file names
#files=(${files.collect { "\"$it\"" }.join(' ')})
# Create an array of unquoted file paths



paste ${files.join(' ')} > allmergedcounts.txt

#tail -n +2 allmergedcounts.txt

paste $ref allmergedcounts.txt >allmergedcounts_withreference.txt

"""
}


workflow
{
ch_cfcounts= CUSTOMFEATURECOUNTS(ch_end, ch_bams, ch_refdata)
ch_mcounts= MERGECOUNTS(ch_cfcounts.countcol.collect(), ch_refdata)

}


