#!/bin/sh                                                                                                                                                                     
####################################################################                                                                                                           
#####Adapted form Song Chen's pipeline from Kun Zhang's lab#########                                                                                                           
####################################################################                                                                                                           

##change directory paths as necessary
#Where the jar files are located 
dropseq_dir=/mnt/silencer1/home/amraman/Sebastian/Dropseq_intron
data_dir=/location/of/your/data
#where picard tools is located
picard_dir=/mnt/silencer1/home/amraman/Sebastian/Dropseq_intron/jar

#where STAR aligner is located (used for mapping)
star_dir=/mnt/silencer1/home/amraman/STAR-2.5.2b/bin/Linux_x86_64/

#change to mouse or human or other species 
ref_dir=humGenomeIndex_99/gencode_27  

# this plot tells you about the library quality  (tells you what fraction of read x number of cells contain)
kneeplot_dir=/mnt/silencer1/home/amraman/Sebastian/Sebastian_dropseq                                                          

### When you run this file you can enter ./file samplename species exononly cells 
while getopts ":n:s:c:" options; do
    case $options in
        n ) sampleName=$OPTARG;; # what you want to call the samples (the prefix before the fastq.gz file) 
        s ) species=$OPTARG;;    # m or h (for mouse/human --you can always introduce other options)
        c ) cells=$OPTARG;;      # (# of barcodes you wish to use in your final output) 
    esac
done
shift $(($OPTIND - 1))


################################################################################################                                                                                                                 
#You can also hardcode these options if you want
#sampleName='samplename'  #change sample name  (you can call it whatever you want)                                                                                                    
#species='h'   # change species name  (either m or h)                                                                                                                           
#cells=5000   #change based on # of cells used in the experiment                                                                        
################################################################################################ 

read1=$data_dir/$sampleName.fastq.gz
read2=$data_dir$sampleName.fastq.gz
#################################################################################################
mkdir -p Reports
mkdir -p Tmp

#extract cell barcodes and UMI from read1, convert pair-end fastq files to bam file, and generate fastq file for alignment                                                    \
                                                                                                                                                                               
java -Xmx4g -jar $picard_dir/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 SAMPLE_NAME=$sampleName OUTPUT=/dev/stdout TMP_DIR='pqd'/Tmp |\
java -Xmx4g -jar $picard_dir/picard.jar SortSam I=/dev/stdin O=$sampleName'.unaligned.sorted.bam' SORT_ORDER=queryname TMP_DIR=`pwd`/Tmp

#checking the cell barcodes for low quality bases (marking all cells with more than 1 base < QC threshold) 
$dropseq_dir/TagBamWithReadSequenceExtended I=$sampleName'.unaligned.sorted.bam' O=/dev/stdout SUMMARY=`pwd`/Reports/$sampleName'.cell_tag_report.txt' BASE_RANGE=1-14 BASE_QU\
ALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \

#checking the umi barcodes for low quality bases (marking all cells with more than 1 base < QC threshold) 
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=`pwd`/Reports/$sampleName'.molecule_tag_report.txt' BASE_RANGE=15-24 BASE_QUALITY=10 BARCODED_R\
EAD=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \

#Rejecting the XQ tags (low quality bases) 
$dropseq_dir/FilterBAM TAG_REJECT=XQ I=/dev/stdin O=/dev/stdout | \


$dropseq_dir/TrimStartingSequence INPUT=/dev/stdin OUTPUT=/dev/stdout OUTPUT_SUMMARY=`pwd`/Reports/$sampleName'.adapter_trimming_report.txt' SEQUENCE='TCGTCGGCAGCGTCAGATGTGTA\
TAAGAGACAG' MISMATCHES=0 NUM_BASES=5 | \

#trimming the poly A tail (since we are only counting transcripts from the 3' end) 
$dropseq_dir/PolyATrimmer INPUT=/dev/stdin OUTPUT=/dev/stdout MISMATCHES=0 NUM_BASES=6 | \
tee $sampleName'.unaligned.tagged.bam' | \
java -Xmx8g -jar $picard_dir/picard.jar SamToFastq INPUT=/dev/stdin FASTQ=`pwd`/$sampleName'.unaligned.tagged.fastq'

#########MAPPING INSTRUCTIONS ###########################################
### Only change this is you are not using mouse or human genome 
if [ $species = 'h' ]; then
        refSTAR=$ref_dir
        refFasta=GRCh38.primary_assembly.genome.fa
                refGTF=gencode.v27.primary_assembly.annotation.gtf
else
        refSTAR=../musGenomeIndex
        refFasta=mm10_real.fa
                refGTF=gencode.vM8.primary_assembly.annotation.gtf
fi

############The actual mapping  ##################                                                                                                                                                                         
$star_dir/STAR --runThreadN 16 --genomeDir $ref_dir --readFilesIn `pwd`/$sampleName'.unaligned.tagged.fastq' --outFileNamePrefix $sampleName. --outSAMunmapped Within

#merge aligned sam file with cell barcode/UMI tagged bam file, correct barcode synthesis error, and generate digital expression matrix                                        \
####Instructions for tagging or not tagging reads mapping to introns                                                                                                                                                                               

###################
mkdir -p DGE                                                                                                                                                                                 
###################                                                                                                                                                                                 

## take the mapped reads and map them to exon                                                                                                                                                \                                                                                                                                                                                              
java -Xmx4g -jar $picard_dir/picard.jar SortSam I=$sampleName'.Aligned.out.sam' O=/dev/stdout SO=queryname TMP_DIR=/mnt/thumper/home/amraman/Tmp |\
java -Xmx4g -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout PAIRED_R\
UN=false |\
$dropseq_dir/TagReadWithGeneExon I=/dev/stdin O=$sampleName'.aligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE

## take the unmapped reads and map them to gene                                                                                                                                                                                                 
java -Xmx4g -jar $picard_dir/picard.jar FastqToSam F1=$sampleName'.Unmapped.out.mate1' O=$sample_name'.unmapped.sam' SAMPLE_NAME=$sample_name
$star_dir/STAR --runThreadN 16 --genomeDir $ref_dir --readFilesIn `pwd`/$sampleName'.Unmapped.out.mate1' --outFileNamePrefix $sampleName'.map2' --outSAMunmapped Within                      \
java -Xmx4g -jar $picard_dir/picard.jar SortSam I=$sampleName'.map2.Aligned.out.sam' O=/dev/stdout SO=queryname TMP_DIR=/mnt/thumper/home/amraman/Tmp |\
java -Xmx4g -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout PAIRED_R\
UN=false |\
$dropseq_dir/TagReadWithGene I=/dev/stdin O=$sampleName'.intronaligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE

######### Merge reads that mapped just to exons and those that mapped to introns 
java -Xmx4g -jar $picard_dir/picard.jar MergeSamFiles I=$sampleName'.intronaligned.gene.bam' I=$sampleName'.aligned.gene.bam' O=$sampleName'.aligned_both.gene.bam'
$dropseq_dir/DetectBeadSynthesisErrors I=$sampleName'.aligned_both.gene.bam' O=$sampleName'.aligned.clean.bam' OUTPUT_STATS=`pwd`/Reports/$sampleName'.synthesis_stats.txt' SUMMARY=`pwd`/Rep\
orts/$sampleName'.synthesis_stats_summary.txt' NUM_BARCODES=$cells PRIMER_SEQUENCE=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG

#Make a histogram of barcodes 
$dropseq_dir/BAMTagHistogram I=$sampleName'.aligned_both.gene.bam' O=`pwd`/DGE/$sampleName'.intron_fixreadsByBarcode.txt.gz' TAG=XC

#calaulate the expression matrix 
$dropseq_dir/DigitalExpression I=$sampleName'.aligned_both.gene.bam' O=`pwd`/DGE/$sampleName'.intron_fix_counts.tsv' SUMMARY=`pwd`/Reports/$sampleName'.inton_fix.count_summary.txt' NUM_CORE_\
BARCODES=$cells EDIT_DISTANCE=1


##########################################                                                                                                                                     
mv $sampleName'.Log.out' Reports/
mv $sampleName'.Log.progress.out' Reports/
mv $sampleName'.SJ.out.tab' Reports/
mv $sampleName'.Log.final.out' Reports/

mv $sampleName'.unaligned.sorted.bam' Tmp/
mv $sampleName'.aligned.gene.bam' Tmp/
mv $sampleName'.aligned.clean.bam' Tmp/
mv $sampleName'.unaligned.tagged.fastq' Tmp/
mv $sampleName'.unaligned.tagged.bam' Tmp/
mv $sampleName'.Aligned.out.sam' Tmp/

cd DGE
## For better estimating how many barcodes (cells)  should be analyzed
$kneeplot_dir/KneePlot.R $sampleName $((cells*5))
