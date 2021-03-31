# Assessing-Intratumor-Genetic-Heterogeneity

* Author: Jieting Hu

# Overview of project:

* This is a project based on Whole-Exome Sequencing and the assessment of its reliability on predicting true intratumor genetic heterogeneity. 

* The data and assessment is based on : https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31635-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124718316358%3Fshowall%3Dtrue

* Intratumor genetic heterogeneity (ITGH) is referred to as the cells with different molecular and phenotypical profiles within the same tumor specimen. This study aims to evaluate the reliability of WES analysis of somatic variant detection
 in the context of Intratumor Genetic Heterogeneity. The data is obtained from WES on DNA from biopsies obtained from three anatomically distinct regions of six primary breast cancers (6 × 3 biopsies) and from matched peripheral blood leukocytes.
 
 * Basically, There are 6 surgically resected tumors, each sampled at 3 distinct regions, a,b and c biological replicates. There are also one germline replicate and one technical replicate so a total of 30 samples. 
 
 
# step1: get data for analysis:

* ## make wes_cancer file as working directory
```
mkdir ~/wes_cancer
cd ~/wes_cancer
```
in ~/wes_cancer make biosoft project data 3 directories
biosoft for sotfware
project for intermediate files
data for data
```
mkdir biosoft project data
cd project
```
in project make intermediate files.
```
mkdir -p 0.sra 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/{qualimap,flagstat,stats} 5.gatk/gvcf 6.mutect 7.annotation/{vep,annovar,funcotator,snpeff} 8.cnv/{gatk,cnvkit,gistic,facet} 9.pyclone 10.signature
```

* The software installation is omited. 

* download the annotation file for reference genome(bed)
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
cat CCDS.current.txt | grep  "Public" | perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}'|sort -u |bedtools sort -i |awk '{if($3>$2) print "chr"$0"\t0\t+"}'  > hg38.exon.bed
```

* Using Aspera to download gene files from https://sra-explorer.info/
```
$ cd ~/wes_cancer/project/1.raw_fq
$ cat >fq_download.sh
#!/usr/bin/env bash
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182418/SRR3182418.fastq.gz . && mv SRR3182418.fastq.gz SRR3182418_exome_sequencing_of_case5_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182418/SRR3182418_1.fastq.gz . && mv SRR3182418_1.fastq.gz SRR3182418_exome_sequencing_of_case5_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182418/SRR3182418_2.fastq.gz . && mv SRR3182418_2.fastq.gz SRR3182418_exome_sequencing_of_case5_germline_2.fastq.gz
......
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182446/SRR3182446_2.fastq.gz . && mv SRR3182446_2.fastq.gz SRR3182446_exome_sequencing_of_case5_techrep_2_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182447/SRR3182447.fastq.gz . && mv SRR3182447.fastq.gz SRR3182447_exome_sequencing_of_case5_biorep_C.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182447/SRR3182447_1.fastq.gz . && mv SRR3182447_1.fastq.gz SRR3182447_exome_sequencing_of_case5_biorep_C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182447/SRR3182447_2.fastq.gz . && mv SRR3182447_2.fastq.gz SRR3182447_exome_sequencing_of_case5_biorep_C_2.fastq.gz
^C

nohup  bash  fq_download.sh  &
```

delete unnecessary files
```
find ./*gz -size 1M | xargs rm
 ##change names of the files 
ls *gz | while read id; do mv ${id} ${id:31:100}; done
```

 The output becomes:
 ```
~/wes_cancer/project/1.raw_fq$ ls
case3_germline_1.fastq.gz     
case3_germline_2.fastq.gz     
case4_germline_1.fastq.gz 
...
...
nohup.out		   

```
* Fastqc and Trimming Adapters
```
conda activate wes
cd ~/wes_cancer/project

(wes) ubuntu@VM-0-17-ubuntu:~/wes_cancer/project$ conda activate wes
(wes) ubuntu@VM-0-17-ubuntu:~/wes_cancer/project$ cd ~/wes_cancer/project
(wes) ubuntu@VM-0-17-ubuntu:~/wes_cancer/project$ 
(wes) ubuntu@VM-0-17-ubuntu:~/wes_cancer/project$ cat config | while read id  do fastqc --outdir ./3.qc/raw_qc/ --threads 16 ./1.raw_fq/${id}*.fastq.gz >>       ./3.qc/raw_qc/${id}_fastqc.log 2>&1 
done

multiqc  ./3.qc/raw_qc/*zip  -o ./3.qc/raw_qc/multiqc

cat config | while read id
do
	fastqc --outdir ./3.qc/raw_qc/ --threads 16 ./1.raw_fq/${id}*.fastq.gz >> ./3.qc/raw_qc/${id}_fastqc.log 2>&1 
done 

multiqc  ./3.qc/raw_qc/*zip  -o ./3.qc/raw_qc/multiqc
```
```
# Using trim_galore to trim bad quality fastq and remove adapters. 
## trim_galore.sh
cat config | while read id
do
	fq1=./1.raw_fq/${id}_1.fastq.gz
	fq2=./1.raw_fq/${id}_2.fastq.gz
	trim_galore  --paired -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 8 -o ./2.clean_fq  $fq1  $fq2 >> ./2.clean_fq/${id}_trim.log 2>&1
done

nohup bash trim_galore.sh &
```

```
# The data is so big and I am sharing this HPC with others and some of my samples get deleted bu others to allocate more space. This is just a practice so 
not a big deal to just analyze this much data. 
```
```
cd ~/wes_cancer/data
gunzip Homo_sapiens_assembly38.fasta.gz
time bwa index -a bwtsw -p gatk_hg38 ~/wes_cancer/data/Homo_sapiens_assembly38.fasta
[bwa_index] Pack FASTA... 23.77 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=6434693834, availableWord=464768632
[BWTIncConstructFromPacked] 10 iterations done. 99999994 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 199999994 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 299999994 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 399999994 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 499999994 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 599999994 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 699999994 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 799999994 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 899999994 characters processed.
[BWTIncConstructFromPacked] 100 iterations done. 999999994 characters processed.
[BWTIncConstructFromPacked] 110 iterations done. 1099999994 characters processed.
[BWTIncConstructFromPacked] 120 iterations done. 1199999994 characters processed.
[BWTIncConstructFromPacked] 130 iterations done. 1299999994 characters processed.
[BWTIncConstructFromPacked] 140 iterations done. 1399999994 characters processed.
[BWTIncConstructFromPacked] 150 iterations done. 1499999994 characters processed.
[BWTIncConstructFromPacked] 160 iterations done. 1599999994 characters processed.
[BWTIncConstructFromPacked] 170 iterations done. 1699999994 characters processed.
```

```
[bwt_gen] Finished constructing BWT in 711 iterations.
[bwa_index] 2519.46 seconds elapse.
[bwa_index] Update BWT... 16.50 sec
[bwa_index] Pack forward-only FASTA... 15.06 sec
[bwa_index] Construct SA from BWT and Occ... 1038.62 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -a bwtsw -p gatk_hg38 /home/ubuntu/wes_cancer/data/Homo_sapiens_assembly38.fasta
[main] Real time: 3712.628 sec; CPU: 3613.414 sec

real	61m52.635s
user	59m57.243s
sys	0m16.172s
```

```
mv *_val_2.fq.gz ~/wes_cancer/project/2.clean_fq
mv *_val_1.fq.gz ~/wes_cancer/project/2.clean_fq

## bwa.sh
INDEX=~/trybwa/gatk_hg38/gatk_hg38
cat config | while read id
do 
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}
bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2 |samtools sort  -@ 2 -o $sample.bam -
done 

bash bwa.sh
```
```
ubuntu@VM-0-17-ubuntu:~/wes_cancer/project/2.clean_fq$ bash bwa.sh
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 5 -R @RG\tID:case1_biorep_B_1_val_1.fq.gz\tSM:case1_biorep_B_1_val_1.fq.gz\tLB:WGS\tPL:Illumina /home/ubuntu/trybwa/gatk_hg38/gatk_hg38 case1_biorep_B_2_val_2.fq.gz
[main] Real time: 2.895 sec; CPU: 2.896 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 5 -R @RG\tID:case2_biorep_B_techrep_1_1_val_1.fq.gz\tSM:case2_biorep_B_techrep_1_1_val_1.fq.gz\tLB:WGS\tPL:Illumina /home/ubuntu/trybwa/gatk_hg38/gatk_hg38 case2_biorep_A_2_val_2.fq.gz

...
[main] CMD: bwa mem -t 5 -R @RG\tID:case6_biorep_C_1_val_1.fq.gz\tSM:case6_biorep_C_1_val_1.fq.gz\tLB:WGS\tPL:Illumina /home/ubuntu/trybwa/gatk_hg38/gatk_hg38 case6_biorep_C_2_val_2.fq.gz
[main] Real time: 2.775 sec; CPU: 2.775 sec
ubuntu@VM-0-17-ubuntu:~/wes_cancer/project/2.clean_fq$ ls
case1_biorep_B_1_val_1.fq.gz
case1_biorep_B_1_val_1.fq.gz.bam
case1_biorep_B_2_val_2.fq.gz
case2_biorep_A_2_val_2.fq.gz
case2_biorep_B_techrep_1_1_val_1.fq.gz
case2_biorep_B_techrep_1_1_val_1.fq.gz.bam
case2_biorep_C_2_val_2.fq.gz
case2_germline_1_val_1.fq.gz
case2_germline_1_val_1.fq.gz.bam
case2_germline_2_val_2.fq.gz
case3_biorep_A_1_val_1.fq.gz
case3_biorep_A_1_val_1.fq.gz.bam
case3_biorep_A_2_val_2.fq.gz
case3_biorep_B_1_val_1.fq.gz
case3_biorep_B_1_val_1.fq.gz.bam
case3_biorep_B_2_val_2.fq.gz
case3_biorep_C_techrep_1_1_val_1.fq.gz
case3_biorep_C_techrep_1_1_val_1.fq.gz.bam
case3_germline_2_val_2.fq.gz
case4_biorep_A_1_val_1.fq.gz
case4_biorep_A_1_val_1.fq.gz.bam
case4_germline_2_val_2.fq.gz
case4_techrep_2_1_val_1.fq.gz
case4_techrep_2_1_val_1.fq.gz.bam
case5_biorep_A_1_val_1.fq.gz
case5_biorep_A_1_val_1.fq.gz.bam
case5_biorep_A_2_val_2.fq.gz
case5_biorep_B_techrep_1_1_val_1.fq.gz
case5_biorep_B_techrep_1_1_val_1.fq.gz.bam
case5_biorep_B_techrep_1_2_val_2.fq.gz
case6_biorep_A_techrep_1_1_val_1.fq.gz
case6_biorep_A_techrep_1_1_val_1.fq.gz.bam
case6_biorep_A_techrep_1_2_val_2.fq.gz
case6_biorep_C_1_val_1.fq.gz
case6_biorep_C_1_val_1.fq.gz.bam
case6_biorep_C_2_val_2.fq.gz
case6_germline_1.fastq.gz_trim.log
```

```
mv *bam ~/wes_cancer/project/4.align
ubuntu@VM-0-17-ubuntu:~/wes_cancer/project/2.clean_fq$ cd ~/wes_cancer/project/4.align
ubuntu@VM-0-17-ubuntu:~/wes_cancer/project/4.align$ ls
case1_biorep_B_1_val_1.fq.gz.bam
case2_biorep_B_techrep_1_1_val_1.fq.gz.bam
case2_germline_1_val_1.fq.gz.bam
case3_biorep_A_1_val_1.fq.gz.bam
case3_biorep_B_1_val_1.fq.gz.bam
case3_biorep_C_techrep_1_1_val_1.fq.gz.bam
case4_biorep_A_1_val_1.fq.gz.bam
case4_techrep_2_1_val_1.fq.gz.bam
case5_biorep_A_1_val_1.fq.gz.bam
case5_biorep_B_techrep_1_1_val_1.fq.gz.bam
case6_biorep_A_techrep_1_1_val_1.fq.gz.bam
case6_biorep_C_1_val_1.fq.gz.bam
flagstat
qualimap
stats
```
| col | Name | description | example |
| --- | --- | --- | --- |
| 1 | QNAME | The name of the read | case6_biorep_A_tech |
| 2 | FLAG | binary description | 129 
| 3 | RNAME | Reference sequence NAME, CHR | CHR17 |
| 4 | POS | 1-based leftmost mapping POSition on CHR | 42852401 |
| 5 | MAPQ | MAPping Quality | 60|
| 6 | CIGAR | \*|([0-9]+[MIDNSHPX=])+ M I D | 76M3D |  
| 7 | RNEXT | Reference name of the mate/next read | Chr11 |
| 8 | PNEXT |  observed Template LENgth | 96427099 |
| 9 | TLEN  | observed Template LENgth | 0 |
| 10 | SEQ | The second line of the file seq | CGTTC...TCT |
| 11 | QUAL | ASCII of Phred-scaled base QUALity+33 | @@@?D..CA


```
#extract a small bam for IGV visualizaiton. 
ubuntu@VM-0-17-ubuntu:~/wes_cancer/project/4.align$ samtools index case1_biorep_B_1_val_1.fq.gz.bam
samtools view -h case1_biorep_B_1_val_1.fq.gz.bam chr17 | samtools view -Sb - > small.bam
ubuntu@VM-0-17-ubuntu:~/wes_cancer/project/4.align$ samtools index small.bam
```

```
## stats.sh
cat config | while read id
do
	bam=./4.align/${id}.bam
	samtools stats -@ 16 --reference ~/wes_cancer/data/Homo_sapiens_assembly38.fasta ${bam} > ./4.align/stats/${id}.stat

	plot-bamstats -p ./4.align/stats/${id} ./4.align/stats/${id}.stat
done

#qualimap
cat config | while read id
do
	qualimap bamqc --java-mem-size=10G -gff ~/wes_cancer/data/hg38.exon.bed -nr 100000 -nw 500 -nt 16 -bam ./4.align/${id}.bam -outdir ./4.align/qualimap/${id}
done

#GATK
dbsnp_146.hg38.vcf.gz 
dbsnp_146.hg38.vcf.gz.tbi 
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi 
Homo_sapiens_assembly38.fasta
Homo_sapiens_assembly38.fasta.gz 
Homo_sapiens_assembly38.fasta.fai 
Homo_sapiens_assembly38.dict 
1000G_phase1.snps.high_confidence.hg38.vcf.gz 
1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi




#markduplicate.
GATK=~/wes_cancer/biosoft/gatk-4.1.4.1/gatk

cat config  | while read id
do
	BAM=./4.align/${id}.bam
	if [ ! -f ./5.gatk/ok.${id}_marked.status ]
	then
		echo "start MarkDuplicates for ${id}" `date`
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
		-I ${BAM} \
		--REMOVE_DUPLICATES=true \
		-O ./5.gatk/${id}_marked.bam \
		-M ./5.gatk/${id}.metrics \
		1>./5.gatk/${id}_log.mark 2>&1 
		
		if [ $? -eq 0 ]
		then
			touch ./5.gatk/ok.${id}_marked.status
		fi
		echo "end MarkDuplicates for ${id}" `date`
		samtools index -@ 16 -m 4G -b ./5.gatk/${id}_marked.bam ./5.gatk/${id}_marked.bai
	fi
done

samtools view case1_biorep_B_1_1.fq.gz.bam | wc -l
93168816

samtools view case1_biorep_B_1_1.fq.gz_marked.bam | wc -l
93168816



#Base quality score recalibration
GATK=~/wes_cancer/biosoft/gatk-4.1.4.1/gatk
snp=~/wes_cancer/data/dbsnp_146.hg38.vcf.gz
indel=~/wes_cancer/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
cat config  | while read id
do
	if [ ! -f ./5.gatk/${id}_bqsr.bam ]
	then
		echo "start BQSR for ${id}" `date`
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
		-R $ref  \
		-I ./5.gatk/${id}_marked.bam  \
		--known-sites ${snp} \
		--known-sites ${indel} \
		-O ./5.gatk/${id}_recal.table \
		1>./5.gatk/${id}_log.recal 2>&1 
		
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  ApplyBQSR \
		-R $ref  \
		-I ./5.gatk/${id}_marked.bam  \
		-bqsr ./5.gatk/${id}_recal.table \
		-O ./5.gatk/${id}_bqsr.bam \
		1>./5.gatk/${id}_log.ApplyBQSR  2>&1 
		
		echo "end BQSR for ${id}" `date`
	fi
done



















































