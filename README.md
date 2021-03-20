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
case3_germline_1.fastq.gz      SRR3182430_2.fastq.gz.partial
case3_germline_2.fastq.gz      SRR3182430.fastq.gz.partial
case4_germline_1.fastq.gz      SRR3182431_1.fastq.gz.partial
fq_download.sh		       SRR3182431_2.fastq.gz.partial
nohup.out		       SRR3182431.fastq.gz.partial
SRR3182418_1.fastq.gz.partial  SRR3182432_1.fastq.gz.partial
SRR3182418_2.fastq.gz.partial  SRR3182432_2.fastq.gz.partial
SRR3182418.fastq.gz.partial    SRR3182432.fastq.gz.partial
SRR3182419_1.fastq.gz.partial  SRR3182433_1.fastq.gz.partial
SRR3182419_2.fastq.gz.partial  SRR3182433_2.fastq.gz.partial
SRR3182419.fastq.gz.partial    SRR3182433.fastq.gz.partial
SRR3182420_2.fastq.gz.partial  SRR3182434_1.fastq.gz.partial
SRR3182422_1.fastq.gz.partial  SRR3182434_2.fastq.gz.partial
SRR3182422_2.fastq.gz.partial  SRR3182434.fastq.gz.partial
SRR3182422.fastq.gz.partial    SRR3182435_1.fastq.gz.partial
SRR3182423_1.fastq.gz.partial  SRR3182435_2.fastq.gz.partial
SRR3182423_2.fastq.gz.partial  SRR3182435.fastq.gz.partial
SRR3182423.fastq.gz.partial    SRR3182436_1.fastq.gz.partial
SRR3182424_1.fastq.gz.partial  SRR3182436_2.fastq.gz.partial
SRR3182424_2.fastq.gz.partial  SRR3182436.fastq.gz.partial
SRR3182424.fastq.gz.partial    SRR3182437_1.fastq.gz.partial
SRR3182425_1.fastq.gz.partial  SRR3182437_2.fastq.gz.partial
SRR3182425_2.fastq.gz.partial  SRR3182438_1.fastq.gz.partial
SRR3182425.fastq.gz.partial    SRR3182438_2.fastq.gz.partial
SRR3182426_1.fastq.gz.partial  SRR3182438.fastq.gz.partial
SRR3182426_2.fastq.gz.partial  SRR3182439_1.fastq.gz.partial
SRR3182426.fastq.gz.partial    SRR3182439_2.fastq.gz.partial
SRR3182427_1.fastq.gz.partial  SRR3182439.fastq.gz.partial
SRR3182427_2.fastq.gz.partial  SRR3182440_1.fastq.gz.partial
SRR3182427.fastq.gz.partial    SRR3182440_2.fastq.gz.partial
SRR3182428_1.fastq.gz.partial  SRR3182440.fastq.gz.partial
SRR3182428_2.fastq.gz.partial  SRR3182441_1.fastq.gz.partial
SRR3182428.fastq.gz.partial    SRR3182441_2.fastq.gz.partial
SRR3182429_1.fastq.gz.partial  SRR3182441.fastq.gz.partial
SRR3182429_2.fastq.gz.partial  SRR3182443_1.fastq.gz.partial
SRR3182429.fastq.gz.partial    SRR3182443_2.fastq.gz.partial
SRR3182430_1.fastq.gz.partial  SRR3182443.fastq.gz.partial
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





















