# plasma
**1. Read Processing**
* _Trimming Adapters and low quality reads_

`trimmomatic PE -phred33 <SAMPLE_RAW_R1.fastq> <SAMPLE_RAW_R2.fastq> <SAMPLE_PAIRED_R1.fastq> <SAMPLE_UNPAIRED_R1.fastq>  <SAMPLE_PAIRED_R2.fastq> <SAMPLE_UNPAIRED_R2.fastq> ILLUMINACLIP:<ADAPTER_REFERENCE.txt>:2:40:15: SLIDINGWINDOW:4:30 LEADING:30 TRAILING:30 MINLEN:35`

* _Trimming to specific length_

`fastx_trimmer -l 145 -i <SAMPLE_PAIRED_R1.fastq> -o <SAMPLE_TRIMMED_R1.fastq>`

`fastx_trimmer -l 145 -i <SAMPLE_PAIRED_R2.fastq> -o <SAMPLE_TRIMMED_R2.fastq>`

* _Alignment_

`bwa mem -t 8 <REFERENCE_hg19.fasta> <SAMPLE_TRIMMED_R1.fastq> <SAMPLE_TRIMMED_R2.fastq> | samtools view -b -S -F 0x800 | samtools sort <SAMPLE_FILTERED_SORTED.bam> `

* _Mark Duplicates_

`picard MarkDuplicates INPUT=<SAMPLE_FILTERED_SORTED.bam> OUTPUT=<SAMPLE.bam> ASSUME_SORTED=true`

**2. Somatic variant analysis with VarScan2**

* _Generating mpileup_

`samtools mpileup -q 2 -f <REFERENCE_hg19.fasta> <SAMPLE.bam> > <SAMPLE.pileup>`

* _VarScan2_

`java -jar VarScan.v2.4.4.jar somatic <GERMLINE_SAMPLE.pileup> <TUMOUR_SAMPLE.pileup> <SAMPLE.snp> --min-var-freq 0.01`

`java -jar VarScan.v2.4.4.jar processSomatic <SAMPLE.snp> --min-tumor-freq 0.01 --max-normal-freq 0.00`

`awk 'BEGIN{OFS="\t"} NR!=1 {print $1, $2, $2}' <SAMPLE.snp.Somatic.hc> > <SAMPLE.snp.Somatic.hc.list>`

`bam-readcount -q 2 -f <REFERENCE_hg19.fasta> -l <SAMPLE.snp.Somatic.hc.list> <SAMPLE.bam> > <SAMPLE.somatic_hc.readcount>`

* _Filtering somatic calls from Varscan2_

`java -jar VarScan.v2.4.4.jar fpfilter <SAMPLE.snp.Somatic.hc> <SAMPLE.somatic_hc.readcount> --output-file <SAMPLE_fpfiltered> --min-var-freq 0.01 --min-strandedness 0.0`

`awk '{if($10 >= 5 && $6 == 0 && $5 + $6 >= 10 && $9 + $10 >= 10) {print $0}}' <SAMPLE_fpfiltered> > <SAMPLE_fpfiltered_somatic.snp>`

* _Extraction of shared variants (i.e. variants shared between plasma and tumour)_

`awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' <SAMPLE_PLASMA_fpfiltered_somatic.snp> <SAMPLE_TUMOUR_fpfiltered_somatic.snp> > <SAMPLE_PLASMA_fpfiltered_somatic_shared.snp>`

* _Extraction of unique variants (i.e. variants not shared between plasma and tumour)_

`awk 'FNR==NR{a[$1,$2]++}FNR!=NR && !a[$1,$2]{print}' <SAMPLE_TUMOUR_fpfiltered_somatic.snp> <SAMPLE_PLASMA_fpfiltered_somatic.snp> > <SAMPLE_PLASMA_fpfiltered_somatic_unique.snp>`

**3. Extraction of somatic reads**

These reads were extracted using JAPSA (https://github.com/mdcao/japsa) and somatic reads extraction tool can be deployed using script name jsa.hts.aareads.

`jsa.hts.aareads --input <SAMPLE.bam> --reference <REFERENCE_hg19.fasta>  --vcf <SAMPLE_fpfiltered_somatic.snp> --output <SAMPLE_somatic_reads.sam>`

`samtools view -b -S <SAMPLE_somatic_reads.sam> | samtools sort -@ 8 -o <SAMPLE_somatic_reads_sorted.bam> && samtools index <SAMPLE_somatic_reads_sorted.bam> <SAMPLE_somatic_reads_sorted.bai> ` 

**4. Annotation of variants with Annovar**

`awk '{print $1,$2,$2,$3,$4 }' <SAMPLE_fpfiltered_somatic.snp> > <SAMPLE.annovarlist.txt>`

`perl annotate_variation.pl -out <SAMPLE.annotation> -build hg19 <SAMPLE.annovarlist.txt> <DIRECTORY_PATH_TO_humandb_hg19>`
