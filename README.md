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

`java -jar VarScan.v2.3.9.jar somatic --min-var-freq 0.01 <GERMLINE_SAMPLE.pileup>  <TUMOUR_SAMPLE.pileup>  <SAMPLE.vcf.snp>`

* _Filtering somatic calls from Varscan2_

`awk '{ if(0 >= 5 &&  == 0 &&  +  >= 10 &&  + 0 >= 10 && 8 >= 2 && 9 >= 2) { print  }}' <SAMPLE.vcf.snp> | grep -i somatic > <SAMPLE.filtered.vcf.snp>`

* _Extraction of shared variants (i.e. variants shared between plasma and tumour)_

`awk 'FNR==NR{a[,]=;next}{if(b=a[,]){print b}}' <SAMPLE_PLASMA.filtered.vcf.snp> <SAMPLE_TUMOUR.filtered.vcf.snp>  > <SAMPLE_PLASMA.shared.vcf.snp>`

* _Extraction of unique variants (i.e. variants not shared between plasma and tumour)_

`awk 'FNR==NR{a[,]++}FNR!=NR && !a[,]{print}' <SAMPLE_PLASMA.filtered.vcf.snp> <SAMPLE_TUMOUR.filtered.vcf.snp>  > <SAMPLE_PLASMA.unique.vcf.snp>`

**3. Extraction of somatic reads**

These reads were extracted using JAPSA (https://github.com/mdcao/japsa) and somatic reads extraction tool can be deployed using script name jsa.hts.aareads.

`jsa.hts.aareads --input <SAMPLE.bam> --reference <REFERENCE_hg19.fasta>  --vcf <SAMPLE.filtered.vcf.snp> --output <SAMPLE.shared.somatic_reads.sam>`

`samtools view -b -S <SAMPLE.shared.somatic_reads.sam> | samtools sort <SAMPLE.shared.somatic_reads.sam.bam>` 

**4. Annotation of variants with Annovar**

`awk '{print ,,,,}' <SAMPLE.filtered.vcf.snp> > <SAMPLE.annovarlist.txt>`

`annotate_variation.pl -out <SAMPLE.annotation> -build hg19 <SAMPLE.annovarlist.txt> <DIRECTORY_PATH_TO_humandb_hg19>`
