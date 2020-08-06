echo Virus mapping, variant calling, and analysis pipeline by Daniel Ward \(daniel.ward1@lshtm.ac.uk\)

for f in $(ls *_R1_*.fastq | sed 's/_R1_001.fastq//' | sort -u); do  

echo ------------NOW TRIMMING AND MAPPING------------	

java -jar ../trimmomattic/trimmomatic-0.38.jar PE -phred33 -trimlog ../trimmomattic_logs/${f}_001.fastq.log.txt ${f}_R1_001.fastq ${f}_R2_001.fastq  ../trimmomattic_output/${f}_R1_001.fastq ../trimmomattic_output/unpaired/${f}_R1_001.unpaired.fastq ../trimmomattic_output/${f}_R2_001.fastq ../trimmomattic_output/unpaired/${f}_R2_001.unpaired.fastq ILLUMINACLIP:../trimmomattic/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 HEADCROP:30 &> /dev/null

bwa index ../reference_sequence/ref.fasta &> /dev/null

bwa mem ../reference_sequence/ref.fasta ../trimmomattic_output/${f}_R1_001.fastq ../trimmomattic_output/${f}_R2_001.fastq  | samtools view -q 15 -b -S  -  | samtools sort - -o ../samtools_output/${f}.sorted.bam

samtools index ../samtools_output/${f}.sorted.bam 

samtools faidx ../reference_sequence/ref.fasta

samtools stats -c 1,100,10  ../samtools_output/${f}.sorted.bam | grep 'COV\|SN' > ../samtools_stats/${f}.stats

echo -------

echo --------PERFORMING LOFREQ ANALYSIS----------

#../lofreq/lofreq_star-2.1.2/bin/lofreq call -f ../reference_sequence/ref.fasta -o ../lofreq_output/${f}.lofreq.vcf  ../samtools_output/${f}.sorted.bam

#bgzip ../lofreq_output/${f}.lofreq.vcf

#bcftools index ../lofreq_output/${f}.lofreq.vcf.gz

#bcftools stats ../lofreq_output/${f}.lofreq.vcf.gz > ../bcftools_stats/${f}.bcfstats

#bedtools genomecov -bga -ibam ../samtools_output/${f}.sorted.bam -g ../reference_sequence/ref.fasta.fai | awk '$4 < 10' > ../bedtools_stats/${f}.zero.bed

#bcftools consensus -f ../reference_sequence/ref.fasta -m ../bedtools_stats/${f}.zero.bed  ../lofreq_output/${f}.lofreq.vcf.gz -o ../consensus_sequence_output/${f}.fasta 

echo SAMPLE ${f} COMPLETE; done

echo ----SAMPLE PROCESSING COMPLETE----
echo now run 'data_analysis.sh'
