#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=128GB
#SBATCH -c 8

# Module
module load gatk
module load samtools

FILE_DIR=$(echo $1 | awk 'BEGIN{FS=OFS="/"}{NF--; print}')
FILE_NAME=$(echo $1 |  awk -F/ '{print $NF}')
PREFIX=${FILE_NAME%.bam}

# echo ${FILE_DIR}
# echo ${FILE_NAME}

cd ${FILE_DIR}
# echo $(pwd)
# echo ${PREFIX}

# PathSeq filtering
ulimit -c unlimited

gatk --java-options "-Xmx80g"  PathSeqFilterSpark \
        --input ${FILE_NAME} \
        --paired-output ${PREFIX}_PathSeq_filtered_paired.bam \
        --unpaired-output ${PREFIX}_PathSeq_filtered_unpaired.bam \
        --min-clipped-read-length 60 \
        --kmer-file /mnt/rstor/genetics/LaFramboiseLab/dxn150/PathSeq/ref_file/GRCh37-lite.hss \
        --filter-bwa-image /mnt/rstor/genetics/LaFramboiseLab/dxn150/PathSeq/ref_file/GRCh37-lite.fa.img

echo "Pass PathSeq"

# Fastq conversion
samtools sort -n ${PREFIX}_PathSeq_filtered_paired.bam | samtools fastq - -1 ${PREFIX}_PathSeq_filtered_paired_R1.fastq -2 ${PREFIX}_PathSeq_filtered_paired_R2.fastq
samtools sort -n ${PREFIX}_PathSeq_filtered_unpaired.bam | samtools fastq - > ${PREFIX}_PathSeq_filtered_unpaired.fastq

# Bowtie2 Unmapped CHM13
bowtie2 -x /home/dxn150/ncbi_dataset/data/GCF_009914755.1/CHM13_genome \
  -1 ${PREFIX}_PathSeq_filtered_paired_R1.fastq \
  -2 ${PREFIX}_PathSeq_filtered_paired_R2.fastq \
  --un-conc-gz ${PREFIX}_CHM13_unmapped_R%.fastq.gz \
  -p 8 > /dev/null

bowtie2 -x /home/dxn150/ncbi_dataset/data/GCF_009914755.1/CHM13_genome \
  -U ${PREFIX}_PathSeq_filtered_unpaired.fastq \
  --un ${PREFIX}_CHM13_unmapped_unpaired.fastq \
  -p 8 > /dev/null

# KrakenUniq filter CHM13
krakenuniq --db /scratch/pioneer/users/dxn150/krakenuniq_CHM13 --threads 8 \
    --output ${PREFIX}_pe_chm13.kraken \
    --report ${PREFIX}_pe_chm13.report \
    --paired ${PREFIX}_PathSeq_filtered_paired_R1.fastq ${PREFIX}_PathSeq_filtered_paired_R2.fastq \
    --unclassified-out ${PREFIX}_kraken_unclassified_paired_R#.fastq	

krakenuniq --db /scratch/pioneer/users/dxn150/krakenuniq_CHM13 --threads 8 \
    --output ${PREFIX}_unpaired_chm13.kraken \
    --report ${PREFIX}_unpaired_chm13.report \
    ${PREFIX}_CHM13_unmapped_unpaired.fastq \
    --unclassified-out ${PREFIX}_kraken_unclassified_unpaired.fastq

echo "Pass Kraken filter human"

# KrakenUniq standard database
krakenuniq --db /scratch/pioneer/users/dxn150/krakenuniq_standard --threads 8 --paired \
    --output ${PREFIX}_kraken_paired_standard_db.kraken \
    --report ${PREFIX}_kraken_paired_standard_db.report \
    ${PREFIX}_kraken_unclassified_paired_R1.fastq ${PREFIX}_kraken_unclassified_paired_R2.fastq

krakenuniq --db /scratch/pioneer/users/dxn150/krakenuniq_standard --threads 8 \
    --output ${PREFIX}_kraken_unpaired_standard_db.kraken \
    --report ${PREFIX}_kraken_unpaired_standard_db.report \
    ${PREFIX}_kraken_unclassified_unpaired.fastq

# KrakenUniq remove contaminants
python --version
python /mnt/rstor/genetics/LaFramboiseLab/dxn150/Scripts/KrakenTools-1.2/extract_kraken_reads.py --taxid 0 9606 519051 28198 28200 82135 307456 646534 1816686 1386088 1898961 1580 1613 1589 1590 1736525 1245 1777 146021 197461 721133 60136 1622075 1662865 86188 1167641 1510458 283734 220697 46123 33951 222 522 12916 469 153265 2040 1033 59753 114627 1347386 92793 1663 35823 76890 182269 12960 146937 191 55087 1386 958 532 84756 1678 34098 169215 85413 374 55080 1696 41275 2755 32008 2747 75 79328 59732 283 334107 1716 77583 106589 2034 281915 42 978 56112 1298 80865 36853 79206 46913 37914 75654 120831 212791 106591 547 561 66831 237 204456 1860 940550 501022 1742989 397456 475794 963 303379 274591 1649459 1004300 135575 53457 188905 29580 1470540 32257 57493 53452 424207 47251 1907117 88 131079 28073 400634 149698 68287 407 407 184923 133 16 378210 33882 48073 1269 44471 29404 58050 423349 64001 354354 914 1543704 165696 528 376469 1158 846 44249 227873 1822464 265 249411 1621534 1004302 361607 122277 84567 47494 335058 302485 657 28100 52972 289201 1743 53246 286 83618 497 67572 48736 379 75309 1661425 1827 1068 212743 29574 323620 28067 28065 65047 644355 318147 2843398 1436289 36862 1912216 391952 1279 29407 1073 34008 125216 215579 1203 620 207599 165695 13687 165697 40323 1301 1054211 114248 2425 99479 85079 111782 158851 2060 401469 34072 364316 343873 338 626 39492 1443441 1177154 95605 1545443 518 29459 199 2751 394935 546 237258 1502 1268241 1871336 413503 28141 636 31973 1736502 859 1314686 397457 727 1736490 1629124 137716 86102 171674 548 571 1463165 1596 1363 1358 1736496 1639 242606 69966 1126 1549810 34062 28045 1834067 1834082 1370121 495 489 28449 204799 37329 494023 549 553 92647 1629550 1786003 554 51663 587 54291 1736528 416169 95607 2047 28901 615 1410620 95606 1280 29388 1282 246432 1283 1292 538381 307486 11757 11780 11885 12177 12235 12267 99182 342409 369960 1891754 11788 11807 11786 11801 44561 368797 61673 168238 35269 --exclude --include-children -k ${PREFIX}_kraken_paired_standard_db.kraken -r ${PREFIX}_kraken_paired_standard_db.report -s1 ${PREFIX}_kraken_unclassified_paired_R_1.fastq -s2 ${PREFIX}_kraken_unclassified_paired_R_2.fastq --fastq-output -o ${PREFIX}_paired_non_human_non_contaminant_R1.fastq -o2 ${PREFIX}_paired_non_human_non_contaminant_R2.fastq

python /mnt/rstor/genetics/LaFramboiseLab/dxn150/Scripts/KrakenTools-1.2/extract_kraken_reads.py --taxid 0 9606 519051 28198 28200 82135 307456 646534 1816686 1386088 1898961 1580 1613 1589 1590 1736525 1245 1777 146021 197461 721133 60136 1622075 1662865 86188 1167641 1510458 283734 220697 46123 33951 222 522 12916 469 153265 2040 1033 59753 114627 1347386 92793 1663 35823 76890 182269 12960 146937 191 55087 1386 958 532 84756 1678 34098 169215 85413 374 55080 1696 41275 2755 32008 2747 75 79328 59732 283 334107 1716 77583 106589 2034 281915 42 978 56112 1298 80865 36853 79206 46913 37914 75654 120831 212791 106591 547 561 66831 237 204456 1860 940550 501022 1742989 397456 475794 963 303379 274591 1649459 1004300 135575 53457 188905 29580 1470540 32257 57493 53452 424207 47251 1907117 88 131079 28073 400634 149698 68287 407 407 184923 133 16 378210 33882 48073 1269 44471 29404 58050 423349 64001 354354 914 1543704 165696 528 376469 1158 846 44249 227873 1822464 265 249411 1621534 1004302 361607 122277 84567 47494 335058 302485 657 28100 52972 289201 1743 53246 286 83618 497 67572 48736 379 75309 1661425 1827 1068 212743 29574 323620 28067 28065 65047 644355 318147 2843398 1436289 36862 1912216 391952 1279 29407 1073 34008 125216 215579 1203 620 207599 165695 13687 165697 40323 1301 1054211 114248 2425 99479 85079 111782 158851 2060 401469 34072 364316 343873 338 626 39492 1443441 1177154 95605 1545443 518 29459 199 2751 394935 546 237258 1502 1268241 1871336 413503 28141 636 31973 1736502 859 1314686 397457 727 1736490 1629124 137716 86102 171674 548 571 1463165 1596 1363 1358 1736496 1639 242606 69966 1126 1549810 34062 28045 1834067 1834082 1370121 495 489 28449 204799 37329 494023 549 553 92647 1629550 1786003 554 51663 587 54291 1736528 416169 95607 2047 28901 615 1410620 95606 1280 29388 1282 246432 1283 1292 538381 307486 11757 11780 11885 12177 12235 12267 99182 342409 369960 1891754 11788 11807 11786 11801 44561 368797 61673 168238 35269 --exclude --include-children -k ${PREFIX}_kraken_unpaired_standard_db.kraken -r ${PREFIX}_kraken_unpaired_standard_db.report -s ${PREFIX}_kraken_unclassified_unpaired.fastq --fastq-output -o ${PREFIX}_unpaired_non_human_non_contaminant.fastq

echo "Pass Kraken remove contaminant"

# KrakenUniq classified
krakenuniq --db /scratch/pioneer/users/dxn150/krakenuniq_standard --threads 8 --paired \
    --output ${PREFIX}_kraken_paired_non_human_non_contaminant_standard_db.kraken \
    --report ${PREFIX}_kraken_paired_non_human_non_contaminant_standard_db.report \
    ${PREFIX}_paired_non_human_non_contaminant_R1.fastq ${PREFIX}_paired_non_human_non_contaminant_R2.fastq

krakenuniq --db /scratch/pioneer/users/dxn150/krakenuniq_standard --threads 8 \
    --output ${PREFIX}_kraken_unpaired_non_human_non_contaminant_standard_db.kraken \
    --report ${PREFIX}_kraken_unpaired_non_human_non_contaminant_standard_db.report \
    ${PREFIX}_unpaired_non_human_non_contaminant.fastq 

echo "Pass Kraken classification"

# Combine KrakenUniq report
python /mnt/rstor/genetics/LaFramboiseLab/dxn150/Scripts/KrakenTools-1.2/combine_kreports.py --no-headers --only-combined -r ${PREFIX}_paired_non_human_non_contaminant_standard_db.report ${PREFIX}_unpaired_non_human_non_contaminant_standard_db.report -o ${PREFIX}_combined_non_human_non_contaminant_standard_db.report

