
The lists of the genes that are upregulated due to protoplasting are:

1) Ox_leaf_protoplast_v12_final_table_DEGs_2reps_final_August_2022.csv : This list will be used for Cardamine hirsuta single cell analyses. I have used the following RNASeq files:
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/4826/4826.A
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/4826/4826.B
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/4867/4867.A
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/4867/4867.B

The scripts used for the mapping and raw-read count are the following:
hisat2 -q -x /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/ChirsutaOX.fa -1 4826_A_run672_GGCTACAG_S81_L006_R1_001.fastq.gz,4826_A_run676_GGCTACAG_S19_L002_R1_001.fastq.gz -2 4826_A_run672_GGCTACAG_S81_L006_R2_001.fastq.gz,4826_A_run676_GGCTACAG_S19_L002_R2_001.fastq.gz -S 4826_A.sam
hisat2 -q -x /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/ChirsutaOX.fa -1 4826_B_run672_TCTGCTGT_S82_L006_R1_001.fastq.gz,4826_B_run676_TCTGCTGT_S20_L002_R1_001.fastq.gz -2 4826_B_run672_TCTGCTGT_S82_L006_R2_001.fastq.gz,4826_B_run676_TCTGCTGT_S20_L002_R2_001.fastq.gz -S 4826_B.sam
hisat2 -q -x /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/ChirsutaOX.fa -1 4867_A_run676_TCATTGAG_S32_L002_R1_001.fastq.gz,4867_A_run677_TCATTGAG_S160_L006_R1_001.fastq.gz -2 4867_A_run676_TCATTGAG_S32_L002_R2_001.fastq.gz,4867_A_run677_TCATTGAG_S160_L006_R2_001.fastq.gz -S 4867_A.sam
hisat2 -q -x /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/ChirsutaOX.fa -1 4867_B_run676_TACCGAGC_S33_L002_R1_001.fastq.gz,4867_B_run677_TACCGAGC_S161_L006_R1_001.fastq.gz -2 4867_B_run676_TACCGAGC_S33_L002_R2_001.fastq.gz,4867_B_run677_TACCGAGC_S161_L006_R2_001.fastq.gz -S 4867_B.sam
samtools view -b /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/4826_A.sam -o /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_leaf.bam
samtools view -b /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/4826_B.sam -o /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_protoplast.bam
samtools view -b 4867_A.sam -o 4867_A.bam
samtools view -b 4867_B.sam -o 4867_B.bam
samtools sort /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_leaf.bam -o /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_leaf.sorted.bam
samtools sort /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_protoplast.bam -o /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_protoplast.sorted.bam
samtools sort 4867_A.bam -o 4867_A.sorted.bam
samtools sort 4867_B.bam -o 4867_B.sorted.bam
samtools index /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_leaf.sorted.bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_leaf.sorted.bam.bai
samtools index /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_protoplast.sorted.bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/Ox_WT_protoplast.sorted.bam.bai
samtools index 4867_A.sorted.bam 4867_A.sorted.bam.bai
samtools index 4867_B.sorted.bam 4867_B.sorted.bam.bai
htseq-count -s reverse -r name --type=gene --idattr=ID -f bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/with_NEW_genome/4867_A.sorted.bam /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/hirsutaOX.v12.gff > 4867_A_raw_count_v12.txt
htseq-count -s reverse -r name --type=gene --idattr=ID -f bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/with_NEW_genome/4867_B.sorted.bam /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/hirsutaOX.v12.gff > 4867_B_raw_count_v12.txt
htseq-count -s reverse -r name --type=gene --idattr=ID -f bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/with_NEW_genome/Ox_WT_leaf.sorted.bam /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/hirsutaOX.v12.gff > /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/with_NEW_genome/Ox_WT_leaf_raw_count_v12.txt
htseq-count -s reverse -r name --type=gene --idattr=ID -f bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/with_NEW_genome/Ox_WT_protoplast.sorted.bam /biodata/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/CHRISTOS/GENOMES/CURRENT/hirsutaOX.v12.gff > /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4826/with_NEW_genome/Ox_WT_protoplast_raw_count_v12.txt

The DEGs were found by using EdgeR.

In total I have identified 1180 genes (>=2-fold; q < 0.05 and upregulated in protoplasts) to be excluded from the Cardamine hirsuta single cell analyses


2) Col0.leaf.vs.protoplast_table_2pseudorep_final_August_2022.csv : This list will be used for Arabidopsis thaliana single cell analyses. I have used the following RNASeq files:
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/4867/4867.C
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/4867/4867.D
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/4867/4867.E
/netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Project_4867/4867/4867.F

The scripts used for the mapping and raw-read count are the following:
hisat2 -q -x /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.fa -1 4867_C_run676_TGTGAAGA_S34_L002_R1_001.fastq.gz,4867_C_run677_TGTGAAGA_S181_L007_R1_001.fastq.gz -2 4867_C_run676_TGTGAAGA_S34_L002_R2_001.fastq.gz,4867_C_run677_TGTGAAGA_S181_L007_R2_001.fastq.gz -S 4867_C.sam
hisat2 -q -x /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.fa -1 4867_D_run676_GCAACATT_S35_L002_R1_001.fastq.gz,4867_D_run677_GCAACATT_S182_L007_R1_001.fastq.gz -2 4867_D_run676_GCAACATT_S35_L002_R2_001.fastq.gz,4867_D_run677_GCAACATT_S182_L007_R2_001.fastq.gz -S 4867_D.sam
hisat2 -q -x /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.fa -1 4867_E_run676_GCTCCTTG_S36_L002_R1_001.fastq.gz,4867_E_run677_GCTCCTTG_S183_L007_R1_001.fastq.gz -2 4867_E_run676_GCTCCTTG_S36_L002_R2_001.fastq.gz,4867_E_run677_GCTCCTTG_S183_L007_R2_001.fastq.gz -S 4867_E.sam
hisat2 -q -x /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.fa -1 4867_F_run676_GCACTGTC_S37_L002_R1_001.fastq.gz,4867_F_run677_GCACTGTC_S184_L007_R1_001.fastq.gz -2 4867_F_run676_GCACTGTC_S37_L002_R2_001.fastq.gz,4867_F_run677_GCACTGTC_S184_L007_R2_001.fastq.gz -S 4867_F.sam
samtools view -b 4867_C.sam -o 4867_C.bam
samtools view -b 4867_D.sam -o 4867_D.bam
samtools view -b 4867_E.sam -o 4867_E.bam
samtools view -b 4867_F.sam -o 4867_F.bam
samtools sort 4867_C.bam -o 4867_C.sorted.bam
samtools sort 4867_D.bam -o 4867_D.sorted.bam
samtools sort 4867_E.bam -o 4867_E.sorted.bam
samtools sort 4867_F.bam -o 4867_F.sorted.bam
samtools index 4867_C.sorted.bam 4867_C.sorted.bam.bai
samtools index 4867_D.sorted.bam 4867_D.sorted.bam.bai
samtools index 4867_E.sorted.bam 4867_E.sorted.bam.bai
samtools index 4867_F.sorted.bam 4867_F.sorted.bam.bai
htseq-count -s reverse -r name --type=gene --idattr=gene_id -f bam 4867_C.sorted.bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.gtf > 4867_C_raw_count.txt
htseq-count -s reverse -r name --type=gene --idattr=gene_id -f bam 4867_D.sorted.bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.gtf > 4867_D_raw_count.txt
htseq-count -s reverse -r name --type=gene --idattr=gene_id -f bam 4867_E.sorted.bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.gtf > 4867_E_raw_count.txt
htseq-count -s reverse -r name --type=gene --idattr=gene_id -f bam 4867_F.sorted.bam /netscratch/dep_tsiantis/grp_tsiantis/1_CURRENT_LAB_MEMBERS/Christos/RNASeq_RCO_RELATED/Arabidopsis_RCO.gtf > 4867_F_raw_count.txt

In total I have identified 2022 genes (>=2-fold; q < 0.05 and upregulated in protoplasts)to be excluded from the Arabidopsis single cell analyses


3) Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv : This list will be used for the integration analyses between Arabidopsis and Cardamine
This list was resulted by merging "Col0.leaf.vs.protoplast_table_2pseudorep_final_August_2022.csv" and "Ox_leaf_protoplast_v12_final_table_DEGs_2reps_final_August_2022.csv" and keeping the genes with orthologues in both species. In total I have identified 2270 genes (>=2-fold; q < 0.05 and upregulated in protoplasts)to be excluded from the integration analyses.


