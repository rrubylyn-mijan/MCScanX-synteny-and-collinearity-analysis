# MCScanX-synteny-and-collinearity-analysis
This workflow identifies syntenic relationships and collinear gene blocks between wheat genomes using MCScanX.
The analysis combines GFF3 annotations, protein sequence alignments (BLASTP), and collinearity visualization via Circos or SynVisio.

## Install MCScanX
```bash
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanX
make

# Verify Installation
./MCScanX

# Access Downstream Tools
cd downstream_analyses
ls
java -version
java dot_plotter
```

## 2. Prepare Protein Sequences
```bash
# Extract Coding Sequences (CDS)
conda create -n myenv -c bioconda gffread
conda activate myenv
gffread --version

# Translate CDS → Proteins
conda create -n emboss_env -c bioconda emboss
conda activate emboss_env
transeq -sequence wheat.subject.cds.fa -outseq wheat-subject-protein-sequences.fasta
transeq -sequence wheat.query.cds.fa -outseq wheat-query-protein-sequences.fasta
```

## 3. Combine Protein Files
```bash
cat wheat-subject-protein-sequences.fasta wheat-query.protein.fa > combined-wheats-subject-query.fasta
```

## 4. Prepare GFF Files
```bash
awk '$3 == "gene" {split($9, a, ";"); print $1, a[1]"_1", $4, $5}' OFS='\t' wheat-subject.high.gff3 > fixed-wheat-subject.gff

# Fix Gene Naming for Compatibility
awk '{split($2, a, "="); print $1"\t"$3"\t"$4"\t"a[2]}' wheat-query.gff > wheat-query-fixed.gff

# Combine Annotations
cat wheat-subject.gff wheat-query-fixed.gff > combined-wheats-subject-query.gff
```

## 5. Run BLASTP Alignments
```bash
# Create BLAST Databases
ml blast/2.16.0

makeblastdb -in combined-wheats-subject-query.fasta -dbtype prot -out combined_db
makeblastdb -in wheat-query.protein.fa -dbtype prot -out wheat_db

# Run Self and Cross BLASTP
blastp -query wheat-subject-protein-sequences.fasta -db wheat-subject -out wheat-subject.blast -evalue 1e-5 -outfmt 6 -num_threads 4
blastp -query wheat-query-protein-sequences.fasta -db wheat-query_db -out wheat-subject-query-blast-results.txt -evalue 1e-5 -outfmt 6 -num_threads 4

# Example PBS Script
#!/bin/bash
#PBS -q default
#PBS -N blastp_sumai3_chinese
#PBS -l select=1:mem=100gb:ncpus=1
#PBS -l walltime=24:00:00
#PBS -M firstname.lastname@ndsu.edu
#PBS -m abe
#PBS -W group_list=x-ccast-prj-name

ml blast/2.16.0
cd /directory/this/saved/mcscanx_results_wheat-subject-wheat-query
blastp -query wheat-subject-protein-sequences.fasta -db wheat-query_db -out wheat-subject-query-blast-results.txt -outfmt 6 -evalue 1e-5 -num_threads 4
```

## 6. Run MCScanX
```bash
# Example PBS Job
#!/bin/bash
#PBS -q bigmem
#PBS -N mcscanx_sumai3
#PBS -l select=1:mem=300gb:ncpus=8
#PBS -l walltime=128:00:00
#PBS -M firstname.lastname@ndsu.edu
#PBS -m abe
#PBS -W group_list=x-ccast-prj-name

export PATH=$PATH:/directory/this/saved/MCScanX
cd /directory/this/saved/MCScanX

./MCScanX /directory/this/saved/mcscanx_result_combined_wheat-subject_wheat-query/combined -k 80 -g -2 -s 2 -e 1e-10 -m 0
```

## 7. Generate Circos Links
```bash
# Extract Collinear Gene Pairs
awk '
NR == FNR {
    gsub(/_1$/, "", $2);
    pos[$2] = $1 "\t" $3 "\t" $4;
    next;
}
(/^[[:space:]]*[0-9]+-[[:space:]]*[0-9]+:/) {
    gene1 = $3; gene2 = $4;
    gsub(/_1$/, "", gene1); gsub(/_1$/, "", gene2);
    if (gene1 in pos && gene2 in pos) {
        print "link" ++n, pos[gene1];
        print "link" n, pos[gene2];
    }
}' wheat.gff wheat.collinearity > circos_links.txt

# Rename Chromosomes for Circos
awk '{gsub(/^chr/, "ta", $2); print}' circos_links.txt > collinearity_links.txt

# Assign Colors per Chromosome
awk '
BEGIN {
    color["1"]="#C379FF"; color["2"]="#EDDE74"; color["3"]="#E97E7E";
    color["4"]="#E89DAF"; color["5"]="#F2DACB"; color["6"]="#93B07"; color["7"]="#77DD77";
}
{
    if (match($2, /ta([1-7])[ABD]/, m))
        print $1, $2, $3, $4, "color="color[m[1]];
    else print $0;
}' collinearity_links.txt > colored_output.txt
```

## 8. Run MCScanX on Individual Genomes
```bash
ml blast-plus/2.14.1
blastn -query wheat-subject.fasta -db wheat-query_db -evalue 1e-10 -outfmt 6 -out wheat-subject.blast -num_threads 8

(Same format for other genomes; replace genome file names accordingly).
```

## 9. Prepare SynVisio Input
```bash
gene_id,chr,start,end,block_id
Gene1,Chr01,100,200,Block1
Gene2,Chr02,300,400,Block1

(Then upload to SynVisio to explore syntenic blocks interactively).
```

Maintainer:

Ruby Mijan
