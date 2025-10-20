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
## 3. Run BLASTP Alignments
```bash
# Create BLAST Databases
ml blast/2.16.0

makeblastdb -in wheat-query.protein.fa -dbtype prot -out wheat_db

# Run BLASTP
blastp -query wheat-subject-protein-sequences.fasta -db wheat-query_db -out wheat-subject-query.blast -evalue 1e-5 -outfmt 6 -num_threads 4

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

## 4. Prepare GFF Files
```bash
awk -F'\t' '$3=="gene"{
  split($9,a,";");
  id="";
  for(i in a)
    if(a[i]~/^ID=/){
      id=a[i];
      sub(/^ID=/,"",id)
    }
  if(id!="")
    print $1, id, $4, $5
}' OFS='\t' wheat-subject.high.gff3 > wheat-subject-chr.gff # Do the same for wheat-query

# Change chr to Ta
awk 'BEGIN{OFS="\t"}{$1 = gensub(/^chr/, "Ta", 1, $1); print}' wheat-subject-chr.gff > wheat-subject-Ta.gff

# Change Chr to Cs
awk 'BEGIN{OFS="\t"}{$1 = gensub(/^Chr/, "Cs", 1, $1); print}' wheat-query-Chr.gff > wheat-query-Cs.gff

# Filter out all lines where the 1st column (chromosome name) doesn’t start with "chr" and keep only the true chromosomes (chr1A–chr7D)
awk 'BEGIN{FS=OFS="\t"} $1 ~ /^Ta[1-7][ABD]?$/ {print}' wheat-subject-Ta.gff > wheat-subject.gff

awk 'BEGIN{FS=OFS="\t"} $1 ~ /^Cs[1-7][ABD]?$/ {print}' wheat-query-Cs.gff > wheat-query.gff
```

## 5. Normalize BLAST IDs to match your GFF
```bash
awk 'BEGIN{OFS="\t"}{
  # Fix query (col 1)
  gsub(/\.[0-9]+_/, "_", $1)        # ... .1_  ->  _
  gsub(/(LC|HC)\.[0-9]+_/, "_", $1) # ... LC.1_ or HC.1_ -> _
  gsub(/(LC|HC)_/, "_", $1)         # fallback: LC_ or HC_ -> _

  # Fix subject (col 2)
  gsub(/\.[0-9]+_/, "_", $2)
  gsub(/(LC|HC)\.[0-9]+_/, "_", $2)
  gsub(/(LC|HC)_/, "_", $2)

  print
}' glenn-csv21.blast > glenn-csv21.norm.blast
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
