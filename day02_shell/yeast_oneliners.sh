   #!/usr/bin/env bash
   # Day 2: one-liners on the yeast genome (Ensembl release 63, assembly R64-1-1)
   # Run from the repo root: bash day02_shell/yeast_oneliners.sh

   FA=data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
   GTF=data/Saccharomyces_cerevisiae.R64-1-1.63.gtf

   echo "1) number of sequences in the FASTA"
   grep -c "^>" "$FA"

   echo "2) length of each sequence"
   awk '/^>/ { if (name != "") print name, total; name = $1; total = 0 }
        !/^>/ { total += length($0) }
        END { print name, total }' "$FA"

   echo "3) genes per chromosome"
   awk '$3=="gene"' "$GTF" | cut -f1 | sort | uniq -c

   echo "4) total number of genes"
   awk '$3=="gene"' "$GTF" | wc -l