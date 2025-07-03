annotation=reference/*gff
reference=reference/*fasta

#degeneracy annotation for coding transcripts
degenotate/degenotate.py -a $annotation -g $reference -o out.degenotate/ -d --overwrite

#extract 4-fold degenerate sites
awk '{if($5==4) print $1":"$2}' out.degenotate/degeneracy-all-sites.bed > 4-fold.sites
