fasta=cell-line-sequence.fasta
name=cell-line-db

makeblastdb -in $fasta -input_type fasta -dbtype nucl -out $name
