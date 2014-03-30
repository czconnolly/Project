#Script to set the values for AMPS alignments
output_file=SLC45A3alignment.out
print_horizontal=
mode=pairwise
matrix_file=md.mat
pairwise_random=100,100,1
gap_penalty=8.0
constant =8
seq_file=SLC45A3.pir
pairwise_align_file=SLC45A3_alignment.fasta
