The .c and .h files in this directory are needed to construct the programs
enhancer and makechrmutr.

makechrmutr is used to convert the .gbk files for the chromosomes of an
organism into a compressed form which is used by the program enhancer. You
can construct the program with the makefile make.makechrmutr. 

enhancer is used to search a genome for clusters of motifs. The motifs are
specified using IUPAC codes for the nucleotides. After the program searches
for all occurrances of the motifs, the user is invited to describe the
parameters of a cluster: how many of the specified motifs must occur in the 
cluster, and the maximum size of the cluster in terms of nucleotides.

After enhancer terminates, the file "hashresults" contains details of the
execution of enhancer. Save this file elsewhere if you want to save the results
for later analysis. Every execution of enhancer overwrites "hashresults".

enhancer may be built by using the makefile make.enhancer. The executable
enhancer.exe will run on a MAC.
