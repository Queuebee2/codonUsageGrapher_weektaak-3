### codonUsageGrapher
A simple script to create nested pie charts displaying the codon usage ratio for each amino acid, alongside the amino acid ratio within the coding region of a gene.


### how to use:
you don't use it, it's not done.


                  
#### hard-coded values 
Since this is just a simple script, some values are hard-coded. some can be found at the top of the script under #hardcoded. The following features are (still) hard-coded:

 - the folder name that holds sequences
 - "```.fasta```" extentions on filenames get stripped
 - a filter keyword ( "```CDS```" by default)
 - a delimiter ( "```_```" by default)
 - the order in which files are named
 	- ```filename[0]``` is assumed as naming the part
 	- ```filename[1]``` is assumed to be the descriptive part
 	
