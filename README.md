### codonUsageGrapher.py
deprecated mess

#### how to use:
not viable

### getcodons.py
visualise the codon usage from  DNA sequences recorded in fasta format
with sunburst diagrams

#### how to use:
create a list of ```string``` format filenames to parse and
throw it into the hardcoded ```list``` variable ```files```
at the top under the ```#HARDCODED``` section. Make sure the files are
within the same directory as the script or make sure the script knows
how to find them.

#### features:

- find startcodon
- create sliding window loop iterating withing frame from startcodon index
- keep track of:
	- how many times a codon occurs
	- the amount of aminoacids when all synonymous codons are added up
	- amount of errorrs/unknown codons
                  
#### hard-coded values 
Since this is just a simple script, some values are hard-coded. some can be found at the top of the script under #hardcoded. The following features are (still) hard-coded:

```files```
