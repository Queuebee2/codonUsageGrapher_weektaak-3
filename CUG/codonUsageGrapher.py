#Author:    Milain Lambers | Queuebee2
#Date:      05-03-2018
version = "0.5"


# rough idea
# - using 2 dictionaries: one to look up, one to update
# - find start for frame
# - start parsing, counting, keeping statistic values in dicts uptodate
# - find end and check if frame was % 3 == 0 length
# - create graph.



# changelog

#NEW
# and codon2Amino to translate
# added dictionaries 'emptyData' to keep track of self.codonUsage
#? 
# pasted dictionaries from  blok2, afvink3 (temporarily?) (U->T)
# added function:   createFastaobjects
# added function:   getFilenames
# started Fastafile class
# added function:   setCWD
# added function:   spacer
# added function:   main    (initial)
# added function:   findUp
# 5-3-2018
#OLD

#imports
import os
from pathlib import Path


#hardcoded globals
# delimiter = '_'
#
pathKeyword = 'sequences'
recursionDepth = 3
specialFilenames = ['arabidopsis-thaliana_chr2.fasta',
                    'escherichia-coli_genome.fasta',
                    'salmonella-enterica_genome.fasta',
                    'homo-sapiens_chr12.fasta']
                    

def main(verbosity=False):
    """
    """
    # find fasta files
    sequencesPath = findUp(pathKeyword)
    print("Found presumable path with sequences")
    print(sequencesPath)
    spacer()
    setCWD(sequencesPath)
    availableFastas =  getFilenames()

    # create objects to hold information
    objectsToGraph = createFastaObjects(availableFastas, 'CDS', specialFilenames)

    if verbosity:
        print("printing object attributes")
        for o in objectsToGraph:
            print('printing name:',o.name)
            o.printSummary()
            o.createDictionary()


    print("hopefully I did that right...")

# lots of errors because I realised this too late
"""
A = adenine
C = cytosine
G = guanine
T = thymine
R = G A (purine)
Y = T C (pyrimidine)
K = G T (keto)
M = A C (amino)
S = G C (strong bonds)
W = A T (weak bonds)
B = G T C (all but A)
D = G A T (all but C)
H = A C T (all but G)
V = G C A (all but T)
N = A G C T (any)
"""

emptyCodonUsage = {'TOTAL': 0,
                   'ERRORS':0,
 'Ala ': {'Ala _TOTAL': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0},
 'Arg': {'Arg_TOTAL': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0, 'AGG': 0},
 'Asn': {'Asn_TOTAL': 0, 'AAT': 0, 'AAC': 0},
 'Asp': {'Asp_TOTAL': 0, 'GAT': 0, 'GAC': 0},
 'Cys': {'Cys_TOTAL': 0, 'TGT': 0, 'TGC': 0},
 'Gln': {'Gln_TOTAL': 0, 'CAA': 0, 'CAG': 0},
 'GlT': {'GlT_TOTAL': 0, 'GAA': 0, 'GAG': 0},
 'Gly': {'Gly_TOTAL': 0, 'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0},
 'His': {'His_TOTAL': 0, 'CAT': 0, 'CAC': 0},
 'Ile': {'Ile_TOTAL': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0},
 'LeT': {'LeT_TOTAL': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0},
 'Lys': {'Lys_TOTAL': 0, 'AAA': 0, 'AAG': 0},
 'Met': {'Met_TOTAL': 0, 'ATG': 0},
 'Phe': {'Phe_TOTAL': 0, 'TTT': 0, 'TTC': 0},
 'Pro': {'Pro_TOTAL': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0},
 'Ser': {'Ser_TOTAL': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0},
 'Thr': {'Thr_TOTAL': 0, 'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0},
 'Trp': {'Trp_TOTAL': 0, 'TGG': 0},
 'Tyr': {'Tyr_TOTAL': 0, 'TAT': 0, 'TAC': 0},
 'Val': {'Val_TOTAL': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0},
 'Start': {'Start_TOTAL': 0, 'ATG': 0, 'CTG': 0, 'TTG': 0, 'GTG': 0, 'ATT': 0},
 'Stop': {'Stop_TOTAL': 0, 'TAG': 0, 'TGA': 0, 'TAA': 0}}
    
codon2Amino = {'GCT': 'Ala ', 'GCC': 'Ala ',
               'GCA': 'Ala ', 'GCG': 'Ala ',
               'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg',
               'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
               'AAT': 'Asn', 'AAC': 'Asn', 'GAT': 'Asp',
               'GAC': 'Asp', 'TGT': 'Cys', 'TGC': 'Cys',
               'CAA': 'Gln', 'CAG': 'Gln', 'GAA': 'GlT',
               'GAG': 'GlT', 'GGT': 'Gly', 'GGC': 'Gly',
               'GGA': 'Gly', 'GGG': 'Gly', 'CAT': 'His',
               'CAC': 'His', 'ATT': 'Ile', 'ATC': 'Ile',
               'ATA': 'Ile', 'TTA': 'LeT', 'TTG': 'LeT',
               'CTT': 'LeT', 'CTC': 'LeT', 'CTA': 'LeT',
               'CTG': 'LeT', 'AAA': 'Lys',
               'AAG': 'Lys', 'ATG': 'Met',
               'TTT': 'Phe', 'TTC': 'Phe',
               'CCT': 'Pro', 'CCC': 'Pro',
               'CCA': 'Pro', 'CCG': 'Pro',
               'TCT': 'Ser', 'TCC': 'Ser',
               'TCA': 'Ser', 'TCG': 'Ser',
               'AGT': 'Ser', 'AGC': 'Ser',
               'ACT': 'Thr', 'ACC': 'Thr',
               'ACA': 'Thr', 'ACG': 'Thr',
               'TGG': 'Trp', 'TAT': 'Tyr',
               'TAC': 'Tyr', 'GTT': 'Val',
               'GTC': 'Val', 'GTA': 'Val',
               'GTG': 'Val',
               'TAG': 'Stop', 'TGA': 'Stop', 'TAA': 'Stop'}


        
def createFastaObjects(filesToParse, keyword=None, special=None,delimiter='_' ,verbose=False,):
    """
    args
        keyword: a keyword that has to be in the name

        special: a list of strings containing exact names (for hardcoding purposes
        if parsing is too complicated)
    """
    createdStuff = [] # list of objects to return
    
    for fastaFilename in filesToParse:
        parts = fastaFilename.strip('.fasta').split(delimiter)
        name = parts[0]
        description = parts[1]
        if not keyword:
            print('creating Fastafile object with name:',name,'\tdesc:',description)
            fastaObject = Fastafile(name, description, fastaFilename)
            createdStuff.append(fastaObject)
        elif keyword in fastaFilename:
            print('creating Fastafile object with name:',name,'\tdesc:',description)
            fastaObject = Fastafile(name, description, fastaFilename)
            createdStuff.append(fastaObject)
        elif special:
            if fastaFilename in special:
                print('creating Fastafile object with name:',name,'\tdesc:',description)
                fastaObject = Fastafile(name, description, fastaFilename)
                createdStuff.append(fastaObject)
        else:
            if verbose: print("no keyword match with filename, no object.")
        
    spacer()
    return createdStuff


def getFilenames():
    """lists available files in cwd"""
    return [f for f in os.listdir('.') if os.path.isfile(f)] # list comprehension stolen from SO

def spacer():
    print("-------------------------")


def setCWD(path):
    os.chdir(path)
    print("working in \"", os.getcwd(), "\" now")
    spacer()
    
def findUp(searchterm, thispath=False, depth=recursionDepth ,verbose=False):
    """version 3, march 5, 2018
    changelog
    
        okt 8, 2017
        forgot to add 'return' before starting recursive findUP
        
        okt 9, 2017
        try to recursively search parent folders for the folder thats named as searchterm

        mar 5, 2018
        added: oschar for windows/linux compatibility
    """

    # get operating system type
    if os.name == 'nt':
        oschar = '\\'
    else:
        oschar = '/'


    # get path if none given
    if not thispath:
        thispath = os.getcwd()

    # end recursion
    if depth <= 0:
        raise RecursionError("recusion limit (set to "+str(recursionDepth)+") exceeded")
    print('looking for',searchterm,'in',thispath)
    for dirname, nested_dirnames, absolute_filenames in os.walk(thispath):
        if dirname.split(oschar)[-1] == searchterm:
            result = str(dirname)
            if verbose: print("\nfindUp found",result)
            # return filepath in string format
            return result
        
    if verbose: print('findUp going up a directory!')
    return findUp(searchterm, str(Path(thispath).parent), depth-1)


class Fastafile():
    """ docstring """
    def __init__(self, name, description, filename):
        self.name = name
        self.description = description
        self.filename = filename
        self.codonUsage = emptyCodonUsage.copy()

    def getName(self):
        return self.name

    def getDescription(self):
        return self.name

    def getFilename(self):
        return self.filename
        
    def printSummary(self):
        print('name:',self.name,'\tdesc:',self.description, '\thasfilename:', str(bool(self.filename)))


    def createDictionary(self, startCodons=False, stopCodons=False, verbose=False):
        """finds the start of a coding frame by using
        the sliding window algorithm until it finds the index
        of the first start codon and then creates a dictionary with
        total:count and aminoacid:{total:count,codonsynonym:count} key-value pairs

        changelog

        args

            sequenceFile: iOtextwrapper object of a fasta format
            DNA coding sequence.

            startcodons(default=False) bool/list: list of used
            codons to filter by for the start of the frame.
                False: uses all known startcodons

            stopcodons(default=False) bool/list: list of used
            codons to filter by for finding the end of the frame.
                False: uses all known startcodons

        returns

        errors
        
        """    
        # set defaults if none given
        if not startCodons:
            startCodons = ["ATG" , "CTG" , "TTG" , "GTG" , "ATT"]
        if not stopCodons:
            stopCodons = ["TAG" , "TGA" , "TAA"]

        file = open(self.filename, 'r')

        def findStart(string):
            first = -1
            for codon in startCodons:
                index = string.find(codon)
                print(codon, first, index)
                if index < first and (index != -1):
                    first = index
                elif (first == -1) and (index > -1):
                    first = index

            print('smallest index', first)
            return first

        def isStopCodon(codon):
            if codon in stopCodons:
                print("stopcodon found", codon)
                return True
      

        self.codonUsage = emptyCodonUsage.copy()

        headersCounted = 0
        framesCounted = 0
        prevline = next(file)
        frameFound = False

        # no clue if the frame is correct at this point.
        for line in file:
            
            line = line.strip() # remove '\n'
            
            if verbose: print(line)

            # find a Frame
            if not frameFound:
                if prevline.startswith('>'):
                    headersCounted+=1
                    print("found header #"+str(headersCounted))
                    prevline=line
                else:
                    window = prevline+line[:2]
                    frameStart = findStart(window)
                    if frameStart > -1:
                        startCodon = window[frameStart] + window[frameStart+1] + window[frameStart+2]
                        framesCounted +=1
                        print('found frame #'+str(framesCounted)+' with frameStart matching',startCodon,'at',frameStart)
                        frameshift = frameStart % 3
                        frameFound = True
                    # finish this by starting to count
                    for n in range(frameStart, len(window)-2, 3):
                        codon = codon = window[n]+window[n+1]+window[n+2]
                        if isStopCodon(codon):
                            frameFound = False
                        else:
                            try:
                                if 'N' in codon:
                                    self.codonUsage['TOTAL'] += 1
                                else:
                                    if verbose: print(codon)
                                    aminoacid = codon2Amino[codon]
                                    if verbose: print(aminoacid)
                                    self.codonUsage[aminoacid][aminoacid+'_TOTAL'] += 1
                                    self.codonUsage[aminoacid][codon] += 1
                                    self.codonUsage['TOTAL'] += 1
                            except KeyError:
                                self.codonUsage['TOTAL']+=1
                                self.codonUsage['ERRORS']+=1
                                print(codon,'gave error')
                            
                            
                        if verbose: print('vprint codon1:', codon, n)
                    prevline = line
                    
            # if frame found
            elif frameFound:
                if not line.startswith('>'):
                    window = prevline+line[:2]
                    for i in range(0+frameshift, len(window)-2,3):
                        codon = window[i]+window[i+1]+window[i+2]
                        if verbose: print('vprint codon2:',codon, i)
                        try:
                            if 'N' in codon:
                                self.codonUsage['TOTAL'] += 1
                            else:
                                if verbose: print(codon)
                                aminoacid = codon2Amino[codon]
                                if verbose: print(aminoacid)
                                self.codonUsage[aminoacid][aminoacid+'_TOTAL'] += 1
                                self.codonUsage[aminoacid][codon] += 1
                                self.codonUsage['TOTAL'] += 1
                        except KeyError:
                            self.codonUsage['ERRORS']+=1
                            self.codonUsage['TOTAL']+=1
                            print(codon,'gave error')
                    prevline = line
                
                else:
                    frameFound = False
                    prevline = line

        print("codon usage for",self.filename)     
        print(self.codonUsage)




if __name__ == '__main__':
    print("running from main", version)
    print("-------------------------")
    main(True) # developer mode: verbose is turned on
