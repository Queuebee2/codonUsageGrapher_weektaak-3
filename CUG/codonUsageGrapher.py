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
# dont try to open pictures as fasta -_-
# grapher function for dictionaries
# and codon2Amino to translate
# added dictionaries 'emptyData' to keep track of codonUsage
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
from matplotlib import pyplot as plot


# hardcoded globals
# delimiter = '_'
#
pathKeyword = 'weektaak_3_relevant'
graphKeyword = 'newgraphs'
recursionDepth = 3

                    

def main(verbosity=False):
    """
    """
    userChoice = input("run automatically? y/n :")
    if userChoice in 'NNoonnoo':
        print("please provide the directory path where to-be-analysed sequence fasta files reside")
        sequencesDir = input()
    else:
        sequencesDir= findUp(pathKeyword)

        
    # find fasta files
    
    print("Found presumable path with sequences")
    print(sequencesDir)
    spacer()
    setCWD(sequencesDir)
    availableFiles =  getFilenames() # lists files in cwd

    # create objects to hold information
    objectsToGraph = createFastaObjects(availableFiles)

    """    if verbosity:
        print("printing object attributes")
        for o in objectsToGraph:
            print('printing name:',o.name)
            o.printSummary()
    """      

    for o in objectsToGraph:
        print('creating codonUsage dicitonary for',o.name)
        o.readFile()
        o.makeCodonUsage()
        print("deleting header-sequence pairs now...")
        o.cleanup()
        print('putting together graph...')
        o.visualise()


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

emptyCodonUsage = {'TOTAL': {'TOTAL_TOTAL':0,'unknownCodon':0},
                   'ERRORS':{'ERRORS_TOTAL':0},
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

startCodons = ["ATG" , "CTG" , "TTG" , "GTG" , "ATT"]
stopCodons = ["TAG" , "TGA" , "TAA"]
        
def createFastaObjects(listOfFilepaths, keyword=None, special=None,delimiter='_' ,verbose=False,):
    """

    changelog
    - removed most parsing features
    - change: filesToParse -> listOfFilepaths
    - deprecated: keyword, special because given directory is assumed to be filled
      only with correct fastas, so none have to be filtered out
    
    args
        listOfFilepaths:list - a list with all the paths to files in a given directory
        
        (deprecated)keyword: a keyword that has to be in the name

        (deprecated)special: a list of strings containing exact names (for hardcoding purposes
        if parsing is too complicated)

        delimiter:default= '_' - assumed to be used as separator for name and description
        of the given filename


    """
    fastaObjects = [] # list of objects to return

    count = 0
    for fastaFilename in listOfFilepaths:
        parts = fastaFilename.strip('.fasta').split(delimiter)
        name = parts[0]
        description = parts[1]

        print("creating fastaObject:",end='')
        print(name, description, fastaFilename)
        FastaObject = Fastafile(name, description, fastaFilename)
        count+=1
        fastaObjects.append(FastaObject)
        

    print('created',str(count),'fasta objects successfully')
    spacer()
    return fastaObjects  


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
        self.headersParsed = 0
        self.headersAndSequences = None

    
    def readFile(self):
        """ takes a file, returns a dictionary of header:string pairs"""
        print('creating headers and sequence dict for',self.name)
        file = open(self.filename, 'r')

        headerSeqDic = dict()
        header = next(file).strip('\n')
        sequence = ''
        headercount = 1
        linecount = 0
        for line in file:
            linecount+=1
            if linecount % 100000 ==0:
                print('line '+str(linecount))
            if line.startswith('>'):
                headercount+=1
                headerSeqDic[header] = sequence
                header = line.strip('\n')
                sequence = ''
            else:
                sequence += line.strip('\n')

        #append last one
        headerSeqDic[header] = sequence

        print('created',headercount,'header-sequence pairs for',self.name)
        file.close()
        self.HeadersParsed = headercount
        self.headersAndSequences = headerSeqDic
        
    def cleanup(self):
        self.headersAndSequences = None
        
    def getName(self):
        return self.name

    def getDescription(self):
        return self.name

    
    def getFilename(self):
        return self.filename

    def getCodonUsage(self):
        return self.codonUsage
        
    def printSummary(self):
        print('name:',self.name,'\tdesc:',self.description,\
              '\thasfilename:', str(bool(self.filename)))

    def visualise(self, graphErrors=True):
        """ visualises the self.codonUsage with a sunburst diagram
        args:
            graphErrors:bool - option to toggle graphing the amount of
            found errors (unknown codons) alongside the other data
            """

        dictionary = self.codonUsage
        
        outer_labels =[]
        inner_labels = []
        outerPie = []
        innerPie = []

        for aminoAcid in sorted(dictionary.keys()):
            # base cases
            if aminoAcid == 'TOTAL':
                continue
            elif aminoAcid == 'ERRORS':
                """
                if graphErrors:
                    # create error slice with unknown codons
                    outer_labels.append('unknown Codons')
                    outerPie.append(dictionary['TOTAL']['unknownCodon'])

                    inner_labels.append("unknown codons")
                    innerPie.append(dictionary['ERRORS']['ERRORS_TOTAL'])
                """
                continue
            else:
                for codon in sorted(dictionary[aminoAcid].keys()):
                    # skip totals, we only use that in the inner piechart.
                    if '_TOTAL' in codon:
                        continue
                    print("subkey:",codon,'\t',dictionary[aminoAcid][codon])

                    if dictionary[aminoAcid][codon] >0:
                    # create outer slice for codon
                        outer_labels.append(codon)
                        outerPie.append(dictionary[aminoAcid][codon])

                # create inner pieslice for aminoAcid
                if dictionary[aminoAcid][aminoAcid+'_TOTAL'] >0:
                    inner_labels.append(aminoAcid)
                    innerPie.append(dictionary[aminoAcid][aminoAcid+'_TOTAL'])



        plot.pie(outerPie, labels=outer_labels, radius=10,labeldistance=0.8, counterclock=False, rotatelabels=True)
        plot.pie(innerPie, labels=inner_labels, radius=7, labeldistance=0.7,counterclock=False, rotatelabels=True)
        
        plot.title(self.name+self.description)
        plot.axis('equal')
        plot.savefig(self.name+self.description+'.png')
        plot.show()



    def makeCodonUsage(self):
        for header, sequence in self.headersAndSequences.items():
            print('complementing dictionary with',header[:15])
            self.createDictionary(sequence)

    def findStart(self, string, verbose=True): # deprecated
        """ finds the first startcodon in the whole string that
        defines the reading frame, from a list of possible startcodons

        Args:
            string: the string to parse

        returns:
            index of startcodon

        """
        first = -1
        for codon in startCodons:
            index = string.find(codon)
            if verbose: print(codon, first, index)
            if index < first and (index != -1):
                first = index
                foundCodon = codon
            elif (first == -1) and (index > -1):
                first = index
                foundCodon = codon

        if index != -1:
            print('found startcodon', foundCodon,'at pos',index)
            print('that looks like this')
            print(string[:index+4])
            return index, codon
        else:
            print('no startcodon found')
            return -1, 'NAN'
            
    def createDictionary(self, string, verbose=False):
        """
        
        """
        
        
        startindex = string.find('ATG')
        # only looking for ATG anymore

        if startindex != -1:
            print("this sequence starts at", startindex)
            #nvm counting them makes the graph obnoxious!
            #aminoacid = codon2Amino[startCodon]
            #self.codonUsage[aminoacid][aminoacid+'_TOTAL'] +=1
            #self.codonUsage[aminoacid][startCodon] +=1
        else:
            print("No readingframe found, cause: lack of starcodon")
            return {'ERROR':{'ERROR':1000000}}

        
        for i in range(startindex+3, len(string)-2, 3):
            codon = string[i:i+3]
            if codon in stopCodons:
                #self.codonUsage['Stop'][codon] += 1
                #self.codonUsage['Stop']['Stop_TOTAL'] += 1
                #self.codonUsage['TOTAL']['TOTAL_TOTAL'] +=1
                print('found stopcodon',codon,'at',str(i))
                break
            else:
                try:
                    # normal case

                    aminoacid = codon2Amino[codon]
                    self.codonUsage[aminoacid][aminoacid+'_TOTAL'] += 1
                    self.codonUsage[aminoacid][codon] += 1
                    self.codonUsage['TOTAL']['TOTAL_TOTAL'] += 1

                # unknown letter case
                except KeyError:
                    self.codonUsage['ERRORS']['ERRORS_TOTAL'] +=1
                    self.codonUsage['TOTAL']['unknownCodon'] +=1
                    self.codonUsage['TOTAL']['TOTAL_TOTAL'] +=1
                    print(codon,'at', str(i), 'gave error')
                    continue
                    
        print("Finished counting")

        # count totals to check if it is right

        foundTotal = self.codonUsage['TOTAL']['TOTAL_TOTAL']
        calcTotals = 0
        calcCodons = 0
        for key, value in self.codonUsage.items():
            for subkey, subvalue in value.items():
                if 'TOTAL' in subkey:
                    if subkey != 'TOTAL_TOTAL':
                        calcTotals += subvalue
                
                else:
                    calcCodons += subvalue
        print("Totals of each: ",end='')
        print(foundTotal, calcTotals, calcCodons)
        if foundTotal == calcTotals and calcTotals == calcCodons:
            print("It seems like it all counts up!")

        print(self.codonUsage)
        return self.codonUsage




if __name__ == '__main__':
    print("running from main", version)
    print("-------------------------")
    main(True) # developer mode: verbose is turned on
