#weektaak_3-blok_3

#IMPORTS
import os
from pathlib import Path
from matplotlib import pyplot as plot

#HARDCODES
###############################################################################
fastaname = 'deprecated'
files = ['SIVmnd2_CDS.fasta']

def main(verbose=True):
    """ criss cross bonanza"""


    # for each file in a list of filenames
    for filename in files:
        print("analysing",filename)

        # create a header:sequence dict
        headerSeqPairs = readFile(filename)
        print("made headers-seq dict ",end='')
        codonUsage = emptyCodonUsage.copy()
        print("counting codons...")

        # for each header:sequence, analyse the sequence, adding to the dictionary of counts.
        for header, seqString in headerSeqPairs.items():
            codonUsage = getCodonUsage(seqString, codonUsage)

        print("graphing dictionary for",filename)

        #finally, visualise the dictionary
        visualise(codonUsage)




def Fastafile(): # deprecated
    """ an object"""
    def __init__(self, name):
        self.name = name
        self.codonUsage = emptyCodonUsage.copy()
        self.filename = None

    def setFilename(self, filename):
        self.filename = filename


def readFile(filename):
    """ takes a file, returns a dictionary of header:string pairs"""

    file = open(filename, 'r')

    headerSeqDic = dict()
    header = next(file).strip('\n')
    sequence = ''
    
    for line in file:
        if line.startswith('>'):
            headerSeqDic[header] = sequence
            header = line.strip('\n')
            sequence = ''
        else:
            sequence += line.strip('\n')

    #append last one
    headerSeqDic[header] = sequence
    
    return headerSeqDic


# aminoacid:codon:count mapping
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


# translation dictionary
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
            
def findStart(string, verbose=True):
    """ finds the first startcodon in the whole string that
    defines the reading frame, from a list of possible startcodons

    Args:
        string: the string to parse

    returns:
        index of startcodon

    """
    startIndex = -1
    startCodon = -1
    for codon in startCodons:
        index = string.find(codon)
        if index > 0:
            print("Found", codon, 'at', index)
            if index < startIndex or startIndex == -1:
                startIndex = index
                startCodon = codon
        

    print('found earliest startcodon', startCodon,'at index =', startIndex)
    return startIndex, startCodon
            

def getCodonUsage(string, dictionary):
    """ iterates over the readingframe within a string counting all
    slices of 3 letters, assuming they are codons within reading frame

    Args
        string - the string to analyse

        dicitonary - the dictionary in which to keep the counts

    returns
        an updated dictionary
        
    errors
        Keyerror - when a key does not exist, add one to the 'ERROR' placeholder
        counter.
    """

    codonUsage = dictionary
    
    startindex, startCodon = findStart(string)

    if startindex != -1:
        print("this sequence starts with",startCodon,'at', startindex)
        aminoacid = codon2Amino[startCodon]
        codonUsage[aminoacid][aminoacid+'_TOTAL'] +=1
        codonUsage[aminoacid][startCodon] +=1
    else:
        print("No readingframe found, cause: lack of starcodon")
        return {'ERRORS':{'ERRORS_TOTAL':10000000000}} # this is to indicate horrible failures

    
    for i in range(startindex+3, len(string)-2, 3):
        codon = string[i:i+3]
        if codon in stopCodons:
            codonUsage['Stop'][codon] += 1
            codonUsage['Stop']['Stop_TOTAL'] += 1
            codonUsage['TOTAL']['TOTAL_TOTAL'] +=1
            print('found stopcodon',codon,'at',str(i))
            break
        else:
            try:
                # normal case
                aminoacid = codon2Amino[codon]
                codonUsage[aminoacid][aminoacid+'_TOTAL'] += 1
                codonUsage[aminoacid][codon] += 1
                codonUsage['TOTAL']['TOTAL_TOTAL'] += 1

            # unknown letter case
            except KeyError:
                codonUsage['ERRORS']['ERRORS_TOTAL'] +=1
                codonUsage['TOTAL']['unknownCodon'] +=1
                codonUsage['TOTAL']['TOTAL_TOTAL'] +=1
                print(codon,'at', str(i), 'gave error')
                continue
                
    print("Finished counting")

    # count totals to check if it is right

    foundTotal = codonUsage['TOTAL']['TOTAL_TOTAL']
    calcTotals = 0
    calcCodons = 0
    for key, value in codonUsage.items():
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
    return codonUsage
            
    
    
def visualise(dictionary,graphErrors=False):
        """ visualises the self.codonUsage with a sunburst diagram
        codons in the outer layer, accompanying the aminoacid they code for
        args:
            graphErrors:bool - option to toggle graphing the amount of
            found errors (unknown codons) alongside the other data

        
            """
        
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
                    #print("subkey:",codon,'\t',dictionary[aminoAcid][codon])

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
        
        plot.title('test')
        plot.axis('equal')
        plot.savefig('test'+'test'+'.png')
        plot.show()



if __name__ == '__main__':
    print("running from main")
    main()
else:
    print("imported",__name__,"successfully")
