#rewriting the code with a different attac.

#hardcoded name for test purposes

fastaname = 'SIVmnd2_CDS.fasta'

def main(verbose=True):
    """ criss cross bonanza"""

    headerSeqPairs = readFile(fastaname)
    if verbose:
        for header, seqString in headerSeqPairs.items():
            print(header[:10], seqString[:10])
            codonUsage = getCodons(seqString)


    # get filenames (read path)
    # list filenames
    # create object for each filename (each organism)
    # create all header-



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
            
def findStart(string, verbose=True):
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
            

def getCodonUsage(self, string):
    """ searches for the reading frame by finding a startCodon
    keeps iterating over each following codon until a stopcodon is found
    or the end of the line is reached
    while iterating, keeps a dictionary with counts for each aminoacid and
    codon
    returns a dictionary with {AminoAcid:{Codon:Count}} nested dictionaries
    """
    
    
    startindex, startCodon = findStart(string)

    if startindex != -1:
        print("this sequence starts with",startCodon,'at', startindex)
        aminoacid = codon2Amino[startCodon]
        codonUsage[aminoacid][aminoacid+'_TOTAL'] +=1
        codonUsage[aminoacid][startCodon] +=1
    else:
        print("No readingframe found, cause: lack of starcodon")
        return {'ERROR':{'ERROR':1000000}}

    
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
                else:
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
    print("Totals of each: "end='')
    print(foundTotal, calcTotals, calcCodons)
    if foundTotal == calcTotals and calcTotals == calcCodons:
        print("It seems like it all counts up!")
    return codonUsage
            
    
    

if __name__ == '__main__':
    print("running from main")
    main()
else:
    print("imported",__name__,"successfully")
