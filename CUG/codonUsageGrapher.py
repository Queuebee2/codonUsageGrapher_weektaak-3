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
# 
# 6-3-2018
# 
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


#hardcoded
# delimiter = '_'
# 

def main(verbose=False):

    # yes, it's hardcoded. 
    sequencesPath = findUp('sequences')
    print("Found presumable path with sequences")
    print(sequencesPath)
    spacer()
    setCWD(sequencesPath)
    
    availableFastas =  getFilenames()
    fastaObjects = createFastaobjects(availableFastas, 'CDS')
    

aa3 = {"Ala ": ["GCT" , "GCC" , "GCA" , "GCG"] ,
       "Arg": ["CGT" , "CGC" , "CGA" , "CGG" , "AGA" , "AGG"] ,
       "Asn": ["AAT" , "AAC"] ,
       "Asp": ["GAT" , "GAC"] ,
       "Cys": ["TGT" , "TGC"] ,
       "Gln": ["CAA" , "CAG"] ,
       "GlT": ["GAA" , "GAG"] ,
       "Gly": ["GGT" , "GGC" , "GGA" , "GGG"] ,
       "His": ["CAT" , "CAC"] ,
       "Ile": ["ATT" , "ATC" , "ATA"] ,
       "LeT": ["TTA" , "TTG" , "CTT" , "CTC" , "CTA" , "CTG"] ,
       "Lys": ["AAA" , "AAG"] ,
       "Met": ["ATG"] ,
       "Phe": ["TTT" , "TTC"] ,
       "Pro": ["CCT" , "CCC" , "CCA" , "CCG"] ,
       "Ser": ["TCT" , "TCC" , "TCA" , "TCG" , "AGT" ,"AGC"] ,
       "Thr": ["ACT" , "ACC" , "ACA" , "ACG"] ,
       "Trp": ["TGG"] ,
       "Tyr": ["TAT" , "TAC"] ,
       "Val": ["GTT" , "GTC" , "GTA" , "GTG"] ,
       "Start": ["ATG" , "CTG" , "TTG" , "GTG" , "ATT"] ,
       "Stop" : ["TAG" , "TGA" , "TAA"]}

        
def createFastaobjects(filesToParse, keyword=None, delimiter='_' ,verbose=False,):
    createdStuff = [] # list of objects to return
    
    for fastaFilename in filesToParse:
        parts = fastaFilename.strip('.fasta').split(delimiter)
        name = parts[0]
        description = parts[1]
        if not keyword:
            print('creating Fastafile object with name:',name,'\tdesc:',description)
            fastaObject = Fastafile(name, description)
            createdStuff.append(fastaObject)
        elif keyword in name or keyword in description:
            print('creating Fastafile object with name:',name,'\tdesc:',description)
            fastaObject = Fastafile(name, description)
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
    print("working in\"", os.getcwd(), "\"now")
    spacer()
    
def findUp(searchterm,thispath=False, depth=3,verbose=False):
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
        return "ran out of depth to go through"
    
    print('looking for',searchterm,'in',thispath)
    for dirname, nested_dirnames, absolute_filenames in os.walk(thispath):
        if dirname.split(oschar)[-1] == searchterm:
            result = str(dirname)
            if verbose: print("\nfindUp found",result)
            # return filepath in string format
            return result
        
    if verbose: print('findUp going up a directory!')
    return findUp(str(Path(thispath).parent), searchterm, depth-1)


class Fastafile():
    """ docstring """
    def __init__(self, name, desc):
        self.name = name
        self.desc = desc
        self.usedCodon








if __name__ == '__main__':
    print("running from main", version)
    print("-------------------------")
    main()
