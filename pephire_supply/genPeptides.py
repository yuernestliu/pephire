"""
Version 3.1,
1st stable release.
written by Yu Ernest Liu, 2023.02.16

Changes from v3.0: 
.genNewPips(noRepetition=True) now prohibits the generation of duplicate sequences.
.getPipPoolBook(limitLadderonSize) & getLadderonAddress(limitLadderonSize)
set the maximum scale of ladderon to prevent large ladderons from always overshadowing other sequences. 
"""
from pephire_supply import ladderpath as lp
import random

def find_all(s, sub):
    """
    Find all occurrences of a substring in a string.
    
    Args:
    s (str): The string to search within.
    sub (str): The substring to find.
    
    Returns:
    list: A list of indices where the substring is found.
    """
    len_s = len(s)
    subfoundALL = []
    i = 0
    while i < len_s:
        subfound = s.find(sub, i)
        if subfound == -1:
            return subfoundALL
        else:
            subfoundALL.append(subfound)
            i = subfound + 1


def getLadderonAddress(strs, strs_lp, limitLadderonSize=None):
    """
    Find the positions of each ladderon in the given strings.

    Args:
    strs (list of str): The list of strings to search.
    strs_lp: The ladderpath object.
    limitLadderonSize (int, optional): The maximum size of the ladderon. Defaults to None.
    
    Returns:
    dict: A dictionary mapping ladderons to their positions.
    """
    ladderonAddress = {} # Dictionary to store the positions of all ladderons
    for ladderon in strs_lp.ladderonBook.keys():
        if (limitLadderonSize is None) or (len(ladderon) <= limitLadderonSize):
            ladderonAddress[ladderon] = []
            for str0 in strs:
                for i in find_all(str0, ladderon):
                    if i not in ladderonAddress[ladderon]:
                        ladderonAddress[ladderon].append(i)
    # Process the positions of basic units
    for ch, _ in strs_lp.POM[0]:
        ladderonAddress[ch] = []
        for str0 in strs:
            for i, thisletter in enumerate(str0):
                if thisletter == ch:
                    if i not in ladderonAddress[ch]:
                        ladderonAddress[ch].append(i)
    return ladderonAddress


def genNewPeptide(peptideLen, listLadderon, listProb, LadderonAddress, disp=True):
    """
    Generate a new peptide sequence.

    Args:
    peptideLen (int): The length of the peptide.
    listLadderon (list): List of ladderons.
    listProb (list): List of probabilities corresponding to each ladderon.
    LadderonAddress (dict): Dictionary of ladderons and their positions.
    disp (bool): If True, display the peptide generation process.

    Returns:
    str: The newly generated peptide.
    """
    newPeptide = '-' * peptideLen
    NtoFill = peptideLen
    while NtoFill > 0:
        toAdd = random.choices(listLadderon, weights=listProb)[0] # Choose ladderon based on probabilities
        toPut = random.choice(LadderonAddress[toAdd]) # Choose position to insert the ladderon
        countQ = newPeptide.count('-', toPut, toPut+len(toAdd))
        if countQ > 0: # Add only if there are empty spots
            if disp:
                print(newPeptide)
            newPeptide = newPeptide[:toPut] + toAdd + newPeptide[toPut+len(toAdd):]
            NtoFill -= countQ
    return newPeptide

def genNewPips(PipPoolBook, pipPool, N=10, noRepetition=False):
    """
    Generate new peptide sequences.

    Args:
    PipPoolBook (tuple): Contains peptide information.
    pipPool (list): List of existing peptides.
    N (int): Number of new peptides to generate.
    noRepetition (bool): If True, avoid generating duplicate sequences.

    Returns:
    list: A list of new peptide sequences.
    """
    newpips = [] 
    for _ in range(N):
        temp = genNewPeptide(PipPoolBook[0], PipPoolBook[1], PipPoolBook[2], PipPoolBook[3], disp=False)
        if noRepetition:
            # Ensure no duplicate sequences are generated
            if temp not in newpips and temp not in pipPool:
                newpips.append(temp)
        else:
            newpips.append(temp)
    return newpips

def getPipPoolBook(pipPool, limitLadderonSize=None):
    """
    Get the pip pool book.

    Args:
    pipPool (list): The pool of pips.
    limitLadderonSize (int, optional): The maximum size of the ladderon. Defaults to None.

    Returns:
    tuple: A tuple containing information about the peptide, ladderons, their probabilities, and positions.
    """
    peptideLen = len(pipPool[0])
    strs_lp = lp.ladderpath(pipPool, CalPOM=True)
    LadderonAddress = getLadderonAddress(pipPool, strs_lp, limitLadderonSize=limitLadderonSize)

    listLadderon, listProb = [], []
    for pom in strs_lp.POM:
        for ladderon, multi in pom:
            if (limitLadderonSize is None) or (len(ladderon) <= limitLadderonSize):
                listLadderon.append(ladderon)
                # Calculate the frequency of each ladderon being chosen
                listProb.append(multi*len(ladderon))
            
    return (peptideLen, listLadderon, listProb, LadderonAddress)
    # listLadderon: the list of all ladderons, get from pipPool. ['W', 'QL', 'RLA'...]
    # listProb: the probability being taken for new pip, by !!! user defined !!! method. [2,6,8...]
    # LadderonAddress: the position of each ladderon can be. {'AGDEFE': [11], 'RIGDE': [10], ...}
