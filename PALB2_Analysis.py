import re
import pandas as pd

"""
This is a python script that analyse FASTQ file of 20000 PALB2 variants library
The script classify which variant is a forward and reverse sequence based on the attB 
position and then identify the reading sequences, barcodes, and twist barcodes.
Finally, the script identifies which mutation are at which indices and bases on both the variant
and the PALB2 sequences, and then names which amino acids are mutated
All the data are exported into a csv
"""

# create an empty dataFrame object to store the information of the sequences
# 'orginal_index' refers to the index on the library variant (the entire sequence)
# 'mutation_index' refers to the index on the PABL2 gene only (not the entire sequence)
df = pd.DataFrame(columns=['Variant', 'Sequence', 'Orientation', 'attB_index',
                           'read_1', 'read_2', 'barcode', 'barcode_length', 'bar_type',
                           'twist1', 'twist2', 'twist_barcode', 'orginal_index',
                           'mutation_index', 'mutation', 'AA_mutationa'])

# The original PALB2 sequence that is used to find mutation
original_PALB2 = "ATGGACGAGCCTCCCGGGAAGCCCCTCAGCTGTGAGGAGAAGGAAAAGTTAAAGGAGAAACTTGCTTTTTTGAAAAGGGAATACAGCAAGACACTAGCCCGCCTTCAGCGTGCCCAAAGAGCTGAAAAGATTAAGCATTCTATTAAGAAAACAGTAGAAGAACAAGATTGTTTGTCTCAGCAGGATCTCTCACCGCAGCTAAAACACTCAGAACCTAAAAATAAAATATGTGTTTATGACAAGTTACACATCAAAACCCATCTTGATGAAGAAACTGGAGAAAAGACATCTATCACACTTGATGTTGGGCCTGAGTCCTTTAACCCTGGAGATGGCCCAGGAGGATTACCTATACAAAGAACAGATGACACCCAAGAACATTTTCCCCACAGGGTCAGTGACCCTAGTGGTGAGCAAAAGCAGAAGCTGCCAAGCAGAAGAAAGAAGCAGCAGAAGAGGACATTTATTTCACAGGAGAGAGACTGTGTCTTTGGCACTGATTCACTCAGATTGTCTGGGAAAAGACTAAAGGAACAGGAAGAAATCAGTAGCAAAAATCCTGCTAGATCACCAGTAACTGAAATAAGAACTCACCTTTTAAGTCTTAAATCTGAACTTCCAGATTCTCCAGAACCAGTTACAGAAATTAATGAAGACAGTGTATTAATTCCACCAACTGCCCAACCAGAAAAAGGTGTTGATACATTCCTAAGAAGACCTAATTTCACCAGGGCGACTACAGTTCCTTTACAGACTCTATCAGATAGCGGTAGTAGTCAGCACCTTGAACACATTCCTCCTAAAGGTAGCAGTGAACTTACTACTCACGACCTAAAAAACATTAGATTTACTTCACCTGTAAGTTTGGAGGCACAAGGCAAAAAAATGACTGTCTCTACAGATAACCTCCTTGTAAATAAAGCTATAAGTAAAAGTGGCCAACTGCCCACAAGTTCTAATTTAGAGGCAAATATTTCATGTTCTCTAAATGAACTCACCTACAATAACTTACCAGCAAATGAAAACCAAAACTTAAAAGAACAAAATCAAACAGAGAAATCTTTAAAATCTCCCAGTGACACTCTTGATGGCAGGAATGAAAATCTTCAGGAAAGTGAGATTCTAAGTCAACCTAAGAGTCTTAGCCTGGAAGCAACCTCTCCTCTTTCTGCAGAAAAACATTCTTGCACAGTGCCTGAAGGCCTTCTGTTTCCTGCAGAATATTATGTTAGAACAACACGAAGCATGTCCAATTGCCAGAGGAAAGTAGCCGTGGAGGCTGTCATTCAGAGTCATTTGGATGTCAAGAAAAAAGGGTTTAAAAATAAAAATAAGGATGCAAGTAAAAATTTAAACCTTTCCAATGAGGAAACTGACCAAAGTGAAATTAGGATGTCTGGCACATGCACAGGACAACCAAGTTCAAGAACCTCTCAGAAACTTCTCTCATTAACTAAAGTCAGCTCTCCCGCTGGGCCCACTGAAGATAATGACTTGTCTAGGAAGGCAGTTGCCCAAGCACCTGGTAGAAGATACACAGGAAAAAGAAAATCAGCCTGCACCCCAGCATCAGATCATTGTGAACCACTTTTGCCAACTTCTAGCCTGTCGATTGTTAACAGGTCCAAGGAAGAAGTCACCTCACACAAATATCAGCACGAAAAATTATTTATTCAAGTGAAAGGGAAGAAAAGTCGTCATCAAAAAGAGGATTCCCTTTCTTGGAGTAATAGTGCTTATTTATCCTTGGATGATGATGCTTTCACGGCTCCATTTCATAGGGATGGAATGCTGAGTTTAAAGCAACTACTGTCTTTTCTCAGTATCACAGACTTTCAGTTACCTGATGAAGACTTTGGACCTCTTAAGCTTGAAAAAGTGAAGTCCTGCTCAGAAAAACCAGTGGAGCCCTTTGAGTCAAAAATGTTTGGAGAGAGACATCTTAAAGAGGGAAGCTGTATTTTTCCAGAGGAACTGAGTCCTAAACGCATGGATACAGAAATGGAGGACTTAGAGGAAGATCTTATTGTTCTACCAGGAAAATCACATCCCAAAAGGCCAAACTCGCAAAGCCAGCATACAAAGACGGGCCTTTCTTCATCCATATTACTTTATACTCCTTTAAATACGGTTGCGCCTGATGATAATGACAGGCCTACCACAGACATGTGTTCACCTGCTTTCCCCATCTTAGGTACTACTCCAGCCTTTGGCCCTCAAGGCTCCTATGAAAAAGCATCTACAGAAGTTGCTGGACGAACTTGCTGCACACCCCAACTTGCTCATTTGAAAGACTCAGTCTGTCTTGCCAGTGATACTAAACAATTCGACAGTTCAGGCAGCCCAGCAAAACCACATACCACCCTGCAAGTGTCAGGCAGGCAAGGACAACCTACCTGTGACTGTGACTCTGTCCCGCCAGGAACACCTCCACCCATTGAGTCATTCACTTTTAAAGAAAATCAGCTCTGTAGAAACACATGCCAGGAGCTGCATAAACATTCCGTCGAACAGACTGAAACAGCAGAGCTTCCTGCTTCTGATAGCATAAACCCAGGCAACCTACAATTGGTTTCAGAGTTAAAGAATCCTTCAGGTTCCTGTTCCGTAGATGTGAGTGCCATGTTTTGGGAAAGAGCCGGTTGTAAAGAGCCATGTATCATAACTGCTTGCGAAGATGTAGTTTCTCTTTGGAAAGCTCTGGATGCTTGGCAGTGGGAAAAACTTTATACCTGGCACTTCGCAGAGGTTCCAGTATTACAGATAGTTCCAGTGCCTGATGTGTATAATCTCGTGTGTGTAGCTTTGGGAAATTTGGAAATCAGAGAGATCAGGGCATTGTTTTGTTCCTCTGATGATGAAAGTGAAAAGCAAGTACTACTGAAGTCTGGAAATATAAAAGCTGTGCTTGGCCTGACAAAGAGGAGGCTAGTTAGTAGCAGTGGGACCCTTTCTGATCAACAAGTAGAAGTCATGACGTTTGCAGAAGATGGAGGAGGCAAAGAAAACCAATTTTTGATGCCCCCTGAGGAGACTATACTAACTTTTGCTGAGGTCCAAGGGATGCAAGAAGCTCTGCTTGGTACTACTATTATGAACAACATTGTTATTTGGAATTTAAAAACTGGTCAACTCCTGAAAAAGATGCACATTGATGATTCTTACCAAGCTTCAGTCTGTCACAAAGCCTATTCTGAAATGGGGCTTCTCTTTATTGTCCTGAGTCATCCCTGTGCCAAAGAGAGTGAGTCGTTGCGAAGCCCTGTGTTTCAGCTCATTGTGATTAACCCTAAGACGACTCTCAGCGTGGGTGTGATGCTGTACTGTCTTCCTCCAGGGCAGGCTGGCAGGTTCCTGGAAGGTGACGTGAAAGATCACTGTGCAGCAGCAATCTTGACTTCTGGAACAATTGCCATTTGGGACTTACTTCTCGGTCAGTGTACTGCCCTCCTCCCACCTGTCTCTGACCAACATTGGTCTTTTGTGAAATGGTCGGGTACAGACTCTCATTTGCTGGCTGGACAAAAAGATGGAAATATATTTGTATACCACTATTCATAA"


# Search for the attB sequence and return its indices
# If forward strand, then attB is forward and vice versea
def attB_match(search_seq):
    attB_forward = "CCGGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCATCC"
    attB_reverse = "GGATGATCCTGACGACGGAGACCGCCGTCGTCGACAAGCCGG"
    forward_match = re.search(attB_forward, search_seq)
    reverse_match = re.search(attB_reverse, search_seq)
    if forward_match:
        return ["forward", forward_match.start()]
    if reverse_match:
        return ["reverse", reverse_match.start()]
    else:
        return ["N/A", "N/A"]


# Search for both the reading 1 and 2 sequences and return their indices
# The reading sequence indices are used to determine the barcode in between
def reading_match(search_seq):
    read_1 = "ACAAATAGTT"
    read_2 = "TGCGAGTAGT"
    read1_index = 0
    read2_index = 0
    read1_match = re.search(read_1, search_seq)
    read2_match = re.search(read_2, search_seq)
    if read1_match:
        read1_index = read1_match.start()
    else:
        read1_index = "N/A"
    if read2_match:
        read2_index = read2_match.start()
    else:
        read2_index = "N/A"
    return [read1_index, read2_index]


# Find the barcode sequence in between the reading 1 and 2 sequences and its length
def find_barcode(sequence, read_1, read_2):
    if isinstance(read_1, int) and isinstance(read_2, int):
        return seq[read_1 + 10:read_2]


# Return the length of the bar code
def bar_length(bar):
    try:
        return len(bar)
    except:
        return "N/A"


# Return the types of each barcode: normal, indel, or unbar
def bar_type(length):
    try:
        if length == 18:
            return "normal"
        if length <= 6:
            return "unbar"
        else:
            return "indel"
    except:
        return "N/A"


# Search for both the twist reading 1 and 2 sequences and return their indices
# The twist reading sequence indices are used to determine the twist barcode in between
def twist_match(search_seq):
    twist_1 = "AGAGTCAGATCAGTCTCGAG"
    twist_2 = "GAATTCCTGACCTCCTTCTC"
    twist1_index = 0
    twist2_index = 0
    twist1_match = re.search(twist_1, search_seq)
    twist2_match = re.search(twist_2, search_seq)
    if twist1_match:
        twist1_index = twist1_match.start()
    else:
        twist1_index = "N/A"
    if twist2_match:
        twist2_index = twist2_match.start()
    else:
        twist2_index = "N/A"
    return [twist1_index, twist2_index]


# Find the twist barcode sequence in between the twist reading 1 and 2 sequences and its length
def find_twist(sequence, twist_1, twist_2):
    if isinstance(twist_1, int) and isinstance(twist_2, int):
        return seq[twist_1 + 20:twist_2]


# Find the start amino acid of the PABL2 sequence to find the mutation
def find_startAA(search_seq):
    startAA = "ATGGACGAGCCT"
    startAA_match = re.search(startAA, search_seq)
    if startAA_match:
        return startAA_match.start()
    else:
        return "N/A"


# Iterate through the PABL2 sequence and compare it with the input sequence
# Only find the mutation if the twist bar and the start amino acid of PABL2 are present
# The indices of the original sequqnce and the PABL sequence are recorded
# The change in bases is also recorded
def findMutation(seq, twist_bar):
    if twist_bar is None:
        return ["N/A", "N/A", "N/A"]
    original_index = []
    mutation_index = []
    mutation_map = []
    j = 0
    startIndex = find_startAA(seq)
    if isinstance(startIndex, int):
        for i in range(startIndex, len(original_PALB2)):
            try:
                if seq[i] != original_PALB2[j]:
                    original_index.append(j)
                    mutation_index.append(i)
                    mutation_map.append(str(original_PALB2[j]) + ">" + str(seq[i]))
            except:
                return "N/A"
            j = j + 1
        return [original_index, mutation_index, mutation_map]
    else:
        return ["N/A", "N/A", "N/A"]


# Translate a codon of 3 bases into the appropriate amino acids
def amino_acids(aa):
    if aa == "TTT" or aa == "TTC":
        return "Phe"
    if aa == "TTA" or aa == "TTG" or aa == "CTT" or aa == "CTC" or aa == "CTA" or aa == "CTG":
        return "Leu"
    if aa == "ATT" or aa == "ATC" or aa == "ATA":
        return "Ile"
    if aa == "ATG":
        return "Met"
    if aa == "GTT" or aa == "GTC" or aa == "GTA" or aa == "GTG":
        return "Val"
    if aa == "TCT" or aa == "TCC" or aa == "TCA" or aa == "TCG" or aa == "AGT" or aa == "AGC":
        return "Ser"
    if aa == "CCT" or aa == "CCC" or aa == "CCA" or aa == "CCG":
        return "Pro"
    if aa == "ACT" or aa == "ACC" or aa == "ACA" or aa == "ACG":
        return "Thr"
    if aa == "GCT" or aa == "GCC" or aa == "GCA" or aa == "GCG":
        return "Ala"
    if aa == "TAT" or aa == "TAC":
        return "Tyr"
    if aa == "TAA" or aa == "TAG" or aa == "TGA":
        return "Stop"
    if aa == "CAT" or aa == "CAC":
        return "His"
    if aa == "CAA" or aa == "CAG":
        return "Gln"
    if aa == "AAT" or aa == "AAC":
        return "Asn"
    if aa == "AAA" or aa == "AAG":
        return "Lys"
    if aa == "GAT" or aa == "GAC":
        return "Asp"
    if aa == "GAA" or aa == "GAG":
        return "Glu"
    if aa == "TGT" or aa == "TGA":
        return "Cys"
    if aa == "TGG":
        return "Trp"
    if aa == "CGT" or aa == "CGC" or aa == "CGA" or aa == "CGG" or aa == "AGA" or aa == "AGG":
        return "Arg"
    if aa == "GGT" or aa == "GGC" or aa == "GGA" or aa == "GGG":
        return "Gly"

      
      
# A helper method to find which mutation index goes with other indices as part of a codon
# For example, a index of 1 will have 0 and 2 as part of the same codon
def helper_findSameCodonIndex(index):
    if index % 3 == 0:
        return [index + 1, index + 2]
    if index % 3 == 1:
        return [index - 1, index + 1]
    if index % 3 == 2:
        return [index - 1, index - 2]


# A helper method that uses the helper_takeSameCodon_2
# This returns a dictionary that list the each indices that belong to the same codon
# For example, AA1': [687, 689], 'AA2': [1602, 1603, 1604], 'AA3': [1605, 1606]
def helper_takeSameCodon_1(mutation_index):
    if mutation_index == "N/A":
        return "N/A"
    if isinstance(mutation_index, str):
        return "N/A"
    mutation_index_copy = mutation_index.copy()
    aa_map = {}
    count = 1
    helper_takeSameCodon_2(mutation_index_copy, aa_map, count)
    return aa_map


# A recursion helper method that uses helper_findSameCodonIndex
def helper_takeSameCodon_2(mutation_index_copy, aa_map, count):
    if len(mutation_index_copy) == 0:
        return
    taken_index = []
    # Use the helper method to find all possible indices that is in the same codon
    findSameCodonIndex = helper_findSameCodonIndex(mutation_index_copy[0])
    # Append the taken indices and then remove it
    taken_index.append(mutation_index_copy[0])
    mutation_index_copy.remove(mutation_index_copy[0])
    # Find the indices that is in the same codon as the the taken index
    for i in range(0, len(mutation_index_copy)):
        if mutation_index_copy[i] in findSameCodonIndex:
            taken_index.append(mutation_index_copy[i])
    # Remove all taken indices from the list of remaining mutation indices
    mutation_index_copy = [i for i in mutation_index_copy if i not in taken_index]
    aa_map["AA" + str(count)] = taken_index
    count = count + 1
    # Recursion
    helper_takeSameCodon_2(mutation_index_copy, aa_map, count)


# Return the remaining indices of the same codon given a list of indices
# For example, input [0,2] will return [0,1,2]
def helper_getCodonIndex(index_list):
    if len(index_list) == 3:
        return index_list
    if len(index_list) == 2:
        if index_list[0] % 3 == 0 and index_list[1] % 3 == 1:
            return [index_list[0], index_list[1], index_list[1] + 1]
        if index_list[0] % 3 == 1 and index_list[1] % 3 == 2:
            return [index_list[0] - 1, index_list[0], index_list[1]]
        if index_list[0] % 3 == 0 and index_list[1] % 3 == 2:
            return [index_list[0], index_list[0] + 1, index_list[1]]
    if len(index_list) == 1:
        if index_list[0] % 3 == 0:
            return [index_list[0], index_list[0] + 1, index_list[0] + 2]
        if index_list[0] % 3 == 1:
            return [index_list[0] - 1, index_list[0], index_list[0] + 1]
        if index_list[0] % 3 == 2:
            return [index_list[0] - 2, index_list[0] - 1, index_list[0]]


# This concantenate a codon of bases from the index list
def helper_getCodonBase(index_list, variant):
    codon_base = ""
    for i in index_list:
        codon_base = codon_base + str(variant[i])
    return codon_base


# Return a dictionary that show which amino acids are mutated and which mutation it is
# For example, return AA1': 'Gln>Ser', 'AA2': 'Thr>Leu', 'AA3': 'Gly>Tyr'
def translate_AA(original_map, mutation_map, seq):
    if mutation_map == "N/A" or original_map == "N/A":
        return "N/A"
    final_map = {}
    # Iterate through both the dictionary of codons in the variant as well as the PALB2 gene
    for (k, v), (k2, v2) in zip(original_map.items(), mutation_map.items()):
        # Get all 3 bases indices of the codon
        get_OriCodon = helper_getCodonIndex(v)
        get_MutCodon = helper_getCodonIndex(v2)
        # Get the bases of the codon
        ori_base = helper_getCodonBase(get_OriCodon, seq)
        mut_base = helper_getCodonBase(get_MutCodon, original_PALB2)
        # Get the amino acids that the codon translate to
        ori_aa = amino_acids(ori_base)
        mut_aa = amino_acids(mut_base)
        # Update the dictionary
        final_map[k] = str(ori_aa) + ">" + str(mut_aa)
    return final_map


# Add row to the dataframe by iterating through the fastq file
with open('Sequel.RunS179-2_S2.002.PALB2lib.ccs.fastq', 'r') as fh:
    count = 1
    lines = []
    for line in fh:
        lines.append(line.rstrip())
        # Only choose lines with length 4 which is how fastq file structures
        if len(lines) == 4:
            # Get the sequence
            seq = lines[1]
            # Find attB
            attB = attB_match(seq)
            # If the sequence is in reverse orientation, then make it forward
            if attB[0] == "reverse":
                seq = seq[::-1]
            # Find the reading 1 and 2 sequence
            reading = reading_match(seq)
            # Find the bar code and its length
            barcode = find_barcode(seq, reading[0], reading[1])
            length = bar_length(barcode)
            # Determine the bar type
            bartype = bar_type(length)
            # Find the twist reading sequence
            twist_read = twist_match(seq)
            # Find the twist bar code
            twist_bar = find_twist(seq, twist_read[0], twist_read[1])
            # Find the index mutation
            mutation_read = findMutation(seq, twist_bar)
            # Find the amino acid mutation
            ## Get the sequence and the PALB2 mutation index
            original_index = mutation_read[0]
            mutation_index = mutation_read[1]
            ## Get the full codon of 3 bases to be translated
            original_map = helper_takeSameCodon_1(original_index)
            mutation_map = helper_takeSameCodon_1(mutation_index)
            ## Translate the codon 3 bases into amino acid
            aa_mutation = translate_AA(original_map, mutation_map, seq)

            # Append row of sequence information above
            df = df.append({'Variant': count, 'Sequence': seq, 'Orientation': attB[0],
                            'attB_index': attB[1], 'read_1': reading[0], 'read_2': reading[1],
                            'barcode': barcode, 'barcode_length': length, 'bar_type': bartype,
                            'twist1': twist_read[0], 'twist2': twist_read[1], 'twist_barcode': twist_bar,
                            'orginal_index': mutation_read[0], 'mutation_index': mutation_read[1],
                            'mutation': mutation_read[2], 'AA_mutationa': aa_mutation},
                           ignore_index=True)

            # Tracking the progress of appending at each 1000 rows by printing
            if len(df.index) % 1000 == 0:
                print(len(df.index))
            lines = []
            count = count + 1

# Export dataframe to csv
df.to_csv(r'PALB2_Analysis.csv', index=False, header=True)
