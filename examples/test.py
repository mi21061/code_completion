from Bio import SeqIO

def Manhattan_distance(vector1: List[int], vector2: List[int]) -> int:
    diff = np.array(vector1) - np.array(vector2)
    distance = sum([abs(x) for x in diff])
    return distance

def Euclidean_distance(vector1: List[int], vector2: List[int]) -> int:
    diff = (np.array(vector1) - np.array(vector2))**2
    distance = np.sqrt(np.sum(diff))
    return distance

def read_data(filePath: str):
    with open(filePath, 'r') as file:
        fasta_sequences_sars1 = SeqIO.parse(file,'fasta')
    return fasta_sequences_sars1

def five_adic_codon_encoding(rna_codon: str) -> int:
    #encode RNA codon (triplet) as one number in 5-adic system
    encoding = ''
    
    for nucleotide in rna_codon:
        if nucleotide == 'C':
            encoding += '1'
        elif nucleotide == 'A':
            encoding += '2'
        elif nucleotide == 'T':
            encoding += '3'
        else: # nucleotide == G
            encoding += '4'
            
    return int(encoding)

def five_adic_codon_distance(codon1: str, codon2: str) -> int:
    encoded_codon1 = five_adic_codon_encoding(codon1)
    encoded_codon2 = five_adic_codon_encoding(codon2)
    
    distance = 0
    if int(encoded_codon1/100) != int(encoded_codon2/100):
        distance = 1

    elif int(encoded_codon1/10) != int(encoded_codon2/10):
        distance = 1/5
        
    elif encoded_codon1 % 10 != encoded_codon2 % 10:
        if abs(encoded_codon1 % 10 - encoded_codon2 % 10) == 2:
            distance = 1/2*1/25
        else:
            distance = 1*1/25
    return distance

def print_one_two_three():
    print('one')
    print('two')
    print('three')

def print_one_until_three():
    print('one')
    print('two')
    print('three')

def print_hello_world():
    print("Hello world!")


