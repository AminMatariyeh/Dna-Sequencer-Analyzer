def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines() 

    sequence =""
    for line in lines:
        if line.startswith(">"): 
            continue  
        sequence +=line.strip().upper() 
    
    
    return sequence

if __name__ == "__main__":


    dna = read_fasta("example.fasta")
    

def count_bases(dna):
    counts = {"A":dna.count("A"),"T":dna.count("T"),"G":dna.count("G"),"C":dna.count("C")}
    return counts
def gc_content(dna):
    g=dna.count("G")
    c=dna.count("C")
    gc_percent= ((g+c)/(len(dna)))*100
    return (gc_percent)
def reverse_complement(dna):
    complement = {"A":"T", "T": "A", "G": "C", "C": "G"}
    reverse_comp= "".join(complement[base] for base in reversed(dna))
    return reverse_comp

codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
    
def translate_dna_to_protein(dna_seq):
    protein = []
    
    for i in range(0, len(dna_seq) -2,3):
        codon= dna_seq[i:i+3]
        amino_acid =codon_table.get(codon.upper(), 'X')  # 'X' for unknown codons
        if amino_acid =='_':  
            break
        protein.append(amino_acid)
    
    return ''.join(protein)



if __name__ == "__main__":

    dna=read_fasta("example.fasta")
    counts=count_bases(dna)


    gc=gc_content(dna)
    reverse_comp=reverse_complement(dna)
    protein_seq = translate_dna_to_protein(dna)



    print("Dna Sequence:",dna)
    print("n\Base Counts:")
    for base, count in counts.items():


        print(f" {base}: {count}")
    print(f"\nGC Content: {gc:.2f}%")

    print("\n Reverse Complement:")
    print(reverse_comp)


    print("\nProtein Translation:")
    print(protein_seq)
    
  