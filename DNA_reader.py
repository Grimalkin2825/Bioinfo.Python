from Bio.Seq import Seq

# ===============================
# USER SEQUENCE VALIDATION
# ===============================
while True:
    adn_input = input("Entre une séquence ADN (A,T,C,G) : ").upper()

    sequence_valide = True
    for base in adn_input:
        if base not in "ATCG":
            print("❌ Séquence invalide ! Veuillez recommencer.")
            sequence_valide = False
            break

    if sequence_valide and adn_input != "":
        break

print("✅ Séquence validée :", adn_input)

# Objet BioPython (pour les fonctions BioPython)
adn = Seq(adn_input)

print("\n Séquence ADN :", adn)
print(" Longueur :", len(adn))
print(" Complément :", adn.complement())
print(" Reverse complément :", adn.reverse_complement())
print(" ARN :", adn.transcribe())
print(" Protéine (frame +1):", adn.translate())

# ===============================
# SIX FRAME TRANSLATION
# ===============================

# Table codons ADN → AA (1 lettre)
codon_table = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W',

'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'CGT':'R','CGC':'R','CGA':'R','CGG':'R',

'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R',

'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# Reverse complement (version string)
def reverse_complement(seq):
    complement = {"A":"T","T":"A","C":"G","G":"C"}
    rev_comp = "".join(complement[base] for base in seq)
    return rev_comp[::-1]

# Traduction d’une frame
def translate_frame(seq, frame):
    protein = ""
    for i in range(frame, len(seq)-2, 3):
        codon = seq[i:i+3]
        protein += codon_table.get(codon, "?")
    return protein

# Traduction des 6 frames
def translate_six_frames(dna):
    dna = dna.upper()
    rev = reverse_complement(dna)

    print("\n==============================")
    print(" SIX FRAME TRANSLATION")
    print("==============================")

    print("\n--- BRIN + ---")
    for frame in range(3):
        prot = translate_frame(dna, frame)
        print(f"Frame +{frame+1}: {prot}")

    print("\n--- BRIN - (reverse complement) ---")
    for frame in range(3):
        prot = translate_frame(rev, frame)
        print(f"Frame -{frame+1}: {prot}")

# ▶️ LANCEMENT
translate_six_frames(adn_input)
