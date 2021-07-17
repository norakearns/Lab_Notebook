def convert_phred(letter):
    QS = ord(letter) - 33
    return QS

def qual_score(phred_score):
    "Takes a string of phred scores, converts them to QS number scores, and then finds the average QS of that string"
    n=0
    for phred in phred_score:
        QS = convert_phred(phred)
        n+=QS
    return((n)/len(phred_score))
    
def validate_base_seq(seq,RNAFlag=False):
    "Takes a string of DNA, makes sure it is entirely A,T,C, and G"
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("U" if RNAFlag else "T") + seq.count("G") + seq.count("C")
    
def gc_content(dna):
    "Returns the % G and C of a DNA string"
    assert validate_base_seq(dna), "String is not DNA"
    DNA = dna.upper()
    gc_content = DNA.count("G") + DNA.count("C")
    return(gc_content/len(dna))

def oneline_fasta(filename):
    "Takes a Fasta file where the sequences are split across multiple lines, and concatenates the sequence lines to a single string"
    file= open(filename, "r")
    all_lines = file.read()
    record_array = all_lines.split('>')
    record_list = []
    for record in record_array[1:]:
        record_line_array = record.split('\n')  
        seq = record_line_array[1:]
        seq_line = ('').join(seq)
        return(seq_line)

DNAbases = "ATGC" 
RNAbases = "ATGC"   
    
def main():  
    assert convert_phred("F") == 37
    print("phred score converted")
    assert qual_score("GHF@EIJ") == 37.42857142857143
    print("quality score is correct")
    assert validate_base_seq("ATGCCCAT") == True
    print("the sequence is all DNA")
    assert oneline_fasta("fasta.test") == 'MSLQMILFFFLLYRVDGQAMLRQKISSTKSQDKTVVIDCDYPSDCYRYIHWYQLKGQTLKRILYAQISGGEPARDAGFELFKIDRKQSNIALKIPELKTEHSAVYYCACWDRSA'
    print("The fasta sequence is all on one line")

if __name__ == "__main__":
    main()

    
