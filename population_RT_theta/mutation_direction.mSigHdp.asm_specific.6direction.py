import sys

def main():
    snv = sys.argv[1] # SV.0328.CE2.INS.A1.alt.consensus.asm.s.bam.read_specific_relativepos_uniq.inside_indel_specific.del0.0bp_length_cutoff.indel_1bp.relativepos
    # ref = sys.argv[2] # /rd/caiya/LTR/population_RT_theta/SV.0328.CE2.INS.$i.alt.consensus.fasta
    sample_name = sys.argv[2] # 

    out = ".".join(snv.split("/")[-1].split(".")[:-1]) + ".csv"

    # ref_fasta_dic = read_ref_fasta(ref)
    # output_mutation_dict(mutation_dict, sample_name, out)

    mutation_dict = mutation_direction_6direction(snv)
    output_mutation_6direction_dict(mutation_dict, sample_name, out)


def initialize_mutation_dictionary():
    nucleotides = ["A", "C", "G", "T"]
    mutation_dict = {
        "C>A": {f"{a}C{b}": 0 for a in nucleotides for b in nucleotides},
        "C>G": {f"{a}C{b}": 0 for a in nucleotides for b in nucleotides},
        "C>T": {f"{a}C{b}": 0 for a in nucleotides for b in nucleotides},
        "T>A": {f"{a}T{b}": 0 for a in nucleotides for b in nucleotides},
        "T>C": {f"{a}T{b}": 0 for a in nucleotides for b in nucleotides},
        "T>G": {f"{a}T{b}": 0 for a in nucleotides for b in nucleotides}
    }
    return mutation_dict

def read_ref_fasta(ref):
    ref_fasta_dic = {}
    with open(ref, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line.lstrip(">").split(":")[0]
                ref_fasta_dic[name] = ""
            else:
                ref_fasta_dic[name] = ref_fasta_dic[name] + line
    return ref_fasta_dic

def get_reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement = ''.join(complement[nuc] for nuc in reversed(dna_sequence))
    return reverse_complement

def read_snv(snv, ref_fasta_dic):
    mutation_dict = initialize_mutation_dictionary()
    with open(snv, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            chrom = ll[0]
            spos = int(ll[3])
            ref = ll[4].upper()
            alt = ll[5].upper()

            if "N" not in ref_fasta_dic[chrom][spos - 2: spos + 1].upper():
                if ref == "C" or ref == "T":
                    Mutation_type = ref + ">" + alt
                    Trinucleotide = ref_fasta_dic[chrom][spos - 2: spos + 1].upper()
                elif ref == "G" or ref == "A":
                    Mutation_type = get_reverse_complement(ref) + ">" + get_reverse_complement(alt)
                    Trinucleotide = get_reverse_complement(ref_fasta_dic[chrom][spos - 2: spos + 1].upper())

                if Mutation_type in mutation_dict and Trinucleotide in mutation_dict[Mutation_type]:
                    mutation_dict[Mutation_type][Trinucleotide] += 1
                else:
                    print("error in " + line)
                    print(Mutation_type + "\t" + Trinucleotide)
        
    return mutation_dict

def mutation_direction_6direction(snv):
    mutation_dict = {}
    with open(snv, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            ref = ll[4].upper()
            alt = ll[5].upper()

            if ref == "C" or ref == "T":
                Mutation_type = ref + ">" + alt + "/" + get_reverse_complement(ref) + ">" + get_reverse_complement(alt)
            elif ref == "G" or ref == "A":
                Mutation_type = get_reverse_complement(ref) + ">" + get_reverse_complement(alt) + "/" + ref + ">" + alt
            
            if Mutation_type in mutation_dict:
                mutation_dict[Mutation_type] += 1
            else:
                mutation_dict[Mutation_type] = 1
    return mutation_dict

def output_mutation_6direction_dict(mutation_dict, sample_name, out):
    with open(out, "w") as final:
        final.write("Mutation type," + sample_name + "\n")
    
        for Mutation_type in sorted(mutation_dict):
            mutaion_count = mutation_dict[Mutation_type]
            out_line = Mutation_type + "," + str(mutaion_count)
            final.write(out_line + "\n")




def output_mutation_dict(mutation_dict, sample_name, out):
    with open(out, "w") as final:
        final.write("Mutation type,Trinucleotide," + sample_name + "\n")
    
        for Mutation_type in sorted(mutation_dict):
            for Trinucleotide in sorted(mutation_dict[Mutation_type]):
                mutaion_count = mutation_dict[Mutation_type][Trinucleotide]
                out_line = Mutation_type + "," + Trinucleotide + "," + str(mutaion_count)
                final.write(out_line + "\n")

if __name__ == "__main__":
    main()

