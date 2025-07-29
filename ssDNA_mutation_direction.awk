awk 'BEGIN {
    FS="\t";
    complement["A"] = "T";
    complement["T"] = "A";
    complement["C"] = "G";
    complement["G"] = "C";
}

{
    # Convert ref and alt bases to uppercase
    ref = toupper($3);
    alt = toupper($4);
    
	# get strand
	split($7, reads, ",");
	split(reads[1], sp, ":");
	strand=sp[3];
	
	if(strand == "+"){
		mutation_type = ref "->" alt
	}
	else if(strand == "-"){
		mutation_type = complement[ref] "->" complement[alt]
	};
	counts[mutation_type]++;    
}

END {
    # Print the summary of mutation types
    for (m in counts) {
        print m "\t" counts[m];
    }
}' 

