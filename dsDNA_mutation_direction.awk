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
    
    # Get the reverse mutation
    reverse = complement[ref] "->" complement[alt];

    # Combine mutation and its reverse in a sorted manner to avoid duplicates
    if (ref alt < complement[ref] complement[alt]) {
        mutation_type = ref "->" alt "/" reverse;
    } else {
        mutation_type = reverse "/" ref "->" alt;
    }

    # Count each combined mutation type
    counts[mutation_type]++;
}

END {
    # Print the summary of mutation types
    for (m in counts) {
        print m "\t" counts[m];
    }
}' 

