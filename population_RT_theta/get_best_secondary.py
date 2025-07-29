import sys
import pysam

def main():
    inbam = sys.argv[1]

    outbam = ".".join(inbam.split("/")[-1].split(".")[:-1]) + ".best_sec.bam"

    best_secondary_alignment(inbam, outbam)


def best_secondary_alignment(input_bam, output_bam):
    samfile = pysam.AlignmentFile(input_bam, "rb")
    filtered = pysam.AlignmentFile(output_bam, "wb", template=samfile)

    best_alignments = {}

    for read in samfile.fetch():
        if read.is_secondary:
            aln_score = read.get_tag('AS')

            if read.query_name not in best_alignments or aln_score > best_alignments[read.query_name][1]:
                best_alignments[read.query_name] = (read, aln_score)
        else:
            filtered.write(read)

    for read, _ in best_alignments.values():
        filtered.write(read)

    samfile.close()
    filtered.close()


if __name__ == "__main__":
    main()
