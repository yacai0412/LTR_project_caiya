import sys
import gzip

def main():
    snv = sys.argv[1] # AA.deepconsensus.cx3.fwd.retrogene.split.gene.fwd.AAC.mm2_splice.F904.shared.hap1.fwd_plus_minus.subpass_5.vaf_0.8.snv
    split_fastq = sys.argv[2]

    out1 = snv + ".readpos.out1"

    fastq_len_dic, fastq_count = read_fastq_get_rlen(split_fastq)
    readpos_ls = read_insnv(snv, fastq_len_dic)
    get_realtive_pos1_ref(readpos_ls, fastq_count, out1)

def read_fastq_get_rlen(fastq):
    fastq_len_dic = {}
    fastq_count = 0
    c = 0
    
    # Check if the file is gzipped
    if fastq.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'  # Read as text
    else:
        open_func = open
        mode = 'r'   # Read as regular text file
    
    with open_func(fastq, mode) as inf:
        for line in inf:
            c += 1
            if c % 4 == 1:
                seqname = line.rstrip().lstrip("@")
                fastq_count += 1
            elif c % 4 == 2:
                rlen = len(line.rstrip())
                fastq_len_dic[seqname] = rlen
                
    return fastq_len_dic, fastq_count


def read_insnv(snv, fastq_len_dic):
    readpos_ls = []
    with open(snv, "r") as inf:
        for line in inf:
            line = line.rstrip("\n")
            ll = line.split("\t")

            readname = ll[6].split(":")[0]
            readpos = int(ll[6].split(":")[1])

            if readname in fastq_len_dic:
                readlength = fastq_len_dic[readname]
                if readpos <= readlength / 2:
                    readpos1 = readpos
                else:
                    readpos1 = -1 * (readlength - readpos + 1)
                readpos_ls.append(readpos1)
            else:
                print("error in " + readname)
            
    return readpos_ls


def get_realtive_pos1_ref(relative_pos1_ls, fastq_count, out1):
    relative_pos1_ref_dic = {}
    max_relativepos1 = max(relative_pos1_ls)
    min_relativepos1 = min(relative_pos1_ls)

    max_abs_relativepos1 = max(abs(max_relativepos1), abs(min_relativepos1))
    for i in range(max_abs_relativepos1):
        relative_pos1_ref_dic[i + 1] = 0
        relative_pos1_ref_dic[-1 * (i + 1)] = 0

    for pos1 in relative_pos1_ls:
        relative_pos1_ref_dic[pos1] += 1
    
    with open(out1, "w") as final1:
        for pos11 in relative_pos1_ref_dic:
            if pos11 == 0:
                outtag1 = "boundary"
            elif pos11 < 0:
                outtag1 = "right"
            elif pos11 > 0:
                outtag1 = "left"

            mutation_rate = relative_pos1_ref_dic[pos11] / fastq_count
            out_line = str(pos11) + "\t" + str(relative_pos1_ref_dic[pos11]) + "\t" + str(fastq_count) + "\t" + str(mutation_rate) + "\t" + outtag1
            final1.write(out_line + "\n")

if __name__ == "__main__":
    main()

