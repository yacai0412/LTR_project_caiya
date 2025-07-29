import sys
import pysam

def main():
    minoverlap = 0.9

    inbam = sys.argv[1]
    outbam = ".".join(inbam.split("/")[-1].split(".")[:-1]) + ".overlap" + str(int(minoverlap * 100)) + ".bam"
    
    read_in_bam(inbam, outbam, minoverlap)


def get_tag0(qtag):
    if qtag == 0 or qtag == 16:
        tag0 = "pri"
    elif qtag == 2048 or qtag == 2064:
        tag0 = "sup"
    return tag0

def get_tag1(qtag):
    if qtag == 0 or qtag == 2048:
        tag0 = "+"
    elif qtag == 16 or qtag == 2064:
        tag0 = "-"
    return tag0

def get_mapping_length(CIGAR):
    mapping_length = 0
    current_num = ''
    for char in CIGAR:
        if char.isdigit():
            current_num += char
        else:
            if char == "D" or char == "X" or char == "N" or char == "=":
                mapping_length += int(current_num)
            current_num = ''
    return mapping_length

def get_overlap(qkey0, qkey1, qvalue0, qvalue1, minoverlap):
    qkey0_ls = qkey0.split("|")
    qkey1_ls = qkey1.split("|")
    qvalue0_ls = qvalue0.split("|")
    qvalue1_ls = qvalue1.split("|")

    readlength0 = int(qvalue0_ls[3])
    readlength1 = int(qvalue1_ls[3])

    overlap = 0
    if qkey0_ls[0] != qkey1_ls[0] and qkey0_ls[1] == qkey1_ls[1] and qkey0_ls[2] != qkey1_ls[2]:
        if qvalue0_ls[0] == qvalue1_ls[0] and ((int(qvalue0_ls[1]) <= int(qvalue1_ls[2]) and int(qvalue0_ls[2]) >= int(qvalue1_ls[1])) or (int(qvalue1_ls[1]) <= int(qvalue0_ls[2]) and int(qvalue1_ls[2]) >= int(qvalue0_ls[1]))):
            overlap = min(int(qvalue0_ls[2]), int(qvalue1_ls[2])) - max(int(qvalue0_ls[1]), int(qvalue1_ls[1])) + 1

    overlap_rate0 = overlap / readlength0
    overlap_rate1 = overlap / readlength1
    if overlap_rate0 >= minoverlap and overlap_rate1 >= minoverlap:
        return True
    else:
        return False

def get_readline(qzmw0, qkey0, qkey1, readname_count_seq_dic, final):
    qname0 = qzmw0 + "/" + qkey0.split("|")[0]
    qname1 = qzmw0 + "/" + qkey1.split("|")[0]    
    c1_0 = int(qkey0.split("|")[3])
    c1_1 = int(qkey1.split("|")[3])

    line0 = readname_count_seq_dic[qname0][c1_0]
    line1 = readname_count_seq_dic[qname1][c1_1]
    final.write(line0)
    final.write(line1)
        
def get_find_overlap_read(qzmw, pri_sup_dic, qkey0, qvalue0, minoverlap, readname_count_seq_dic, qzmw0, final):
    if qzmw in pri_sup_dic:
        c = 0
        for qkey1 in pri_sup_dic[qzmw]:
            qvalue1 = pri_sup_dic[qzmw][qkey1]
            if qvalue1 != 0 and get_overlap(qkey0, qkey1, qvalue0, qvalue1, minoverlap):
                get_readline(qzmw0, qkey0, qkey1, readname_count_seq_dic, final)
                pri_sup_dic[qzmw][qkey1] = 0
                c += 1
        if c == 0:
            pri_sup_dic.setdefault(qzmw, {})[qkey0] = qvalue0
    else:
        pri_sup_dic.setdefault(qzmw, {})[qkey0] = qvalue0

def read_in_bam(inbam, outbam, minoverlap):    
    # chr_dic={"chr2L":1, "chr2R":1, "chr3L":1, "chr3R":1, "chr4":1, "chrX":1, "chrY":1}
    inf = pysam.AlignmentFile(inbam, "rb")
    final = pysam.AlignmentFile(outbam, "wb", template = inf)

    # Get the list of reference sequences (chromosomes)
    references = inf.references

    for ref in references:
        # if ref in chr_dic:
        readname_count_seq_dic = {}
        pri_dic = {}
        sup_dic = {}

        for read in inf.fetch(ref):
            qname = read.query_name
            if qname not in readname_count_seq_dic:
                c1 = 0
            else:
                c1 = max(readname_count_seq_dic[qname].keys()) + 1
            readname_count_seq_dic.setdefault(qname, {})[c1] = read

            qzmw = int(qname.split("/")[1])
            qzmw0 = qname.rstrip("/fwd").rstrip("/rev")
            qtag = int(read.flag)
            qstrand = qname.split("/")[3]
            tag0 = get_tag0(qtag) # pri or sup
            tag1 = get_tag1(qtag) # + or -
            
            qcigar = read.cigarstring
            readlength = get_mapping_length(qcigar)

            target_chr = read.reference_name
            start_pos = read.reference_start # 0-base position
            end_pos = read.reference_end # 0-base position

            qkey0 = qstrand + "|" + tag0 + "|" + tag1 + "|" + str(c1) # fwd|pri|+|1
            qvalue0 =  target_chr + "|" + str(start_pos) + "|" + str(end_pos) + "|" + str(readlength)
            if tag0 == "pri":
                pri_sup_dic = pri_dic
            elif tag0 == "sup":
                pri_sup_dic = sup_dic
            get_find_overlap_read(qzmw, pri_sup_dic, qkey0, qvalue0, minoverlap, readname_count_seq_dic, qzmw0, final)
            
    inf.close()
    final.close()

if __name__ == "__main__":
    main()
