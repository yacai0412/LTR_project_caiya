import sys
import math

def main():
    qv = int(sys.argv[1]) # 20 30 ...
    infastq = sys.argv[2:]
    out = ".".join(infastq[0].split("/")[-1].split(".")[:-2]) + ".q" + str(qv) + ".fastq"

    fastq_dic, rq_dic = read_fastq_get_rq(infastq)    
    filter_qv_fwd_rev(fastq_dic, rq_dic, qv, out)

def read_fastq_get_rq(infastq):
    fastq_dic = {}
    rq_dic = {}

    for fq in infastq:
        with open(fq, "r") as inf:
            while True:
                header = inf.readline().strip().lstrip("@")
                if not header:
                    break  # End of file
                zmw = "/".join(header.split("/")[:2])
                strand = header.split("/")[3]
                sequence = inf.readline().strip()
                plus = inf.readline().strip()
                quality = inf.readline().strip()
                fastq_dic.setdefault(zmw, {})[strand] = sequence + "\n" + plus + "\n" + quality

                error_probabilities = [10**((ord(q) - 33)/-10.0) for q in quality]
                average_error_probability = sum(error_probabilities) / len(error_probabilities)
                qv = -10 * math.log10(average_error_probability)
                rq_dic.setdefault(zmw, {})[strand] = qv

    return fastq_dic, rq_dic


def filter_qv_fwd_rev(fastq_dic, rq_dic, qv, out):
    with open(out, "w") as final:
        for zmw in rq_dic:
            if len(rq_dic[zmw]) == 2 and "fwd" in rq_dic[zmw] and "rev" in rq_dic[zmw] and rq_dic[zmw]["fwd"] >= qv and rq_dic[zmw]["rev"] >= qv:
                fwd_out_line = "@" + zmw + "/deepconsensus/fwd" + "\n" + fastq_dic[zmw]["fwd"]
                rev_out_line = "@" + zmw + "/deepconsensus/rev" + "\n" + fastq_dic[zmw]["rev"]
                final.write(fwd_out_line + "\n")
                final.write(rev_out_line + "\n")



if __name__ == "__main__":
    main()
