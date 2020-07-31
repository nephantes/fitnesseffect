#!/usr/bin/env python
# coding: utf-8

import os, sys, getopt
import math, json

def main(argv):
    #get parameters from command line
    script_dir = os.path.dirname(os.path.realpath(__file__))
    outFile = "db.fa"
    sequence = ""
    counttype = "aa"
    opts, args = getopt.getopt(argv,"s:c:o:",["sequence=","counttype=", "outFile="])    
    for opt, arg in opts:
        if opt in ("-s", "--sequence"):
            sequence = arg
        elif opt in ("-c", "--counttype"):
            counttype = arg
        elif opt in ("-o", "--outFile"):
            outFile = arg

    #nucleotides
    dna = ["A","G","C","T"]
    #nt to aminoacid translation table
    table = ""
    with open(script_dir + "/codontable.json") as codontable:
        table = json.load(codontable)

    with open(outFile, "w") as reads:
        seqid = 1
        # countype = "aa" allows quantifying only aminoacid level changes. It can be multiple point mutations in a codon.
        if (counttype == "aa"):
            for j in range(0, int(math.floor(len(sequence)/3))):
                #get three nt from wt sequence
                threent_wt = sequence[j*3:j*3+3]
                amino_wt = table[threent_wt]
                #write whole WT sequence once
                if (j == 0):
                    reads.write(">type=wt:seqid=-:aapos=-:aa:-:codon:-:ntpos=-\n")
                    reads.write(sequence + "\n")
                # For each codon create a sequence with one aminoacid change at a time 
                # Get the mutation type. If mutation is synonymos = syn and stop = stop otherwise mut
                for codon in table:
                    seqlist = list(sequence)

                    if (seqlist[j*3:j*3+3] != codon):
                        seqlist[j*3:j*3+3] = codon
                        repseq = "".join(seqlist)

                        amino_mut = table[codon]
                        type = "mut"
                        if (amino_wt == amino_mut):
                            type="syn"
                        elif (amino_mut == "*"):
                            type="stop"
                        
                        reads.write(">type="+type +":seqid=" + str(seqid) + ":aapos=" + str(j+1) + ":aa="+amino_mut+":codon="+codon+"\n")
                        reads.write(repseq + "\n")
                        seqid +=1

        # countype = "nt" allows quantifying single nucleotide level changes. It CANNOT be multiple point mutations in a codon.
        elif (counttype == "nt"):
            for j in range(0, len(sequence)):
                #get three nt from wt sequence
                index = int(math.floor(j/3))
                threent_wt = sequence[index*3:index*3+3]
                amino_wt = table[threent_wt]
                #write whole WT sequence once
                if (j ==0):
                    reads.write(">type=wt:seqid=-:aapos=-:aa:-:codon:-:ntpos=-\n")
                    reads.write(sequence + "\n")
                for i in range(0, len(dna)):
                    seqlist = list(sequence)
                    if (seqlist[j] != dna[i]):
                        seqlist[j]=dna[i]
                        repseq = "".join(seqlist)

                        threent_mut = repseq[index*3:index*3+3]
                        amino_mut= table[threent_mut]
                        type = "mut"
                        if (amino_wt == amino_mut):
                            type="syn"
                        elif (amino_mut == "*"):
                            type="stop"
                        
                        reads.write(">type="+type +":seqid=" + str( 4*j + i + 1) +":aapos=" + str(index+1) + ":aa="+amino_mut+ ":codon="+threent_mut+ ":ntpos=" + str(j+1) + "\n")
                        reads.write(repseq + "\n")

if __name__ == "__main__":
   main(sys.argv[1:])
