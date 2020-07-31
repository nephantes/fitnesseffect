#!/usr/bin/env python
# coding: utf-8

import sys, getopt
import random

def main(argv):
    #get parameters from command line
    numReads = 100000
    outFile = "reads.fq"

    opts, args = getopt.getopt(argv,"n:d:o:",["numReads=", "dbfasta=", "outFile="])    
    for opt, arg in opts:
        if opt in ("-n", "--numReads"):
            numReads = arg
        elif opt in ("-d", "--dbfasta"):
            dbfasta = arg
        elif opt in ("-o", "--outFile"):
            outFile = arg

    dbarray = []
    with open(dbfasta) as input:
        for line1 in input:
            line2 = next(input)
            dbarray.append(line1+";"+line2)

    with open(outFile, "w") as reads:
        for i in range(0, int(numReads)):
            randseq = dbarray[random.randint(0, len(dbarray)-1)]
            seqarr = randseq.split(";")
            reads.write(seqarr[0].replace(">", "@").rstrip() + ":read#="+str(i+1) + "\n")
            reads.write(seqarr[1])
            reads.write("+\n||||||||||||||||||||||||||||||\n")

if __name__ == "__main__":
   main(sys.argv[1:])
