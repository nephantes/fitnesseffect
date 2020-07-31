#!/usr/bin/env python
# coding: utf-8

import sys, getopt, os
import math 

def main(argv):
    #get parameters from command line
    inputcounts = ""
    outFile = "freq.tsv"

    opts, args = getopt.getopt(argv,"i:o:",["inputcounts=", "outFile="])    
    for opt, arg in opts:
        if opt in ("-i", "--inputcounts"):
            inputcounts = arg
        elif opt in ("-o", "--outFile"):
            outFile = arg

    arr = []
    wtsyn = 0
    wtsyn_count = 0
    stop = 0
    stop_count = 0
    wt_count = 0
    new_dic = {}
    #Caluclate A value for mutations and count WTSyn and Stop frequencies
    with open(outFile+"_nt.tsv", "w") as reads:
        with open(inputcounts) as input:
            for line in input:
                arr = line.split()

                if(arr[3] == "wt"):
                    wt_count = int(arr[4])
                elif(arr[3] == "syn"):
                    wtsyn += int(arr[4])
                    wtsyn_count += 1
                elif(arr[3] == "stop"):
                    stop += int(arr[4])
                    stop_count += 1
                if (arr[1] != "-"):
                    #print(str(arr[0]) + "\t" + str(arr[1]) + "\t" +  str(arr[2]) +  "\t" + str(arr[3]) +  "\t" +str(arr[4]))
                    if (arr[1] not in new_dic):
                        new_dic[arr[1]]={}
                    if (arr[2] not in new_dic[arr[1]]):
                        new_dic[arr[1]][arr[2]] = int(arr[4])
                    else:
                        new_dic[arr[1]][arr[2]] += int(arr[4])

                    A = round(math.log2(int(arr[4])/wt_count), 4)
                    reads.write(str(arr[0]) + "\t" + str(arr[1]) + "\t" +  str(arr[2]) + "\t" +str(arr[4]) + "\t" + str(A))
                    if (len(arr) == 6):
                        reads.write("\t" + str(arr[5]))
                    reads.write("\n")

    with open(outFile+"_aa.tsv", "w") as reads:
        for aa in new_dic:
            for pos in new_dic[aa]:
                A = 0
                if(new_dic[aa][pos] != 0 and wt_count != 0):
                    A = round(math.log2(int(new_dic[aa][pos])/wt_count), 4)
                reads.write(aa + "\t" + str(pos) + "\t" +  str(new_dic[aa][pos]) + "\t" + str(A)+ "\n")

    reads.close()
    #Calculate A for WTsyn and stop with their average counts normalized with wt_count that they will be used calculating S
    with open(outFile+".txt", "w") as reads:
        reads.write("wtsyn="+str(round(math.log2(wtsyn/(wt_count*wtsyn_count)), 4)) + "\nstop="+str(round(math.log2(stop/(wt_count*stop_count)), 4)) + "\n")
    reads.close()

    #ADD WILD TYPE AA calculations to the matrix

if __name__ == "__main__":
   main(sys.argv[1:])
