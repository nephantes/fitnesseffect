#!/usr/bin/env python
# coding: utf-8

import sys, getopt, os, os.path
import math 

def readFreqToDict(type, filename):
    new_dic = {}
    indexOffset = 1
    if (type == "aa"):
        indexOffset = 0
    with open(filename) as input:
        for line in input:
            arr = line.rstrip().split()
            if (arr[0] not in new_dic):
                new_dic[arr[0]]={}
            if (arr[1+indexOffset] not in new_dic[arr[0]]):
                new_dic[arr[0]][arr[1+indexOffset]] = arr[3+indexOffset]
    return new_dic

def readNormFactToDict(filename):
    new_dic = {}
    with open(filename) as input:
        for line in input:
            arr = line.rstrip().split("=")
            if (arr[0] not in new_dic):
                new_dic[arr[0]]=arr[1]
    return new_dic

def writeTable(type, FS, dic, out):
    with open(out + "_" + FS + "_" + type + ".tsv", "w") as output:
        output.write(type + "\tpos\t" + FS + "\n")
        for key in dic:
            for pos in dic[key]:
                output.write(key + "\t" + str(pos) + "\t" +  str(dic[key][pos]) + "\n")

def calcFS(type, p0, p1, out):
    #Open tsv file for p0 freqs
    freqP0 = readFreqToDict(type, p0+"_"+type+".tsv")
    #Open tsv file for p1 freqs
    freqP1 = readFreqToDict(type, p1+"_"+type+".tsv")
    #Open txt file for normalization factors for p0
    normFactp0 = readNormFactToDict(p0+".txt")
    #Open txt file for normalization factors for p1
    normFactp1 = readNormFactToDict(p1+".txt")

    F_dic = {}
    for key in freqP0.keys():
        if (key not in F_dic):
            F_dic[key] = {}
        for pos in freqP0[key].keys():
            F_dic[key][pos] = float(freqP1[key][pos]) - float(freqP0[key][pos])
    
    Fwtsyn = float(normFactp1["wtsyn"]) - float(normFactp0["wtsyn"])
    Fstop = float(normFactp1["stop"]) - float(normFactp0["stop"])

    writeTable(type, "F", F_dic, out)

    S_dic = {}
    for key in F_dic.keys():
        if (key not in S_dic):
            S_dic[key] = {}
        for pos in F_dic[key].keys():
            S_dic[key][pos] = (F_dic[key][pos] - Fwtsyn) / (Fwtsyn -  Fstop)
    
    writeTable(type, "S", S_dic, out)




def main(argv):
    #get parameters from command line
    inputP0 = ""
    inputP1 = ""
    outFile = "freq.tsv"

    opts, args = getopt.getopt(argv,"p:r:o:",["inputP0=", "inputP1=", "outFile="])    
    for opt, arg in opts:
        if opt in ("-p", "--inputP0"):
            inputP0 = arg
        elif opt in ("-r", "--inputP1"):
            inputP1 = arg
        elif opt in ("-o", "--outFile"):
            outFile = arg

    # calcF for aa
    calcFS("aa", inputP0, inputP1, outFile)
    # calcF for nt
    calcFS("nt", inputP0, inputP1, outFile)

if __name__ == "__main__":
   main(sys.argv[1:])
