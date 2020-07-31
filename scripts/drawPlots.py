#!/usr/bin/env python
# coding: utf-8

import sys, getopt, csv, os, json
import pandas as pd
import math
import numpy as np
#plot the heatmap
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns

def main(argv):
    #get parameters from command line
    script_dir = os.path.dirname(os.path.realpath(__file__))

    inFile = ""
    out = ""
    sequence = ""
    type = ""
    FS = ""

    opts, args = getopt.getopt(argv,"i:o:s:t:f:",["inFile=","out=","sequence=","type=", "FS="])  
    for opt, arg in opts:
        if opt in ("-i", "--inFile"):
            inFile = arg
        elif opt in ("-o", "--outFile"):
            out = arg
        elif opt in ("-s", "--sequence"):
            sequence = arg
        elif opt in ("-t", "--type"):
            type = arg
        elif opt in ("-f", "--FS"):
            FS = arg
    
    #nt to aminoacid translation table
    table = ""
    with open(script_dir + "/codontable.json") as codontable:
        table = json.load(codontable)

    #read in wild-type nt sequence
    length =  int(math.floor(len(sequence)/3))
    
    wt_seq = []
    for j in range(0, int(math.floor(len(sequence)/3))):
        #get three nt from wt sequence
        threent_wt = sequence[j*3:j*3+3]
        if (type == "aa"):
            wt_seq.append(table[threent_wt])
        else:
            wt_seq.append(threent_wt)

    #read tsv file into a data frame
    df = pd.read_csv(inFile, delimiter = '\t')

    #create a new data frame with the fitness effect for every amino acid at each position
    mat = pd.pivot_table(df, index = type, columns = 'pos', values = FS, aggfunc = np.median)

    s = pd.Series(df[FS])
    ax = s.plot.hist()
    ax.figure.savefig(out + "/hist_" + FS + "_"+ type + ".pdf")
    plt.close()

    # Plot heatmap
    sns.set()
    ax = sns.heatmap(mat, square = 'FALSE', xticklabels='auto', yticklabels='auto')
    # prepare wt sequence boxes
    for i in range(0,length):
        ax.add_patch(Rectangle((i, list(mat.index).index(wt_seq[i])), 1, 1, fill=False, edgecolor='blue', lw=3))
    # set the title
    ax.set_title(FS)
    plt.savefig(out + "/heatmap_" + FS + "_"+ type + ".pdf", bbox_inches='tight')



if __name__ == "__main__":
   main(sys.argv[1:])
