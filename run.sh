#wget https://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
#bash ~/miniconda.sh -b -p $HOME/miniconda
#export PATH="$HOME/miniconda/bin:$PATH"
#install bowtie 
#conda install -c bioconda bowtie 
#conda install python.app

seq=TTTAATCAACCCTCAGGAGGGGACCCAGAA
#create db_aa.fa and db_nt.fa files for given nt sequence
#countype = "aa" allows quantifying only aminoacid level changes. It can be multiple point mutations in a codon.
mkdir -p test/index
python scripts/generateFasta.py -s $seq -c aa -o test/index/db.fa
bowtie-build test/index/db.fa test/index/db

#generate some random sequences to test the counts
mkdir -p test/input
python scripts/generateRandomReads.py -n 10000 -d test/index/db.fa -o test/input/test10K_p0.fq
python scripts/generateRandomReads.py -n 10000 -d test/index/db.fa -o test/input/test10K_p1.fq

#align the reads
mkdir -p test/bed
bowtie test/index/db -q test/input/test10K_p0.fq -v 0 -n 0 --best test/bed/test10K_p0.bed
bowtie test/index/db -q test/input/test10K_p1.fq -v 0 -n 0 --best test/bed/test10K_p1.bed

#count the reads per sequence
mkdir -p test/counts
awk '{gsub(":", "\t", $3); gsub("=", "\t", $3); print $3}' test/bed/test10K_p0.bed | sort | uniq -c | sort -k7n -k11b | awk '{print $11"\t"$9"\t"$7"\t"$3"\t"$1}' > test/counts/counts_p0.tsv
awk '{gsub(":", "\t", $3); gsub("=", "\t", $3); print $3}' test/bed/test10K_p1.bed | sort | uniq -c | sort -k7n -k11b | awk '{print $11"\t"$9"\t"$7"\t"$3"\t"$1}' > test/counts/counts_p1.tsv

#count frequencies 
mkdir -p test/freqs
python scripts/calcA.py -i test/counts/counts_p0.tsv -o test/freqs/freqs_p0
python scripts/calcA.py -i test/counts/counts_p1.tsv -o test/freqs/freqs_p1

#calc FS scores
mkdir -p test/scores/
python scripts/calcFS.py -p test/freqs/freqs_p0 -r test/freqs/freqs_p1 -o test/scores/score

#draw heatmap
mkdir -p test/plots/
pythonw scripts/drawPlots.py -i test/scores/score_F_aa.tsv -o test/plots -s $seq -t aa -f F
pythonw scripts/drawPlots.py -i test/scores/score_S_aa.tsv -o test/plots -s $seq -t aa -f S 
pythonw scripts/drawPlots.py -i test/scores/score_F_nt.tsv -o test/plots -s $seq -t nt -f F
pythonw scripts/drawPlots.py -i test/scores/score_S_nt.tsv -o test/plots -s $seq -t nt -f S


