##runlocalblast
import sys
import argparse
import os
import matplotlib.pyplot as plt

#from Bio.Blast.Applications import NcbiblastnCommandline
def seqPlot(lengthlist,min,max):
    plt.hist(lengthlist,bins = range(min,max +1,1))
    plt.show()

def seqlength(seq):
    with open(seq,'r') as s:
        seqLengths = []
        for line in s:
            if line.startswith(">") or line in ['\n']:
                continue
            else:
                seqLengths.append(len(line))
    s.close()
    return seqLengths


def makeShortlist(blastout,pid):
    #extract all sequences with %identity > 95%
    f = open(blastout,'r')
    lines = f.readlines()

    f.close()
    f = open(blastout,'w')

    for line in lines:
        entries = line.split("\t")
        if float(entries[2]) > pid:
            f.write(line)
    f.close()

def makeBlastDB(virdbfastas,blastvirdb):
    os.system("makeblastdb -in " + virdbfastas + " -dbtype nucl -out " + blastvirdb + " -parse_seqids")

def runBlast(query,blastvirdb,output,pid):
    #run blastn
    print(" ".join([query,blastvirdb,output,str(pid)]))
    os.system("blastn -query " + query + " -word_size 7 -db " + blastvirdb + " -evalue 10 -outfmt \"6 std pident nident\" -out " + output)
    #os.system("blastn -query " + query + " -db " + blastvirdb + " -evalue 100 -out " + output)
    makeShortlist(output,pid)

###for testing purposes. implement in main program later

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', help = 'input query for blast', required = True)
        parser.add_argument('-db', '--database', help = 'input viral database for blast', default = 'nr')
        parser.add_argument('-dbout', '--databaseout', help = 'output viral database location for blast', default = 'nr')
        parser.add_argument('-o', '--output', help = 'output files location fasta', required = True)
        parser.add_argument('-p', '--percidentity', help = 'percentage sequence identity BLAST', default = 95)

        args = parser.parse_args()
        #print(args)
        if args.db != 'nr':
            makeBlastDB(args.database,args.databaseout)
        runBlast(args.input,args.databaseout,args.output,args.percidentity)
    except:
        print("See -h for help")
        sys.exit(1)
