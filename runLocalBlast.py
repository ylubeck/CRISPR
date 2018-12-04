##runlocalblast
import sys
import argparse
import os
#from Bio.Blast.Applications import NcbiblastnCommandline

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
    # run contigs vs viral database.
    #run blastn
    print(" ".join([query,blastvirdb,output,str(pid)]))
    #niet de goede outputformat!
    #blastn = NcbiblastnCommandline(query = query, db = blastvirdb, out = output, perc_identity = pid, outfmt = "6")
    #stdout, stderr = blastn()

    os.system("blastn -query " + query + " -db " + blastvirdb + " -outfmt \"6 std pident nident qseq sseq\" -out " + output)
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
