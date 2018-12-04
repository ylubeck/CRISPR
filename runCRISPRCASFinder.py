'''
Main program:
- format inputfiles: must handle batch runs of contigs
- format outputfiles: Must be handled by BLAST
- arguments passed to CRISPRCASfinder
- arguments passed to blast
- split sequences
TODO:
- CRISPRTOOL install: wait for install of dependencies
- test .gff fetcher -> output to blast (new script file?)
'''

import sys
import argparse
import re
import runLocalBlast
import os
import errno
from os import listdir
from os.path import isfile, join, splitext

def getContigs(seqnames):
    contiglist = list()
    prog = re.compile("contig_\d+")

    for n in seqnames:
        m = prog.search(n)
        result = m.group(0)
        contiglist.append(result)
    return contiglist

def getSequences(mypath):
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    seqnamelist = list()

    for filename in files:
        filelocation = '/'.join([mypath, filename])
        with open(filelocation, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    seqnamelist.append(line)
                    #print(line)
    return seqnamelist


def getReverseComplement(contigs):
    '''
    get the reverse complement of a sequence
    writes to the same folder where the contigs are and adds suffix 'reversed'
    '''

    lines = ''
    with open(contigs, 'r') as r:
        lines = r.readlines()
    r.close()

    reversedsequence = ''
    with open(contigs,'a') as a:
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if reversedsequence != '':
                    a.write(reversedsequence)
                a.write("\n")
                a.write(line + "_reversed\n")
                reversedsequence = ''
            else:
                nucs = list(line)
                revnucs = ''
                for i in nucs[::-1]:
                    if i == 'A':
                        revnucs += 'T'
                    elif i == 'T':
                        revnucs += 'A'
                    elif i == 'G':
                        revnucs += 'C'
                    elif i == 'C':
                        revnucs += 'G'
                    else:
                        revnucs += 'N'
                reversedsequence = revnucs + "\n" + reversedsequence
        a.write(reversedsequence)
    a.close()

def addBuildName(file):
    # add the filename in front of each fasta sequence of the input.
    # copies inputfile and returns new file and path
    buildname = file.strip('.fasta').split('/')[-1] #get filename without extension
    filepathlist = file.split('/')[:-2] #get path

    newfilename = 'named_' + buildname + '.fasta'
    filepathlist.append(newfilename)
    newpath = "/".join(filepathlist)
    with open(file, 'r') as f:
        with open(newpath, 'w') as w:
            for line in f:
                if line.startswith(">"):
                    line = line.strip('>')
                    line = ">" + buildname + '_' + line
                    w.write(line)
                else:
                    w.write(line)
        w.close()
    f.close()
    return newpath,buildname

def runCrisprCasFinder(input, output, min, max):
    #more options can be added later.
    #TODO: wait until CCF dependencies are installed
    try:
        os.system("perl /usr/bin/CRISPRCasFinder.pl -i " + str(input) + " -outdir " + str(output) + " -minDR " + str(min) + " -maxDR " + str(max))
        #os.system("CRISPRCasFinder.pl -help")
        print("Done with CCF")
    except Exception as e:
        print("fuck")
        print(str(e))
        sys.exit(2)


def getSpacers(GFFfiles):
    '''
    Each strain has it's own folder with the crisprcasfinderresults
    Format output from CRISPRCASfinder to suitable fasta files for a blast query.
    Each queryfile contains a fasta sequence with ID for each contig and spacer present.
    Each unique queryfile is a unique organism.
    '''

    GFFfiles += "/GFF"
    prog1 = re.compile("CRISPRspacer")
    prog2 = re.compile("sequence=\w+;")

    queryfile = GFFfiles.split("/")[-2] + ".fasta"

    with open(queryfile,'w') as w:
        for file in GFFfiles:
            fn = file.split("/")[-1].strip(".gff")
            spacerID = 0
            with open(file,'r') as f:
                for line in f:
                    m = prog1.search(line)
                    result = m.group(0)
                    if result is not None:
                        n = prog2.search(line)
                        seq = n.group(0).split("=").strip(";")
                        fastaseq = ">" + fn + "_spacerID_" + str(spacerID) + "\n" + seq + "\n"
                        w.write(fastaseq)
                        spacerID += 1
            f.close()
    w.close()
    return queryfile

def main(args):
    '''
    run CRISPRCasFinder
    run blastn script
    '''
    #check if output path exists. If not, create it
    try:
        os.makedirs(args.output)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #runs batch of contigs/scaffolds if dir is given
    fastalist = []
    if os.path.isdir(args.input):
        #get all files that end with .fasta
        fastalist = [os.path.abspath(f) for f in os.listdir(args.input) if f.endswith('.fasta')]
        for f in fastalist:
            print(f)
    else:
        #just get the single file that is given
        if not f.endswith('.fasta'):
            sys.exit("Not a valid fasta file. Check input argument.")
        else:
            fastalist.append(args.input)
    if len(fastalist) == 0:
        sys.exit("Given directory is empty or does not exist. Check input argument.")

    #loop to run single or multiple files. do this in parallel (maybe) later
    #run CRISPRCasFinder. Creates output in outdir (if given, otherwise created)
    for fastafile in fastalist:
        newpath, buildname = addBuildName(fastafile)
        absInput = os.path.abspath(newpath)
        absOutput = os.path.abspath(args.output)
        if args.reverse:
            try:
                print("Fetching reverse complementary fastas ...")
                getReverseComplement(absInput)
            except Exception as e:
                print(str(e))

        print("Start CCF ...")
        runCrisprCasFinder(absInput,absOutput,args.minimum,args.maximum)

        #create filepath for blast output
        blastout = absOutput + "/" + buildname + "blastout.out"
        print("Run blast ...")
        runLocalBlast.runBlast(newpath,args.blastviraldb,blastout,args.percidentity)
        print("Done with blast")

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', help = 'input folder of query fasta sequences of bacteria', required = True)
        parser.add_argument('-o', '--output', help = 'output directory location', default = 'output')
        parser.add_argument('-r', '--reverse', help = 'add this argument if you want to generate the reverse complement', default = False, action = 'store_true')
        parser.add_argument('-vg','--viralfastas', help = 'viral fasta genomes for blast db', default ='nr')
        parser.add_argument('-bvdb', '--blastviraldb', help = 'blast database with phage genomes', default = 'nr')
        parser.add_argument('-min','--minimum', help = 'minimal length crispr repeat', default = 23)
        parser.add_argument('-max','--maximum', help = 'maximal length crispr repeat', default = 55)
        parser.add_argument('-pid','--percidentity', help = 'minimum sequence identity blastn', default = 95)

        args = parser.parse_args()

        main(args)
    except:
        print("See -h for help")
        sys.exit(1)
    #main(args)
