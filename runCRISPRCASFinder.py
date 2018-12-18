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
import pandas as pd
import matplotlib.pyplot as plt
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
                    if i == 'A' or i == 'a':
                        if i.islower():
                            revnucs += 't'
                        else:
                            revnucs += 'T'
                    elif i == 'T' or i == 't':
                        if i.islower():
                            revnucs += 'a'
                        else:
                            revnucs += 'A'
                    elif i == 'G' or i == 'g':
                        if i.islower():
                            revnucs += 'c'
                        else:
                            revnucs += 'C'
                    elif i == 'C' or i == 'c':
                        if i.islower():
                            revnucs += 'g'
                        else:
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
    newfilename = 'fwd_and_rev_' + buildname + '.fasta'
    filepathlist.append(newfilename)
    newpath = "/".join(filepathlist)
    with open(file, 'r') as f:
        with open(newpath, 'w') as w:
            for line in f:
                if line.startswith(">"):
                    line = line.strip('>')
                    line = ">" + buildname + '_' + line
                    #apparently, CCF can\'t handle parentheses
                    line = re.sub('[().]',"",line)
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
        os.system("perl /usr/bin/CRISPRCasFinder.pl -i " + str(input) + " -q -so /opt/vmatch-2.3.0/SELECT/sel392.so -outdir " + str(output) + " -minDR " + str(min) + " -maxDR " + str(max))

        print("Done with CCF")
    except Exception as e:
        print("nope")
        print(str(e))
        sys.exit(2)

def getSpacers(GFFfiles, hits, buildname):
    '''
    Each strain has it's own folder with the crisprcasfinderresults
    Format output from CRISPRCASfinder to suitable fasta files for a blast query.
    Each queryfile contains a fasta sequence with ID for each contig and spacer present.
    Each unique queryfile is a unique organism.
    '''
    hitlist = []
    for h in hits:
        hitname = buildname + GFFfiles + str(h) + '.gff'
        hitlist.append(hitname)

    prog1 = re.compile("CRISPRspacer")
    prog2 = re.compile("sequence=\w+;")

    queryfile = "../spacers/" + buildname + ".fasta"

    with open(queryfile,'w') as w:
        for file in hitlist:
            fn = file.split("/")[-1].strip(".gff")
            spacerID = 0
            with open(file,'r') as f:
                for line in f:
                    m = prog1.search(line)
                    if m is not None:
                        n = prog2.search(line)
                        seq = n.group(0).split("=")[1].strip(";")
                        fastaseq = ">" + fn + "_spacerID_" + str(spacerID) + "\n" + seq + "\n\n"
                        w.write(fastaseq)
                        spacerID += 1
            f.close()
    w.close()
    return queryfile

def getCCFhits(summaryfile):
    # nog toe te voegen aan de main. extractie van hits in CCF
    summary = pd.read_csv(summaryfile,delimiter = '\t')
    hits = summary[summary['Nb CRISPRs'] > 0]
    hits = hits['Sequence(s)']
    return hits

def makeDirs(outputlocation):
    if not os.path.exists(outputlocation):
        os.makedirs(str(outputlocation + "/blastout"))
        os.makedirs(str(outputlocation + "/spacers"))
    else:
        print("Directories already exist")

def getInputfiles(curwd, input):

    fastalist = []
    if os.path.isdir(input):
        #for correct abspath, wd has to be set where the input files are
        #os.chdir(input)
        for file in os.listdir(input):
            #print(os.path.join(input,file))
            if file.endswith(".fasta"):
                fastalist.append(os.path.join(input,file))
    else:
        #just get the single file that is given
        if not input.endswith('.fasta'):
            sys.exit("Not a valid fasta file. Check input argument.")
        else:
            fastalist.append(input)
    return fastalist

def main(args):
    '''
    run CRISPRCasFinder
    run blastn script
    '''
    #check if output path exists. If not, create it
    makeDirs(args.output)

    #runs batch of contigs/scaffolds if dir is given
    curcwd = os.getcwd()
    fastalist = getInputfiles(curcwd,args.input)

    if len(fastalist) == 0:
        sys.exit("Given directory is empty or does not exist. Check input argument.")

    #loop to run single or multiple files. do this in parallel (maybe) later
    absBVDB = os.path.abspath(args.blastviraldb)
    absOutput = os.path.abspath(args.output)
    outpath = absOutput + '/CCF/'
    os.makedirs(outpath)

    for fastafile in fastalist:
        os.chdir(curcwd)
        newpath, buildname = addBuildName(fastafile)
        absInput = os.path.abspath(newpath)
        if args.reverse:
            try:
                print("Fetching reverse complementary fastas ...")
                getReverseComplement(absInput)
            except:
                print('Oops.')

        print("Start CCF ...")
        os.chdir(outpath)

        runCrisprCasFinder(absInput,buildname,args.minimum,args.maximum)
        CCFhits = getCCFhits(str(buildname + '/TSV/CRISPR-Cas_summary.tsv'))
        spacerfasta = getSpacers('/GFF/', CCFhits, buildname)
        #create filepath for blast output
        absspacer = os.path.abspath(os.path.join('.',spacerfasta))
        blastout = absOutput  + "/blastout/" + buildname + "blastout.out"
        print("Run blast ...")
        runLocalBlast.runBlast(absspacer,absBVDB,blastout,args.percidentity)
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
