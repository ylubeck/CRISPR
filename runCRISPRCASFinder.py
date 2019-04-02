'''
Main program:
- format inputfiles: must handle batch runs of contigs
- format outputfiles: Must be handled by BLAST
- arguments passed to CRISPRCASfinder
- arguments passed to blast
- split sequences
TODO:
- start making database
- wait for install prodigal
'''
import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import re
import runLocalBlast
import os
import json

def getContigNode(seqname):
    #works for spades and clc output (scaffolds/contigs)
    prog1 = re.compile("contig_\d+")
    prog2 = re.compile("NODE_\d+")
    m1 = prog1.search(seqname)
    m2 = prog2.search(seqname)

    if m1 is not None:
        s = m1
    elif m2 is not None:
        s = m2
    else:
        prog = re.compile(">\d+_")
        m = prog.search(seqname)
        s = m.group(0)
    result = s.group(0)
    return result

#TO BE DELETED
def getSequences(mypath):
    files = [f for f in os.listdir(mypath) if os.isfile(os.join(mypath, f))]
    seqnamelist = list()

    for filename in files:
        filelocation = '/'.join([mypath, filename])
        with open(filelocation, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    seqnamelist.append(line)
    return seqnamelist

def readCCFjson(jsonfile,evidence_threshold,buildname):
    '''
    read json file from CCF output and convert it to a fasta file with the correct order and array + spacer IDs
    send the spacers to reconstructArray for phylogeny purposes later on
    '''

    conservationDRs = []
    ATcontent = []
    spacerfasta = '../spacers/' + buildname + '.fasta'
    finalspacerdata = []
    onlyspacers = []
    prevOrientation = "ND"
    with open(jsonfile,'r') as j:
        data = json.load(j)
        for seq in data["Sequences"]:
            if len(seq["Crisprs"]) > 0:
                counter = 0
                arrayLength = len(seq["Crisprs"])
                for crispr in seq["Crisprs"]:
                    if crispr["Evidence_Level"] >= evidence_threshold:
                        conservationDRs.append(crispr["Conservation_DRs"])
                        ID = crispr["Name"]
                        if crispr["Potential_Orientation"] == '+':
                            nspacers = 1
                            narrays = arrayLength - counter
                            prevOrientation = '+'
                        elif crispr["Potential_Orientation"] == '-':
                            nspacers = crispr["Spacers"]
                            narrays = 1 + counter
                            prevOrientation = '-'
                        else:
                            if prevOrientation == "+":
                                nspacers = 1
                                narrays = arrayLength - counter
                            elif prevOrientation == '-':
                                nspacers = crispr["Spacers"]
                                narrays = 1 + counter
                            else:
                                nspacers = 1
                                narrays = arrayLength - counter
                        for reg in crispr["Regions"]:
                            if reg["Type"] == "LeftFLANK" and reg["Leader"] == 1:
                                ATcontent.append(reg["AT"])
                            elif reg["Type"] == "RightFLANK" and reg["Leader"] == 1:
                                ATcontent.append(reg["AT"])
                            elif reg["Type"] == "Spacer":
                                if crispr["Potential_Orientation"] == '+' or prevOrientation == '+':
                                    #sequence reverse complementair maken en schrijven
                                    revseq = reverseComplement(reg["Sequence"])
                                    line = '>' + ID + '_arrayID_' + str(narrays) + '_spacerID_' + str(nspacers) + '\n'+ revseq + '\n'
                                    finalspacerdata.insert(0,line)
                                    onlyspacers.insert(0,revseq)
                                    nspacers += 1
                                else:
                                    line = '>' + ID + '_arrayID_' + str(narrays) + '_spacerID_' + str(nspacers) + '\n'+ reg["Sequence"] + '\n'
                                    finalspacerdata.append(line)
                                    onlyspacers.append(reg["Sequence"])
                                    nspacers -= 1
                        counter += 1

    j.close()

    with open(spacerfasta, 'w') as w:
        for line in finalspacerdata:
            w.write(line)
    w.close()

    #make the superspacers for phylogeny
    reconstructArray(onlyspacers, buildname)

    return conservationDRs, ATcontent, spacerfasta

def reconstructArray(spacers, buildname):
    loc = '../reconstructedArray/' + buildname + '_CRISPRarray.fasta'
    #DR = 'N' * 30
    fas = '>' + buildname + "_crispr\n"
    for s in spacers:
        fas += s #+ DR
    with open(loc,'w') as w:
        w.write(fas)
    w.close()

def makeSinglefile():
    outloc = "../phylogeny/allsequences.fasta"
    infiles = getInputfiles("../reconstructedArray/")

    with open(outloc,'w') as w:
        for f in infiles:
            with open(f,'r') as r:
                data = r.read()
                w.write(data + "\n")
            r.close()
    w.close()

def reverseComplement(sequence):
    '''
    get the reverse complement of a sequence
    '''
    reversedsequence = ''
    nucs = list(sequence)
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
    return revnucs

## can be removed
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
##
def addBuildName(file,newInPath):
    # add the filename in front of each fasta sequence of the input.
    # copies inputfile and returns new file and path
    buildname = file.strip('.fasta').split('/')[-1] #get filename without extension
    newfilename = 'modified_' + buildname + '.fasta'
    newpath = newInPath + newfilename
    with open(file, 'r') as f:
        with open(newpath, 'w') as w:
            for line in f:
                if line.startswith(">"):
                    id = getContigNode(line)
                    line = ">" + buildname + '_' + id + "\n"
                    #apparently, CCF can\'t handle parentheses and dots
                    line = re.sub('[().]',"",line)
                    w.write(line)
                else:
                    w.write(line)
        w.close()
    f.close()
    return newpath,buildname

def runCrisprCasFinder(input, output, min, max,cas):
    if cas:
        try:
            os.system("perl /usr/bin/CRISPRCasFinder.pl -i " + str(input) + " -cas -q -so /opt/vmatch-2.3.0/SELECT/sel392.so -outdir " + str(output) + " -minDR " + str(min) + " -maxDR " + str(max))
            print("Done with CCF")
        except Exception as e:
            print("nope")
            print(str(e))
            sys.exit(2)
    else:
        try:
            os.system("perl /usr/bin/CRISPRCasFinder.pl -i " + str(input) + " -q -so /opt/vmatch-2.3.0/SELECT/sel392.so -outdir " + str(output) + " -minDR " + str(min) + " -maxDR " + str(max))
            print("Done with CCF")
        except Exception as e:
            print("nope")
            print(str(e))
            sys.exit(2)

##can be deleted if JSON works
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

##can be deleted if JSON works
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
        os.makedirs(str(outputlocation + "/reconstructedArray"))
        os.makedirs(str(outputlocation + "/phylogeny"))
    else:
        print("Directories already exist")

def getInputfiles(input):
    fastalist = []
    if os.path.isdir(input):
        for file in os.listdir(input):
            if file.endswith(".fasta"):
                fastalist.append(os.path.join(input,file))
    else:
        #just get the single file that is given
        if not input.endswith('.fasta'):
            sys.exit("Not a valid fasta file. Check input argument.")
        else:
            fastalist.append(input)
    return fastalist

def makePlots(DRcons,ATstats,spacerlengths):
    fig, axs = plt.subplots(3, 1, constrained_layout=True)
    axs[0].hist(DRcons,bins = range(0,100,1))
    axs[0].set_title('Direct repeat conservation')
    axs[0].set_xlabel('Percentage conserved (%)')
    axs[0].set_ylabel('Counts')
    fig.suptitle('CRISPR result', fontsize=16)

    axs[1].hist(ATstats,bins = range(0,100,1))
    axs[1].set_xlabel('Percentage AT (%)')
    axs[1].set_title('AT content of leader sequence')
    axs[1].set_ylabel('Counts')

    axs[2].hist(spacerlengths,bins = range(0,100,1))
    axs[2].set_xlabel('Length of spacers')
    axs[2].set_title('Length of spacers')
    axs[2].set_ylabel('Counts')

    plt.show()

def makeTsv(buildnames, DRcons, ATstats, spacerlengths,outpath):
    dat = {'build': buildnames, 'DRconservation': DRcons, 'ATcontent': ATstats, 'spacerlength': spacerlengths}
    df = pd.DataFrame(dat)
    fn = outpath + '/summary.csv'
    df.to_csv(fn, sep = '\t')


def updateCRISPRdb(speciesid, narrays, nspacers, spacerlist):
    con = db.db_connect()
    cur = con.cursor()




def main(args):
    '''
    run CRISPRCasFinder
    run blastn script
    '''
    #check if output path exists. If not, create it
    makeDirs(args.output)

    #runs batch of contigs/scaffolds if dir is given
    curcwd = os.getcwd()
    fastalist = getInputfiles(args.input)
    if len(fastalist) == 0:
        sys.exit("Given directory is empty or does not exist. Check input argument.")

    #loop to run single or multiple files. do this in parallel (maybe) later
    absBVDB = os.path.abspath(args.blastviraldb)
    absOutput = os.path.abspath(args.output)
    outpath = absOutput + '/CCF/'
    revInvPath = absOutput + '/rev_inv/'
    os.makedirs(outpath)
    os.makedirs(revInvPath)

    #empty list for hist sequence lengths
    ATstats = []
    DRcons = []
    spacerlengths = []
    arraylengths = []
    buildnames = []

    for fastafile in fastalist:
        os.chdir(curcwd)
        newpath, buildname = addBuildName(fastafile,revInvPath)
        absInput = os.path.abspath(newpath)
        if args.reverse:
            try:
                print("Fetching reverse complementary fastas ...")
                getReverseComplement(absInput)
            except Exception as e:
                print('Oops. See error:\n' + str(e))
        print("\nStart CCF ...")
        os.chdir(outpath)
        runCrisprCasFinder(absInput,buildname,args.minimum,args.maximum,args.runcas)
        conservationDRs, ATcontent, spacerfasta = readCCFjson(buildname + '/result.json',args.evidencethreshold,buildname)
        DRcons += conservationDRs
        ATstats += ATcontent
        buildnames = buildnames + [buildname] * len(conservationDRs)
        absspacer = os.path.abspath(os.path.join('.',spacerfasta))
        seqlengths, arrayLength = runLocalBlast.seqlength(absspacer)
        spacerlengths += seqlengths
        arraylengths.append(arrayLength)
        blastout = absOutput  + "/blastout/" + buildname + "blastout.out"
        print("\nRun blast ...")
        runLocalBlast.runBlast(absspacer,absBVDB,blastout,args.percidentity)
        print("\nDone with blast")
    makeSinglefile()

    #make plot of spacerlengths
    if args.statistics:
        #runLocalBlast.seqPlot(spacerlengths,args.minimum,args.maximum)
        #for testing
        print(buildnames)
        print(DRcons)
        print(ATstats)
        print(spacerlengths)
        print(arraylengths)
        #/for testing
        #makeplots(DRcons,ATstats,spacerlengths)
        #makeTsv(buildnames, DRcons, ATstats, spacerlengths,absOutput)

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
        parser.add_argument('-stats','--statistics',help = 'make plots with overview of data', default = False, action = 'store_true')
        parser.add_argument('-et','--evidencethreshold', help = 'minimum evidence level for crispr detection. 1 - 4, higher is stricter.',type = int ,default = 4,choices = range(1,5,1))
        parser.add_argument('-cas', '--runcas' , help = 'additional search for cas genes', default = False, action = 'store_true')
        args = parser.parse_args()

        main(args)
    except Exception as e:
        print("error: " + str(e))
        print("See -h for help")
        sys.exit(1)
