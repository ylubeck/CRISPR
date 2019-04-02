# database
import re
import sqlite3
import os
import argparse
import sys
import csv

DEFAULT_PATH = os.path.join(os.path.dirname(__file__), 'CRISPRCASDB.sqlite3')

CRISPR_sql = """
CREATE TABLE IF NOT EXISTS crispr (
    id integer PRIMARY KEY,
    species_name text NOT NULL,
    n_arrays integer NOT NULL,
    n_spacers integer NOT NULL
)"""

spacers_sql = """
CREATE TABLE IF NOT EXISTS spacers (
    id integer PRIMARY KEY,
    spacer text NOT NULL,
    array_id integer NOT NULL,
    spacer_id integer NOT NULL,
    species_id integer,
    FOREIGN KEY (species_id) REFERENCES crispr (id)
)"""

def db_connect(db_path=DEFAULT_PATH):
    con = sqlite3.connect(db_path)
    return con

def addCRISPR(con, species_id, n_arrays, n_spacers):
    cur = con.cursor()
    insertCRISPR = "INSERT INTO crispr (species_name, n_arrays, n_spacers) VALUES (?, ?, ?)"
    cur.execute(insertCRISPR, (str(species_id),int(n_arrays),int(n_spacers)))
    return cur.lastrowid

def addSpacer(con, spacer, array_id, spacer_id, species_id):
    cur = con.cursor()
    insertSpacer = "INSERT INTO spacers (spacer, array_id, spacer_id, species_id) VALUES (?, ?, ?, ?)"
    cur.execute(insertSpacer, (str(spacer),int(array_id),int(spacer_id), int(species_id)))
    return cur.lastrowid

def setTables(con):
    cur = con.cursor()
    cur.execute(CRISPR_sql)
    cur.execute(spacers_sql)

def readFile(file):
    #read fasta spacer file from runCRISPRCASFinder
    prog1 = re.compile("arrayID_\d+")
    prog2 = re.compile("spacerID_\d+")

    speciesid = file.strip('.fasta').split('/')[-1]
    spacer = ''
    arrayid = 0
    spacerid = 0
    prev_arrayid = 0

    arraycounter = 0
    spacercounter = 0

    datadict = {}
    spacerdict = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                m1 = prog1.search(line)
                m2 = prog2.search(line)
                arrayid = int(m1.group(0).split("_")[1])
                spacerid = int(m2.group(0).split("_")[1])
                spacercounter += 1
                spacer = next(f).strip()

                if arrayid != prev_arrayid:
                    arraycounter += 1
                    prev_arrayid = arrayid
                    spacerdict = {}
                spacerdict[spacerid] = spacer
                datadict[arrayid] = spacerdict
    f.close()

    return speciesid, datadict, arraycounter,spacercounter

def updateDatabase(files, con, checkpoint):
    #loop over files in directory, add to db
    for f in files:
        speciesid, datadict, narrays, nspacers = readFile(f)
        #add checkpoint for double speciesid
        if checkpoint:
            cur = con.cursor()
            cur.execute("SELECT id FROM crispr where species_name = ?", (speciesid,))
            data = cur.fetchone()
            if data is not None:
                answer = input(speciesid + " already exists! Ignore warning and add record? [Y/n/quit]")
                if answer == 'quit':
                    sys.exit("Database update aborted")
                if answer not in ['Y','y']:
                    continue
        print('Adding data from:\n' + str(f) + "\n")
        ID = addCRISPR(con, speciesid, narrays, nspacers)
        for k1, v1 in datadict.items():
            for k2, v2 in v1.items():
                ID2 = addSpacer(con, v2, k1, k2, ID)


#testen!
def exportDB(con):
    cur = con.cursor()

    fspecies = 'species.tsv'
    fspacers = 'spacers.tsv'

    cur.execute("SELECT id, species_name, n_arrays, n_spacers FROM crispr")
    results = cur.fetchall()
    writeTable(results,fspecies)

    cur.execute("SELECT id, spacer, array_id, spacer_id, species_id FROM spacers")
    results = cur.fetchall()
    writeTable(results,fspacers)

#testen!
def writeTable(results,tablefile):
    with open(tablefile, 'w') as w:
        csv.register_dialect("custom", delimiter = '\t', skipinitialspace = True)
        writer = csv.writer(w, dialect = "custom")
        for row in results:
            writer.write(row)
    w.close()

if __name__ ==  '__main__':
    con = db_connect()
    setTables(con)
    cur = con.cursor()
    #for testing
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', help = 'input file or directory to add data. \nThese files should be located in the spacer folder in the output directory.', required = True)
        parser.add_argument('-c', '--checkpoint', help = 'check if species id already exists in database.', default = False, action = 'store_true')
        parser.add_argument('-e', '--export', help = 'export tables to TSV-files', default = False, action = 'store_true')
        args = parser.parse_args()
        filespath = "/data/lubecky/CRISPRCAS/scripts/salmonella_large/spacers/"
        testfiles = [os.path.join(filespath,x) for x in os.listdir(filespath)]
        updateDatabase(testfiles,con,args.checkpoint)
        if args.export:
            exportDB(con)
        con.commit()
    except:
        con.rollback()
        raise RuntimeError("Foek")

    cur.execute("SELECT id, species_name, n_arrays, n_spacers FROM crispr")
    results = cur.fetchall()
    for row in results:
        print(row)

    cur.execute("SELECT id, spacer, array_id, spacer_id, species_id FROM spacers")
    results = cur.fetchall()
    for row in results:
        print(row)
