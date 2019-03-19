#!/usr/bin/env python3
#interactive program run
from tkinter import filedialog
from tkinter import *
from guiargs import GUIargs
import runCRISPRCASFinder

input_path = ''
output_path = ''
DB_path = ''

def selectInput():
    global input_path
    input_path = filedialog.askopenfilenames(title = "Select inputfile(s)")

def selectOutput():
    global output_path
    output_path = filedialog.askdirectory(title = "Select a directory for the output")

def selectBlastDB():
    global DB_path
    DB_path = filedialog.askdirectory(title = "Select BLAST database directory")

def var_states():
    ''' For testing. Can be removed when GUI works '''
    a = str(runblast.get())
    b = str(dostats.get())
    c = str(minDR.get())
    d = str(maxDR.get())
    e = str(evidence.get())

    args = GUIargs(input_path,output_path,False,"vfdb",DB_path,c,d,95,b,e)

    args.printVars()

def run_crispr_cas_finder():
    #handler which will receive vars from the GUI

    pass

def openREADME():
    print("placeholder for readme")

def showProgress():
    window = Tk()
    window.title("Progress...")
    pbar = tkinter.ttk.Progressbar(window, length = 100)


if __name__ == '__main__':
    root = Tk()
    menu = Menu(root)

    #top menu bar
    filemenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = "Start", menu = filemenu)
    filemenu.add_command(label = "Select inputfile(s)...", command = selectInput)
    filemenu.add_command(label = "Select outputfolder...", command = selectOutput)
    filemenu.add_separator()
    filemenu.add_command(label = "Quit", command = root.quit)

    editmenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = "Edit", menu = editmenu)
    editmenu.add_command(label = "Change BLAST Database...", command = selectBlastDB)

    helpmenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = "Help", menu = helpmenu)
    helpmenu.add_command(label = "help...", command = openREADME)

    runblast = IntVar()
    dostats = IntVar()
    Checkbutton(root, text = "Run blast on spacers", variable = runblast).grid(row = 0, sticky = W)
    Checkbutton(root, text = "Histograms of output", variable = dostats).grid(row = 0,column = 1, sticky = W)

    #entry boxes
    Label(root, text = "Minimum length CRISPR DRs").grid(row = 2)
    Label(root, text = "Maximum length CRISPR DRs").grid(row = 3)
    Label(root, text = "Evidence threshold CRISPRs").grid(row = 4)
    Label(root, text = "\'%\' sequence identity spacers BLAST").grid(row = 5)
    minDR = Entry(root)
    maxDR = Entry(root)
    evidence = Entry(root)
    percID = Entry(root)

    minDR.insert(END,23)
    maxDR.insert(END,55)
    evidence.insert(END,4)
    percID.insert(END,95)

    minDR.grid(row = 2, column = 1)
    maxDR.grid(row = 3, column = 1)
    evidence.grid(row = 4, column = 1)
    percID.grid(row = )

    #for now, check input
    Button(root, text = 'Run program', command = var_states).grid(row = 5, sticky = W)

    root.config(menu = menu)
    root.mainloop()
