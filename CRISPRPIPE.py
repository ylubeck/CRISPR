#interactive program run
from tkinter import filedialog
from tkinter import *

def selectInput():
    input = filedialog.askopenfilenames(title = "Select inputfile(s)")
    print(input)
    return input

def selectOutput():
    output = filedialog.askdirectory(title = "Select a directory for the output")
    print(output)
    return output

def selectBlastDB():
    DB = filedialog.askdirectory(title = "Select BLAST database directory")
    print(DB)
    return DB

def var_states():
    a = str(runblast.get())
    b = str(dostats.get())
    c = str(minDR.get())
    d = str(maxDR.get())
    e = str(evidence.get())
    print(" ".join([a,b,c,d,e]))

def run_crispr_cas_finder():
    #handler which will receive vars from the GUI
    pass

def openREADME():
    print("placeholder for readme")

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

    minDR = Entry(root)
    maxDR = Entry(root)
    evidence = Entry(root)

    minDR.insert(END,23)
    maxDR.insert(END,55)
    evidence.insert(END,4)

    minDR.grid(row = 2, column = 1)
    maxDR.grid(row = 3, column = 1)
    evidence.grid(row = 4, column = 1)

    #for now, check input
    Button(root, text = 'Run program', command = var_states).grid(row = 5, sticky = W)

    root.config(menu = menu)
    root.mainloop()
