#class for args CrisprCasFinderGui

class GUIargs:
    def __init__(self,input,output,reverse,viralfastas,blastviraldb,minimum,maximum,percidentity,statistics,evidencethreshold):
        self.input = input
        self.output = output
        self.reverse = reverse
        self.viralfastas = viralfastas
        self.blastviraldb = blastviraldb
        self.minimum = minimum
        self.maximum = maximum
        self.percidentity = percidentity
        self.statistics = statistics
        self.evidencethreshold = evidencethreshold

    def printVars(self):
        print("input: " + str(self.input))
        print("output: " + str(self.output))
        print("reverse: "  + str(self.reverse))
        print("stats: " + str(self.statistics))
