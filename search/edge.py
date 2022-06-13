class edge :
    def __init__(self, name):
        self.name     = name
        self.enable   = True
        # self.tmp_path = []
        self.toRea    = set()  #
        self.toPro    = set()  #
        self.toEnz    = None
        self.reverse  = False
        self.label    = set()
        self.labelC   = []
        self.labelB   = []
        self.visited  = 0
        self.pathA1   = None
        
    def getpro(self):
        return self.toPro
    def getrea(self):
        return self.toRea
    def getenz(self):
        return self.toEnz

    def getall(self):
        rec = set()
        rec.add(self.toEnz)
        for i in self.toRea:
            rec.add(i)
        for i in self.toPro:
            rec.add(i)
        return rec

    def addpro(self, pro):
        (self.toPro).add(pro)
    def addrea(self, rea):
        (self.toRea).add(rea)
    def setenz(self, enz):
        (self.toEnz) =  (enz)

    def activated(self, added_species):
        if self.toEnz not in added_species:
            return False
        for rea in self.toRea:
            if rea not in added_species:
                return False
        return True

    def show(self, c):
        temp = "\n" + "reaction: "
        temp += (str(self.name)) + " \n"
        temp += "k" + str(c)
        for i in self.getrea():
            temp +=  "*" + (str(i.name))
        temp += "*" + self.getenz().name
        temp += " for "
        for i in self.getrea():
            temp +=  (str(i.name)) + " + "
        temp +=  " =[ " + self.getenz().name + " ]=> "
        for i in self.getpro():
            temp += (str(i.name)) + " + "
        temp += "."
        return temp

    def is_reverse(self):
        allnodes = self.toPro.copy()
        allnodes.add(self.toEnz)
        for rec in self.toEnz.getCatedge():
            if rec.activated(allnodes):
                return True
        for node in self.toPro:
            if "output" in node.label:
                return True
        return False
