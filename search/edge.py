class edge :
    def __init__(self, name):
        self.name     = name
        self.enable   = True
        # self.tmp_path = []
        self.toRea    = set()  #
        self.toPro    = set()  #
        self.toEnz    = ""
        self.reverse  = False
        self.label    = set()
        self.labelC   = []
        self.labelB   = []
        self.visited  = 0
        
    def getpro(self):
        return self.toPro
    def getrea(self):
        return self.toRea
    def getenz(self):
        return self.toEnz

    def getall(self):
        toEnz = []
        toEnz.append(self.toEnz)
        return self.toRea + self.toPro + toEnz

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

    def show(self):
        temp = "\n" + "reaction: "
        temp += (str(self.name)) + " \n"
        temp += "k"
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