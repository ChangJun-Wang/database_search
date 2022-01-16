class edge :
    def __init__(self, name):
        self.name     = name
        self.enable   = True
        # self.tmp_path = []
        self.toRea    = set()  #
        self.toPro    = set()  #
        self.toEnz      = []
        self.reverse  = False
        self.label    = set()
        self.labelC   = []
        self.visited  = 0
        # self.tmp_vis  = 0
        # self.recStack = False
    def getpro(self):
        return self.toPro
    def getrea(self):
        return self.toRea
    def getenz(self):
        return self.toEnz

    def addpro(self, pro):
        (self.toPro).add(pro)
    def addrea(self, rea):
        (self.toRea).add(rea)
    def addenz(self, enz):
        (self.toEnz).append(enz)

    def activated(self, added_species):
        if self.toEnz[0] not in added_species:
            return False
        for rea in self.toRea:
            if rea not in added_species:
                return False
        return True

    def show(self):
        temp = "\n" + "reaction: "
        temp += (str(self.name)) + " \n"

        for i in self.getrea():
            temp +=  (str(i.name)) + " "
        for i in self.getenz():
            temp +=  " =[ " + (str(i.name)) + " ]=> "
        for i in self.getpro():
            temp += (str(i.name)) + " "
        return temp