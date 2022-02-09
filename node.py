from path import path

class node :
    def __init__(self, name):
        self.name     = name
        self.enable   = True
        self.upEdge   = []
        self.downEdge = []
        self.catEdge  = []
        self.upEdgePar   = {}
        self.downEdgePar = {}
        self.pin      = 0
        self.pout     = 0
        self.visited  = 0
        self.label    = set()
        self.path     = []
        self.path_tmp = []
        self.labelA   = ""
        self.labelB   = ""
        self.labelC   = []
        self.labelC_side = []
        self.labelD   = ""
        self.labelE   = []
    def getUpedge(self):
        return self.upEdge
    def getDownedge(self):
        return self.downEdge
    def getCatedge(self):
        return self.catEdge
    def getpin(self):
        return self.pin
    def getpout(self):
        return self.pout

    def setpin(self, pin):
        self.pin = pin
    def setpout(self, pout):
        self.pout = pout

    def addpin(self):
        self.pin = self.pin + 1
    def addpout(self):
        self.pout = self.pout + 1
    def addUpedge(self, rec):
        (self.upEdge).append(rec)
    def addDownedge(self, rec):
        (self.downEdge).append(rec)
    def addCatedge(self, rec):
        (self.catEdge).append(rec)

    def show(self):
        return (self.name)

    def add_path(self, startPro, downRec, pathlist):
        for path in pathlist:
            tmp_path = {}
            tmp_path["related"] = path["related"].copy()
            for i in downRec.getrea():
                tmp_path["related"].add(i)
            tmp_path["pathlist"] = path["pathlist"].copy()
            tmp_path["pathlist"].append(downRec)
            tmp_path["pathnode"] = path["pathnode"].copy()
            tmp_path["pathnode"].append(startPro)
            tmp_path["pathenz"] = path["pathenz"].copy()
            tmp_path["pathenz"].append(downRec.getenz()[0])
            self.path.append(tmp_path)
            assert tmp_path["pathlist"] != path["pathlist"]
            assert tmp_path["pathnode"] != path["pathnode"]
            assert tmp_path != path

    def check_RecSpecies(self, pathnode, downRec):
        for species in (downRec.getrea()):
            if species in pathnode:
                return False
        for species in (downRec.getpro()):
            if species in pathnode:
                return False
        for species in (downRec.getenz()):
            if species in pathnode:
                return False
        return True

    def check_EnzSpecies(self, pathnode, downRecs):
        for rec in downRecs:
            for species in (rec.getpro()):
                if species in pathnode:
                    return False
        return True

    def check_downRec(self, downRec):
        tmp = []
        downRecs = []
        downRecs.append(downRec)
        added_species = set()
        for path in self.path:
            pathenz = path["pathenz"].copy()
            pathenz.append(downRec.getenz()[0])

            added_species = path["related"].copy()
            for enz in path["pathenz"]:
                added_species.add(enz)
            for enz in pathenz:
                for rec in enz.getCatedge():
                    if rec.activated(added_species):
                        downRecs.append(rec)

        for path in self.path:
            if self.check_RecSpecies(path["pathnode"], downRec) and self.check_EnzSpecies(path["pathnode"], downRecs):
                tmp.append(path)
        return tmp


    # check these selected pathnode products whether they have been related species already or not.
    def check_pro(self, pathlist):
        tmp = []
        for path in pathlist:
            if self not in path["related"]:
                tmp.append(path)
        return tmp

    def check_downRectmp(self, downRec):
        tmp = []
        for tmp_path in self.path_tmp:
            if self.check_RecSpecies(tmp_path["pathnode"], downRec):
                tmp.append(tmp_path)
        return tmp
