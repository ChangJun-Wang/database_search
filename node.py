import copy

class node :
    def __init__(self, name):
        path_tmp        = []
        tmp             = {}
        tmp["related"]  = set()
        tmp["pathlist"] = []
        tmp["pathenz"]  = []
        tmp["pathnode"] = []
        path_tmp.append(tmp)
        self.name     = name
        self.enable   = True
        self.upEdge   = []
        self.downEdge = []
        self.catEdge  = []
        self.upEdgePar   = {}
        self.downEdgePar = {}
        self.pin      = 0
        self.pout     = 0
        self.level    = 0
        self.label    = set()
        self.path     = path_tmp
        self.sidepath = copy.deepcopy(path_tmp)

        # self.path     = {}
        # self.pathID   = []
        self.recordA  = copy.deepcopy(path_tmp)
        self.recordC  = {}
        # self.path_tmp = []
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

    def AddPath(self, startPro, downRec, pathlist):
        tmpList = []
        for path in pathlist:
            tmp = {}
            tmp["related"] = path["related"].copy()
            for i in downRec.getrea():
                if i != startPro:
                    tmp["related"].add(i)
            tmp["pathlist"] = path["pathlist"].copy()
            tmp["pathlist"].append(downRec)
            tmp["pathnode"] = path["pathnode"].copy()
            tmp["pathnode"].append(startPro)
            tmp["pathenz"] = path["pathenz"].copy()
            tmp["pathenz"].append(downRec.getenz()[0])
            tmpList.append(tmp)
            assert tmp["pathlist"] != path["pathlist"]
            assert tmp["pathnode"] != path["pathnode"]
            assert tmp != path
        self.path = tmpList

    def CheckRecSpecies(self, pathnode, downRec):
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

    def CheckEnz(self, pathnode, downRecs):
        for rec in downRecs:
            for species in (rec.getpro()):
                if species in pathnode: # or species.name == "NAD+":
                    return False
        return True

    def CheckDownRec(self, downRec):
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
            if self.CheckRecSpecies(path["pathnode"], downRec) and self.CheckEnz(path["pathnode"], downRecs):
                tmp.append(path)
        return tmp


    # check these selected pathnode products whether they have been related species already or not.
    def CheckProduct(self, pathlist):
        tmp = []
        for path in pathlist:
            if self not in path["related"]:
                tmp.append(path)
        return tmp

    def CopyToRecord(self, label, num):
        assert label == "A" or label == "C_side"
        tmpList = []
        if label == "A":
            for path in self.path:
                tmp    = {}
                tmp["related"]  = path["related"].copy()
                tmp["pathlist"] = path["pathlist"].copy()
                tmp["pathenz"]  = path["pathenz"].copy()
                tmp["pathnode"] = path["pathnode"].copy()
                tmpList.append(tmp)
            self.recordA = tmpList
        elif label == "C_side":
            pass

    def CopyToPath(self):
        # assert label == "A" or label == "C_side"
        # if label == "A":
        tmpList = []
        for record in self.recordA:
            tmp    = {}
            tmp["related"]  = record["related"].copy()
            tmp["pathlist"] = record["pathlist"].copy()
            tmp["pathenz"]  = record["pathenz"].copy()
            tmp["pathnode"] = record["pathnode"].copy()
            tmpList.append(tmp)
        self.path = tmpList
        # elif label == "C_side":
        #     pass

    def CopyToSide(self):
        tmpList = []
        for path in self.path:
            tmp    = {}
            tmp["related"]  = path["related"].copy()
            tmp["pathlist"] = path["pathlist"].copy()
            tmp["pathenz"]  = path["pathenz"].copy()
            tmp["pathnode"] = path["pathnode"].copy()
            tmpList.append(tmp)
        self.sidepath = tmpList

    def CheckMerge(self, path0, path1):
        for species in path1["pathnode"]:
            if species in path0["related"]:
                return False
        for species in path0["pathnode"]:
            if species in path1["related"]:
                return False
        return True


    def Merge(self, path, record):
        tmp    = {}
        tmp["related"]  = path["related"].copy()
        tmp["pathlist"] = path["pathlist"].copy()
        tmp["pathenz"]  = path["pathenz"].copy()
        tmp["pathnode"] = path["pathnode"].copy()

        for species in record["related"]:
            tmp["related"].add(species)
        for rec in record["pathlist"]:
            tmp["pathlist"].append(rec)
        for enz in record["pathenz"]:
            tmp["pathenz"].append(enz)
        for species in record["pathnode"]:
            tmp["pathnode"].append(species)

        return tmp

    def MergeAll(self, label, num):
        assert label == "A" or label == "C_side"
        if label == "A":
            path_tmp = []
            for path in self.path:
                for record in self.recordA:
                    if self.CheckMerge(path, record):
                        path_tmp.append(self.Merge(path, record))
            self.recordA = path_tmp
            # if path_tmp != []:
            return True
        elif label == "C_side":
            pass
        return False
