from edge import edge

class node :
    def __init__(self, name):
        self.name     = name
        self.c        = 0
        # self.pathnum  = 3
        self.B_sideRec = None
        self.enable   = True
        self.upEdge   = []
        self.downEdge = []
        self.catEdge  = []
        self.upEdgePar   = {}
        self.downEdgePar = {}
        self.pin      = 0
        self.pout     = 0
        self.level    = 0
        self.mark     = 0
        self.pro_m    = [False,False,False,False,False,False,False]
        # self.sidemark = 0
        self.label    = set()
        self.layer    = 0
        # self.paths    = [[]]*(self.pathnum)
        self.pathA1   = []
        self.pathA2   = []
        # self.pathB    = []
        self.pathOut  = []
        self.path     = []

        self.sidepath = []
        self.recordA  = []

        self.labelA      = ""
        self.labelB      = ""
        self.labelB_side = []
        # self.labelC      = []
        # self.labelC_side = []
        # self.labelD      = ""
        # self.labelE      = []

    # def getUpNode(self):
    #     for rec in self.upEdge:
    #         pass
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

    def SetPath(self, startPro, downRec, pathlist):
        tmpList = []
        for path in pathlist:
            tmp = {}
            tmp["struc"] = path["struc"]
            tmp["related"] = path["related"].copy()
            tmp["pathlist"] = path["pathlist"].copy()
            tmp["pathnode"] = path["pathnode"].copy()
            tmp["pathenz"] = path["pathenz"].copy()
            for i in downRec.getrea():
                if i != startPro:
                    tmp["related"].add(i)
            # for i in downRec.getpro():
            #     if i != self:
            #         tmp["related"].add(i)

            tmp["pathnode"].append(startPro)
            tmp["pathenz"].append(downRec.getenz())
            tmp["pathlist"].append(downRec)
            tmpList.append(tmp)
            # assert tmp["pathlist"] != path["pathlist"]
            # assert tmp["pathnode"] != path["pathnode"]
        self.path = tmpList

    def enz_count(self, path):
        count = 0
        for rec in path["pathlist"]:
            count = count + len(rec.getenz().getCatedge())
        return count

    def AddPath(self, startPro, downRec, pathlist):
        path = pathlist[0]
        path_c = self.enz_count(path)
        for path_tmp in pathlist:
            if path_c >= self.enz_count(path_tmp):
                path = path_tmp
                path_c = self.enz_count(path)
        tmp = {}
        tmp["struc"] = path["struc"]
        tmp["related"] = path["related"].copy()
        tmp["pathlist"] = path["pathlist"].copy()
        tmp["pathnode"] = path["pathnode"].copy()
        tmp["pathenz"] = path["pathenz"].copy()
        for i in downRec.getrea():
            if i != startPro:
                tmp["related"].add(i)
        # for i in downRec.getpro():
        #     if i != self:
        #         tmp["related"].add(i)

        tmp["pathnode"].append(startPro)
        tmp["pathenz"].append(downRec.getenz())
        tmp["pathlist"].append(downRec)
        # assert tmp["pathlist"] != []
        # assert tmp["pathnode"] != path["pathnode"]
        # assert tmp != {}
        if len(self.path) >= 10:
            for path in self.path:
                if len(path["related"]) >= len(tmp["related"]):
                    path = tmp
                    break
        else:
            self.path.append(tmp)

    def CheckRecSpecies(self, pathnode, downRec):
        for species in (downRec.getrea()):
            if species in pathnode:
                return False
        for species in (downRec.getpro()):
            if species in pathnode:
                return False
        if downRec.getenz() in pathnode:
            return False
        return True

    def CheckEnz(self, pathnode, downRecs):
        for rec in downRecs:
            for species in (rec.getpro()):
                if species in pathnode:
                    return False
        return True

    def CheckDownRec(self, downRec):
        tmp = []
        if self in downRec.getpro():
            return []

        # if downRec.getenz().mark == 1:
        #     return []

        # for rea in downRec.getrea():
        #     if rea.mark == 1:
        #         return []

        downRecs = []
        downRecs.append(downRec)
        added_species = set()
        for path in self.path:
            pathenz = path["pathenz"].copy()
            pathenz.append(downRec.getenz())

            added_species = path["related"].copy()
            for enz in path["pathenz"]:
                added_species.add(enz)
            for enz in pathenz:
                for rec in enz.getCatedge():
                    if rec.activated(added_species):
                        downRecs.append(rec)

        for path in self.path:
            if (self.CheckRecSpecies(path["pathnode"], downRec) and self.CheckEnz(path["pathnode"], downRecs)):
                if downRec not in path["pathlist"]:
                    tmp.append(path)
        return tmp


    # check these selected pathnode products whether they have been related species already or not.
    def CheckProduct(self, pathlist, up_limit, outNode):
        tmp = []
        for path in pathlist:
            if self not in path["related"] and self.name != "H2O":
                tmp.append(path)
        return tmp

    def CopyToRecord(self, pathtype):
        path_tmp = []
        for path in self.path:
            tmp    = {}
            tmp["struc"] = path["struc"]
            tmp["related"]  = path["related"].copy()
            tmp["pathlist"] = path["pathlist"].copy()
            tmp["pathenz"]  = path["pathenz"].copy()
            tmp["pathnode"] = path["pathnode"].copy()
            path_tmp.append(tmp)
        
        if pathtype == "A1":
            self.pathA1 = path_tmp
        elif pathtype == "A2":
            self.pathA2 = path_tmp
        elif pathtype == "A":
            self.recordA = path_tmp

    def CopyToPath(self, pathtype):
        path = []
        if pathtype == "A":
            for record in self.recordA:
                tmp    = {}
                tmp["struc"] = record["struc"]
                tmp["related"]  = record["related"].copy()
                tmp["pathlist"] = record["pathlist"].copy()
                tmp["pathenz"]  = record["pathenz"].copy()
                tmp["pathnode"] = record["pathnode"].copy()
                path.append(tmp)
        elif pathtype == "A1":
            for record in self.pathA1:
                tmp    = {}
                tmp["struc"] = record["struc"]
                tmp["related"]  = record["related"].copy()
                tmp["pathlist"] = record["pathlist"].copy()
                tmp["pathenz"]  = record["pathenz"].copy()
                tmp["pathnode"] = record["pathnode"].copy()
                path.append(tmp)
        elif pathtype == "A2":
            for record in self.pathA2:
                tmp    = {}
                tmp["struc"] = record["struc"]
                tmp["related"]  = record["related"].copy()
                tmp["pathlist"] = record["pathlist"].copy()
                tmp["pathenz"]  = record["pathenz"].copy()
                tmp["pathnode"] = record["pathnode"].copy()
                path.append(tmp)
        self.path = path

    def CheckMerge(self, path0, path1):
        if path0 == {} or type(path0) != type({}):
            return False
        if path1 == {} or type(path1) != type({}):
            return False
        if path0["pathlist"] == []:
            return False
        if path1["pathlist"] == []:
            return False
        for species in path1["pathnode"]:
            if species in path0["related"]:
                return False
        for species in path0["pathnode"]:
            if species in path1["related"]:
                return False
        enzyme = set()
        if (path0["pathlist"][0].getenz() == path1["pathlist"][0].getenz()):
            return False
        return True


    def Merge(self, path, record):
        tmp    = {}
        tmp["struc"] = path["struc"]
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

    def MergeAll(self, label):
        # assert label == "A" or label == "B_side"
        path_tmp = []
        c = 0
        for path in self.path:
            for record in self.recordA:
                if self.CheckMerge(path, record):
                    # assert path != {}
                    # assert record["pathlist"] != []
                    path_tmp.append(self.Merge(path, record))
        self.recordA = path_tmp

    def ClearRecordA(self):
        self.recordA.clear()