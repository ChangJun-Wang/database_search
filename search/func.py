from node import node

from edge import edge

# from bfs  import BFS

def BFS(startNode, Target_type):
    search_limit = 1
    bfs          = [startNode]
    stype        = Target_type
    candidate    = set()

    # initializtion for starting node
    bfs[0].CopyToPath()

    # start to traverse
    count = 0
    while (bfs != []):
        startPro = bfs.pop(0)
        if (stype in startPro.label):
            if stype == "A" and (startPro.getpin() > 2):
                # print ("find A")
                foundA(startPro, candidate)
            elif stype == "C_side":
                foundC_side(startPro, candidate)
            elif stype == "B":
                foundB(startPro, candidate)
            elif stype == "B_side":
                foundB_side(startPro, candidate)
        if startPro.level + 1 <= search_limit:
            # for path in startPro.path:
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                CheckDownRec = startPro.CheckDownRec(downRec)
                if CheckDownRec != []:
                    for product in downRec.getpro():
                        CheckProduct = product.CheckProduct(CheckDownRec)
                        if CheckProduct != []:
                            product.AddPath(startPro, downRec, CheckProduct)
                            product.level = startPro.level + 1
                            if (product not in bfs and product.level <= search_limit):
                                bfs.append(product)
                            assert product.path != CheckProduct
                            assert len(product.path) != 0
        # startPro.level = 0
        # startPro.path = []
    return candidate


def foundA(product, candidate):
    if product.path == []:
        return False
    else:
        candidate.add((product, product.labelA))
        return True

def foundC_side(product, candidate):
    for label in product.labelC_side:
        candidate.add((product, label))

def foundB(product, candidate):
    assert(len(product.getDownedge()) > 1)
    for downRec in product.getDownedge():
        candidate.add((product, downRec.name))
        downRec.labelB.append(product.labelB)
        for rea in downRec.getrea():
            if rea != product:
                rea.label.add("B_side")
                # rea.labelB_side.append(product.labelB)
                # rea.labelBmap[product.labelB] = downRec.name

def foundB_side(product, candidate):
    for downRec in product.getDownedge():
        if downRec.labelB!=[]:
            candidate.add((product, downRec.name))

class search:
    def __init__(self):
        self.c = 0
        self.forbid_node = set()
        self.mapToNode   = {}
        self.mapToEdge   = {}
        self.mapToClist  = {}
        self.mapToCnode  = {}
        self.mapToCforbid= {}
        self.nodeList    = []
        self.edgeList    = []
        self.HashTable   = {}
        self.HashCollect = {}
        self.reaction    = 0
        self.typeA       = 0
        self.typeB       = 0
        self.typeC       = 0
        self.typeD       = 0
        self.typeE       = 0

        self.count       = 0

    def split(self, inputline):
        temp    = ""
        tmpList = []
        for w in inputline:
            if (w == " " or w == "\n"):
                if temp != "":
                    tmpList.append(temp)
                temp = ""
            else:
                temp += w
        return tmpList

    def initialize(self):
        for nodes in self.nodeList:
            pathtmp             = {}
            sidetmp             = {}
            rectmp              = {}
            pathtmp["related"]  = set()
            pathtmp["pathlist"] = []
            pathtmp["pathenz"]  = []
            pathtmp["pathnode"] = []
            rectmp["related"]   = set()
            rectmp["pathlist"]  = []
            rectmp["pathenz"]   = []
            rectmp["pathnode"]  = []
            sidetmp["related"]  = set()
            sidetmp["pathlist"] = []
            sidetmp["pathenz"]  = []
            sidetmp["pathnode"] = []
            (nodes.path).append(pathtmp)
            (nodes.sidepath).append(sidetmp)
            (nodes.recordA).append(rectmp)
            nodes.path[0]["related"].add(self.mapToNode["H2O"])
            nodes.path[0]["pathenz"].append(self.mapToNode["spontaneous_reaction"])

            if (nodes.getpin() > 1):
                enz = set()
                for rec in nodes.getUpedge():
                    if rec.getenz() != "spontaneous_reaction":
                        enz.add(rec.getenz())
                self.typeA += 1
                nodes.setpin(len(enz))
                nodes.label.add("A")
                nodes.labelA = str(self.typeA)
            if (nodes.getpout() > 1):
                enz = set()
                for rec in nodes.getDownedge():
                    if rec.getenz() != "spontaneous_reaction":
                        enz.add(rec.getenz())
                self.typeB += 1
                nodes.setpout(len(enz))
                nodes.label.add("B")
                nodes.labelB = str(self.typeB)
            if (nodes.getpout() > 0) and (nodes.getpin() > 0):
                self.typeD += 1
                nodes.label.add("D")
                nodes.labelD = str(self.typeD)
            if (len(nodes.getCatedge()) > 1):
                self.typeE += 1
                nodes.label.add("E")
                nodes.labelE.append(str(self.typeE))
        self.ClearVis()

    def parsing(self, inputfile):
        temp = ""
        tmpSpecies = ""
        tmpList = []
        enzyme = ""
        is_product = 0
        HashKey = "" 

        print("in database parsing")
        for line in inputfile:
            tmpList = self.split(line)

            if (tmpList[0] == "RN"):
                enzyme = tmpList[1]
                if (enzyme not in self.mapToNode):
                    self.mapToNode[enzyme] = node(enzyme)
                    self.nodeList.append(self.mapToNode[enzyme])
            else: #if (tmpList[0] == "RE"):
                is_product = 0
                self.reaction = tmpList[0]
                self.mapToEdge[self.reaction] = edge(self.reaction)
                self.edgeList.append(self.mapToEdge[self.reaction])

                self.mapToEdge[self.reaction].setenz(self.mapToNode[enzyme])
                self.mapToNode[enzyme].addCatedge(self.mapToEdge[self.reaction])

                for species in tmpList[1:]:
                    HashKey = HashKey + species[0]
                    if (species == "="):
                        is_product = 1
                    else:
                        par = 1
                        if len(species) > 2:
                            if species[1] == "_" and species[0].isdigit():
                                par = int(species[0])
                                species = species[2:]
                                
                        if (species not in self.mapToNode):
                            self.mapToNode[species] = node(species)
                            self.nodeList.append(self.mapToNode[species])
                        if (is_product == 0):
                            self.mapToEdge[self.reaction].addrea(self.mapToNode[species])
                            self.mapToNode[species].addDownedge(self.mapToEdge[self.reaction])
                            (self.mapToNode[species]).downEdgePar[self.mapToEdge[self.reaction]] = par
                            self.mapToNode[species].addpout()
                        else:
                            self.mapToEdge[self.reaction].addpro(self.mapToNode[species])
                            self.mapToNode[species].addUpedge(self.mapToEdge[self.reaction])
                            (self.mapToNode[species]).upEdgePar[self.mapToEdge[self.reaction]] = par
                            self.mapToNode[species].addpin()

    def parsing_label(self, inputfile):
        print("in record label parsing")
        tmpList = []
        typeC = 0
        forbid = False
        for line in inputfile:
            tmpList = self.split(line)
            if len(tmpList) > 0:
                if tmpList[0] == "TypeC:":
                    self.typeC += 1
                    forbid = False
                    typeC = tmpList[1]
                    self.mapToClist[typeC] = []
                    assert(tmpList[1].isdigit())
                elif tmpList[0] == "reaction:":
                    self.mapToClist[typeC].append(self.mapToEdge[tmpList[1]])
                    self.mapToEdge[tmpList[1]].label.add("C")
                    self.mapToEdge[tmpList[1]].labelC.append(typeC)
                elif tmpList[0] == "CoreNode:":
                    assert len(tmpList) == 5
                    self.mapToCnode[typeC] = []
                    self.mapToCnode[typeC].append(self.mapToNode[tmpList[1]])
                    self.mapToCnode[typeC].append(self.mapToNode[tmpList[2]])
                    self.mapToCnode[typeC].append(self.mapToNode[tmpList[3]])
                    self.mapToCnode[typeC].append(self.mapToNode[tmpList[4]])
                    self.mapToNode[tmpList[1]].label.add("C")
                    self.mapToNode[tmpList[2]].label.add("C")
                    self.mapToNode[tmpList[3]].label.add("C_side")
                    self.mapToNode[tmpList[4]].label.add("C_side")
                    self.mapToNode[tmpList[1]].labelC.append(typeC)
                    self.mapToNode[tmpList[2]].labelC.append(typeC)
                    self.mapToNode[tmpList[3]].labelC_side.append(typeC)
                    self.mapToNode[tmpList[4]].labelC_side.append(typeC)
                    
                elif tmpList[0] == "Forbidnode:":
                    forbid = True
                    self.mapToCforbid[typeC] = []
                else:
                    if forbid:
                        for i in tmpList:
                            tmp = []
                            tmp.append(self.mapToNode[i])
                        self.mapToCforbid[typeC].append(tmp)

        print("start initializing")
        self.initialize()

    def main(self, input_species):
        candidate = []
        result    = []
        self.BuildStart(input_species)
        candidate0 = BFS(input_species[0], "A")
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord()
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0))

        self.BuildStart(input_species)
        candidate1 = BFS(input_species[1], "A")
        # check if there is any common label found between input0 and input1 path
        print("start merging type A")
        for label in candidate0:
            if label in candidate1:
                if label[0].MergeAll("A", label[1]):
                    candidate.append(label[0])
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate1))
        print(len(candidate))

        for i in range(len(candidate)):
            candidate0 = BFS(candidate[i], "B")
            print ("candidate0 : ", len(candidate0))
            for species in candidate0:
                species[0].CopyToSide()
            self.ClearPath()
            self.ClearLevel()
            for j in range(len(candidate)):
                if i != j:
                    candidate1 = BFS(candidate[j], "B_side")
                    print ("candidate1 : ", len(candidate1))
                    print("start merging type B_side")
                    for label0 in candidate0:
                        for label1 in candidate1:
                            labelB = []
                            if (label0[1] == label1[1]) and (label0[0] != label1[0]):
                                labelB.append(label0[0])
                                labelB.append(label1[0])
                                accept = self.MergeBside(label0, label1)
                                if accept != {}:
                                    allnodes = set()
                                    result.append(accept)
                                    for node in accept["pathnode"]:
                                        print(node.show())
                                    for node in labelB:
                                        print("labelB : ", node.show())
                                    for rec in accept["pathlist"]:
                                        print(rec.show())
                                        for species in rec.getrea():
                                            allnodes.add(species)
                                        for species in rec.getpro():
                                            allnodes.add(species)
                                        allnodes.add(rec.getenz())

                                    print("/////////////////all related reaction///////////////")
                                    print(self.check(allnodes, accept["pathlist"]))
                                    for node in allnodes:
                                        print ("present(" + node.show() + ").")
                                    return True
                    self.ClearPath()
                    self.ClearLevel()

        print(len(result))
        return False

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

    def MergeBside(self, label0, label1):
        path_tmp = {}
        for sidepath in label0[0].sidepath:
            for path in label1[0].path:
                if self.CheckMerge(path, sidepath):
                    path0 = self.Merge(path, sidepath)
                    downRec  = self.mapToEdge[label1[1]]
                    pathnode = []
                    pathnode.append(label1[0])
                    related  = []
                    enz      = []
                    enz.append(downRec.getenz())
                    for rea in downRec.getrea():
                        if rea not in pathnode:
                            related.append(rea)
                    path1 = {}
                    path1["pathnode"] = pathnode
                    path1["related"]  = related
                    path1["pathlist"] = []
                    path1["pathlist"].append(downRec)
                    path1["pathenz"]  = enz
                    if self.CheckMerge(path0, path1):
                        path_tmp = self.Merge(path0, path1)
                        return path_tmp
        return {}

    def BuildStart(self, input_species):
        for node in input_species:
            tmp             = {}
            tmp["related"]  = set()
            tmp["pathlist"] = []
            tmp["pathenz"]  = []
            tmp["pathnode"] = []
            node.path.append(tmp)
        for node in input_species:
            node.path[0]["related"].add(self.mapToNode["H2O"])
            node.path[0]["pathenz"].append(self.mapToNode["spontaneous_reaction"])
            
        (input_species[0].path[0])["pathnode"].append(input_species[1])
        (input_species[1].path[0])["pathnode"].append(input_species[0])

    def CheckPathnode(self, path0, path1):
        assert (len(path0) > 2)
        for species in path1["pathnode"]:
            if (species in path0["pathnode"][2:]):
                return False
        return True

    def CheckMerge(self, path0, path1):
        for species in path1["pathnode"]:
            if (species in path0["pathnode"][2:]):
                return False
        for species in path1["pathnode"]:
            if (species in path0["related"]):
                return False
        for species in path0["pathnode"]:
            if species in path1["related"]:
                return False
        return True

    def check(self, allnodes, pathlist):
        temp = ""
        count = 1
        allnode = allnodes.copy()
        allrec  = set()
        while count != 0:
            tmp = len(allrec)
            for enz in allnodes:
                for reaction in enz.getCatedge():
                    if reaction not in pathlist and reaction not in allrec:
                        if reaction.activated(allnode):
                            for pro in reaction.getpro():
                                allnode.add(pro)
                            allrec.add(reaction)
                            temp += (reaction.show())
            count = len(allrec) - tmp
        return temp
    
    def ClearVis(self):
        for nodes in self.nodeList:
            nodes.visited  = 0
            nodes.recStack = False
        for rec in self.edgeList:
            rec.visited  = 0
            rec.recStack = False

    def ClearTmp(self):
        for nodes in self.nodeList:
            nodes.tmp_vis  = 0
            nodes.recStack = False
        for rec in self.edgeList:
            rec.tmp_vis  = 0
            rec.recStack = False

    
    def ClearPath(self):
        for node in self.nodeList:
            node.path       = []

    def ClearRecordA(self):
        for node in self.nodeList:
            node.recordA    = []

    def ClearLevel(self):
        for node in self.nodeList:
            node.level = 0
