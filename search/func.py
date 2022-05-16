from node import node

from edge import edge


def SearchUp(startNode, search_limit):
    startNode.layer = 1
    bfs          = [startNode]

    # start to traverse
    while (bfs != []):
        startNode = bfs.pop(0)

        if startNode.layer <= search_limit:
            for upRec in startNode.getUpedge():
                for rea in upRec.getrea():
                    if (rea.layer == 0):
                        rea.layer = startNode.layer + 1
                        bfs.append(rea)
 
def SearchDown(startNode, Target_type, up_limit, down_limit):
    search_limit = down_limit
    bfs          = [startNode]
    stype        = Target_type
    candidate    = set()

    # start to traverse
    while (bfs != []):
        startPro = bfs.pop(0)
        if (stype in startPro.label):
            found(stype, startPro, candidate)

        if startPro.level + 1 <= search_limit:
            # for path in startPro.path:
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                CheckDownRec = startPro.CheckDownRec(downRec)
                if CheckDownRec != []:
                    for product in downRec.getpro():
                        if product.layer != 0:
                            CheckProduct = product.CheckProduct(CheckDownRec, up_limit)
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

# def found(stype, product, candidate):
#     if stype == "A" and (product.getpin() > 2):
#         if product.path != []:
#             candidate.add(product)
#     elif stype == "C_side":
#         for label in product.labelC_side:
#             candidate.add(product)
#     elif stype == "B" and (product.getpout() > 1):
#         assert(len(product.getDownedge()) > 1)
#         for downRec in product.getDownedge():
#             candidate.add(product)
#             downRec.labelB.append(product.labelB)
#             for rea in downRec.getrea():
#                 if rea != product:
#                     rea.label.add("B_side")
#                     # rea.labelB_side.append(product.labelB)
#                     rea.B_sideRec = downRec
#     elif stype == "B_side":
#         for downRec in product.getDownedge():
#             if downRec.labelB!=[]:
#                 # for path in product.path:
#                 #     path["pathlist"].append(product.B_sideRec)
#                 #     path["pathnode"].append(product)
#                 candidate.add(product)
#     elif stype == "output":
#         candidate.add(product)

def found(stype, product, candidate):
    if stype == "A" and (product.getpin() > 2):
        if product.path != []:
            candidate.add((product, product.labelA))
    elif stype == "C_side":
        for label in product.labelC_side:
            candidate.add((product, label))
    elif stype == "B" and (product.getpout() > 1):
        assert(len(product.getDownedge()) > 1)
        for downRec in product.getDownedge():
            candidate.add((product, downRec.name))
            downRec.labelB.append(product.labelB)
            for rea in downRec.getrea():
                if rea != product:
                    rea.label.add("B_side")
                    # rea.labelB_side.append(product.labelB)
                    # rea.labelBmap[product.labelB] = downRec.name
    elif stype == "B_side":
        for downRec in product.getDownedge():
            if downRec.labelB!=[]:
                candidate.add((product, downRec.name))
    elif stype == "output":
        candidate.add((product, "output"))


class search:
    def __init__(self):
        self.c = 0
        self.up_limit       = 6
        self.input_species  = []
        self.output_species = []
        self.forbid_node    = set()
        self.mapToNode      = {}
        self.mapToEdge      = {}
        self.mapToClist     = {}
        self.mapToCnode     = {}
        self.mapToCforbid   = {}
        self.nodeList       = []
        self.edgeList       = []
        self.HashTable      = {}
        self.HashCollect    = {}
        self.reaction       = 0
        self.typeA          = 0
        self.typeB          = 0
        self.typeC          = 0
        self.typeD          = 0
        self.typeE          = 0

        self.count          = 0

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
                            if species[1] == "_":
                                if species[0].isdigit():
                                    par = int(species[0])
                                    species = species[2:]
                                elif species[0] == "a":
                                    species = species[2:]
                        elif len(species) > 3:
                            if species[2] == "_":
                                if species[0:2].isdigit():
                                    par = int(species[0:2])
                                    species = species[3:]
                                elif species[0:2] == "an":
                                    species = species[3:]

                        if species[0:6] == "NAD(P)":
                            species0 = "NAD" + species[6:]
                            if (species0 not in self.mapToNode):
                                self.mapToNode[species0] = node(species0)
                                self.nodeList.append(self.mapToNode[species0])
                            
                            if (is_product == 0):
                                self.mapToEdge[self.reaction].addrea(self.mapToNode[species0])
                                self.mapToNode[species0].addDownedge(self.mapToEdge[self.reaction])
                                (self.mapToNode[species0]).downEdgePar[self.mapToEdge[self.reaction]] = par
                                self.mapToNode[species0].addpout()
                            else:
                                self.mapToEdge[self.reaction].addpro(self.mapToNode[species0])
                                self.mapToNode[species0].addUpedge(self.mapToEdge[self.reaction])
                                (self.mapToNode[species0]).upEdgePar[self.mapToEdge[self.reaction]] = par
                                self.mapToNode[species0].addpin()
                            species1 = "NADP" + species[6:]
                            if (species1 not in self.mapToNode):
                                self.mapToNode[species1] = node(species1)
                                self.nodeList.append(self.mapToNode[species1])
                            
                            if (is_product == 0):
                                self.mapToEdge[self.reaction].addrea(self.mapToNode[species1])
                                self.mapToNode[species1].addDownedge(self.mapToEdge[self.reaction])
                                (self.mapToNode[species1]).downEdgePar[self.mapToEdge[self.reaction]] = par
                                self.mapToNode[species1].addpout()
                            else:
                                self.mapToEdge[self.reaction].addpro(self.mapToNode[species1])
                                self.mapToNode[species1].addUpedge(self.mapToEdge[self.reaction])
                                (self.mapToNode[species1]).upEdgePar[self.mapToEdge[self.reaction]] = par
                                self.mapToNode[species1].addpin()

                        else:  
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
        # print("in record label parsing")
        # tmpList = []
        # typeC = 0
        # forbid = False
        # for line in inputfile:
        #     tmpList = self.split(line)
        #     if len(tmpList) > 0:
        #         if tmpList[0] == "TypeC:":
        #             self.typeC += 1
        #             forbid = False
        #             typeC = tmpList[1]
        #             self.mapToClist[typeC] = []
        #             assert(tmpList[1].isdigit())
        #         elif tmpList[0] == "reaction:":
        #             self.mapToClist[typeC].append(self.mapToEdge[tmpList[1]])
        #             self.mapToEdge[tmpList[1]].label.add("C")
        #             self.mapToEdge[tmpList[1]].labelC.append(typeC)
        #         elif tmpList[0] == "CoreNode:":
        #             assert len(tmpList) == 5
        #             self.mapToCnode[typeC] = []
        #             self.mapToCnode[typeC].append(self.mapToNode[tmpList[1]])
        #             self.mapToCnode[typeC].append(self.mapToNode[tmpList[2]])
        #             self.mapToCnode[typeC].append(self.mapToNode[tmpList[3]])
        #             self.mapToCnode[typeC].append(self.mapToNode[tmpList[4]])
        #             self.mapToNode[tmpList[1]].label.add("C")
        #             self.mapToNode[tmpList[2]].label.add("C")
        #             self.mapToNode[tmpList[3]].label.add("C_side")
        #             self.mapToNode[tmpList[4]].label.add("C_side")
        #             self.mapToNode[tmpList[1]].labelC.append(typeC)
        #             self.mapToNode[tmpList[2]].labelC.append(typeC)
        #             self.mapToNode[tmpList[3]].labelC_side.append(typeC)
        #             self.mapToNode[tmpList[4]].labelC_side.append(typeC)
                    
        #         elif tmpList[0] == "Forbidnode:":
        #             forbid = True
        #             self.mapToCforbid[typeC] = []
        #         else:
        #             if forbid:
        #                 for i in tmpList:
        #                     tmp = []
        #                     tmp.append(self.mapToNode[i])
        #                 self.mapToCforbid[typeC].append(tmp)

        print("start initializing")
        self.initialize()

    def test(self):
        count = 0
        for node in self.nodeList:
            for downRec in node.getDownedge():
                count = count + len(downRec.getpro())
        print ("average next node count : ", count/len(self.nodeList))

    def main(self, input_species, output_species):
        self.input_species  = input_species
        self.output_species = output_species

        for node in self.output_species:
            SearchUp(node, self.up_limit)
            node.label.add("output")
        for node in self.input_species:
            print (node.name , node.layer)

        candidate = []
        self.BuildStart(input_species)
        input_species[0].CopyToPath()
        candidate0 = SearchDown(input_species[0], "A", self.up_limit, 1)
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord()
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0))

        self.BuildStart(input_species)
        input_species[1].CopyToPath()
        candidate1 = SearchDown(input_species[1], "A", self.up_limit, 1)
        # check if there is any common label found between input0 and input1 path
        print("start merging type A")
        for label in candidate0:
            if label in candidate1:
                if label[0].MergeAll("A"):
                    candidate.append(label[0])
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate1))
        print(len(candidate))

        for i in range(len(candidate)):
            # count = 0
            # for downRec in candidate[i].getDownedge():
            #     count = count + len(downRec.getpro())
            # print ("next node count : ", count)
            candidate[i].CopyToPath()
            candidate0 = SearchDown(candidate[i], "B", self.up_limit, 1)
            print ("candidate0 : ", len(candidate0))
            print ("candidate name : ", candidate[i].name)
            for species in candidate0:
                species[0].CopyToSide()
            self.ClearPath()
            self.ClearLevel()
            c = 0
            for j in range(len(candidate)):
                if i != j:
                    candidate[j].CopyToPath()
                    candidate1 = SearchDown(candidate[j], "B_side", self.up_limit, 1)
                    print ("candidate1 : ", len(candidate1))
                    print ("candidate name : ", candidate[j].name)
                    print ("start merging type B_side")

                    for label0 in candidate0:
                        # if label0 in self.output_species:
                        if label0[0].layer <= 2:
                        # print ("layer <= 2")
                            for label1 in candidate1:
                                if (label0[1] == label1[1]) and (label0[0] != label1[0]):
                                    # print ("label")
                                    if label0[0].MergeBside(label1[0]):
                                        # print ("find result : ", label0[0].layer)
                                        result = SearchDown(label1[0], "output", self.up_limit, 3)
                                        if len(result) != 0:
                                            # print ("result :", len(result))
                                            for x in result:
                                                for path in x[0].path:
                                                    allnodes = set()
                                                    for node in path["pathnode"]:
                                                        print(node.show())
                                                    # for node in labelB:
                                                    #     print("labelB : ", node.show())
                                                    for rec in path["pathlist"]:
                                                        print(rec.show(self.c))
                                                        self.c = self.c+1
                                                        for species in rec.getrea():
                                                            allnodes.add(species)
                                                        for species in rec.getpro():
                                                            allnodes.add(species)
                                                        allnodes.add(rec.getenz()) 

                                                    print("/////////////////all related reaction///////////////")
                                                    print(self.Check(allnodes, path["pathlist"]))
                                                    for node in allnodes:
                                                        print ("present(" + node.show() + ").")
                                                    return True
                    self.ClearPath()
                    self.ClearLevel()
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
            for i in sidepath["pathlist"][-1].getpro():
                if i != label0[0]:
                    sidepath["related"].add(i)
        for path in label1[0].path:
            for i in path["pathlist"][-1].getpro():
                    if i != label1[0]:
                        path["related"].add(i)

        for sidepath in label0[0].sidepath:
            for path in label1[0].path:
                # flag = False
                for downRec0 in label0[0].getDownedge():
                    if downRec0 not in sidepath["pathlist"] and downRec0 not in path["pathlist"] and downRec0.getenz().name != "spontaneous_reaction":
                        sidepath["pathlist"].append(downRec0)
                        break

                if self.CheckMerge(path, sidepath): 
                    label0[0].mark = 1
                    label1[0].mark = 2
                    for rec in path["pathlist"]:
                        if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
                            rec.getenz().mark = 1
                    for rec in sidepath["pathlist"]:
                        if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
                            rec.getenz().mark = 2
                    for node in path["pathnode"]:
                        if (node not in self.input_species):
                            node.mark = 1
                    for node in sidepath["pathnode"]:
                        if (node not in self.input_species):
                            node.mark = 2
                    path0 = self.Merge(path, sidepath)
                    downRec  = self.mapToEdge[label1[1]]
                    pathnode = []
                    pathnode.append(label1[0])
                    pathnode.append(label0[0])
                    path1 = self.CreatePath(pathnode, downRec)

                    if self.CheckMerge(path0, path1):
                        path_tmp = self.Merge(path0, path1)
                        if self.CheckAll(path_tmp):
                            return path_tmp
                    label0[0].mark = 0
                    label1[0].mark = 0
                    for rec in path["pathlist"]:
                        rec.getenz().mark = 0
                    for rec in sidepath["pathlist"]:
                        rec.getenz().mark = 0
                    for node in path["pathnode"]:
                        node.mark = 0
                    for node in sidepath["pathnode"]:
                        node.mark = 0
        return {}

    def CreatePath(self, pathnode, downRec):
        related  = set()
        enz      = []
        enz.append(downRec.getenz())
        for rea in downRec.getrea():
            if rea not in pathnode:
                related.add(rea)
        for pro in downRec.getpro():
            related.add(pro)
        path = {}
        path["pathnode"] = pathnode
        path["related"]  = related
        path["pathlist"] = []
        path["pathlist"].append(downRec)
        path["pathenz"]  = enz

        return path

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

    def CheckAll(self, path):
        # for node in path["pathnode"]:
        #     if node.name == "NAD(P)H" or node.name == "NAD(P)+":
        #         return False
        allnodes = set()
        allnodes.add(self.mapToNode["H2O"])
        allnodes.add(self.mapToNode["spontaneous_reaction"])
        for rec in path["pathlist"]:
            for node in rec.getrea():
                allnodes.add(node)
            for node in rec.getpro():
                allnodes.add(node)
            allnodes.add(rec.getenz())

        allrec = self.CollectAll(allnodes, path["pathlist"])

        for rec in allrec:
            cr1 = 0
            cr2 = 0
            cp1 = 0
            cp2 = 0
            enzMark = rec.getenz().mark

            for pro in rec.getpro():
                if pro.mark == 1:
                    cp1 = cp1 + 1
                elif pro.mark == 2:
                    cp2 = cp2 + 1
            for rea in rec.getrea():
                if rea in self.input_species:
                    if cp1 != 0 and cp2 != 0:
                        return False
                if rea.mark == 1:
                    cr1 = cr1 + 1
                elif rea.mark == 2:
                    cr2 = cr2 + 1
            if cr1+cr2 == 0:
                if enzMark == 0:
                    return False
            if cr1 != 0 and cp2 != 0:
                return False
            if cp1 != 0 and cr2 != 0:
                return False

        return True

    def CheckPathnode(self, path0, path1):
        # assert (len(path0) > 2)
        for species in path1["pathnode"]:
            if (species not in self.input_species) and (species in path0["pathnode"]):
                return False
        return True

    def CheckMerge(self, path0, path1):
        for species in path1["pathnode"]:
            if (species not in self.input_species) and (species in path0["pathnode"]):
                return False
        for species in path1["pathnode"]:
            if species in path0["related"]:
                return False
        for species in path0["pathnode"]:
            if species in path1["related"]:
                return False
        return True

    def Check(self, allnodes, pathlist):
        temp = ""
        count = 1
        allnode = allnodes.copy()
        allrec  = set()
        allnode.add(self.mapToNode["H2O"])
        allnode.add(self.mapToNode["spontaneous_reaction"])
        while count != 0:
            tmp = len(allrec)
            for enz in allnodes:
                for reaction in enz.getCatedge():
                    if reaction not in pathlist and reaction not in allrec:
                        if reaction.activated(allnode):
                            for pro in reaction.getpro():
                                allnode.add(pro)
                            allrec.add(reaction)
                            temp += (reaction.show(self.c))
                            self.c = self.c+1
            count = len(allrec) - tmp
        return temp
    
    def CollectAll(self, allnodes, pathlist):
        temp = ""
        count = 1
        allnode = allnodes.copy()
        allnode.add(self.mapToNode["H2O"])
        allnode.add(self.mapToNode["spontaneous_reaction"])
        allrec  = set()
        while count != 0:
            tmp = len(allrec)
            for enz in allnodes:
                for reaction in enz.getCatedge():
                    if (reaction not in pathlist) and (reaction not in allrec) and (reaction.activated(allnode)):
                        for pro in reaction.getpro():
                            allnode.add(pro)
                        allrec.add(reaction)
            count = len(allrec) - tmp
        return allrec

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
