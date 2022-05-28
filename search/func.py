from node import node

from edge import edge


def aa(path, rec):
    for rea in rec.getrea():
        if rea in path["pathnode"]:
            return False
    return True

def CheckRec(upRec, path):
    path_tmp = []
    for path in path:
        if aa(path, upRec):
            path["pathlist"].append(upRec)
            path_tmp.append(path)
    return path_tmp


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
 
def SearchDown(startNode, Target_type, up_limit, down_limit, outNode):
    search_limit = down_limit
    bfs          = [startNode]
    candidate    = set()

    # startNode.CopyToPath()
    # start to traverse
    while (bfs != []):
        startPro = bfs.pop(0)
        if (Target_type in startPro.label):
            found(Target_type, startPro, candidate)
            # if Target_type == "output":
            #     return candidate

        if startPro.level + 1 <= search_limit:
            # for path in startPro.path:
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                CheckDownRec = startPro.CheckDownRec(downRec)
                if CheckDownRec != []:
                    # flag = True
                    # for node in outNode:
                    if Target_type == "output" and outNode in downRec.getpro():
                        CheckProduct = outNode.CheckProduct(CheckDownRec, up_limit, outNode)
                        if CheckProduct != []:
                            outNode.AddPath(startPro, downRec, CheckProduct)
                            if (outNode not in bfs):
                                bfs.append(outNode)

                    else:
                        for product in downRec.getpro():
                            if product.layer != 0:
                                CheckProduct = product.CheckProduct(CheckDownRec, up_limit, outNode)
                                if CheckProduct != []:
                                    product.AddPath(startPro, downRec, CheckProduct)
                                    product.level = startPro.level + 1
                                    if (product not in bfs and product.level <= search_limit):
                                        bfs.append(product)
                                    assert product.path != CheckProduct
                                    assert len(product.path) != 0
    return candidate

def found(stype, product, candidate):
    if stype == "A" and (product.getpin() > 2):
        if product.path != []:
            candidate.add((product, product.labelA))
    elif stype == "C_side":
        for label in product.labelC_side:
            candidate.add((product, label))
    elif stype == "B" and (product.getpout() > 1):
        assert(len(product.getDownedge()) > 1)
        # downRec.labelB.append(product.labelB).
        MaxPathRec = product.getDownedge()[0]
        MaxCount = 0
        for downRec in product.getDownedge():
            CheckDownRec = product.CheckDownRec(downRec)
            if len(CheckDownRec) > MaxCount:
                MaxCount = len(CheckDownRec)
                MaxPathRec = downRec

        candidate.add((product, MaxPathRec.name))
        if MaxPathRec.pathA1 == "" :
            MaxPathRec.pathA1 = product
        else:
            if product.layer != 0 and product.layer < MaxPathRec.pathA1.layer:
                MaxPathRec.pathA1 = product

        for rea in MaxPathRec.getrea():
            if rea != product:
                rea.label.add("B_side")
                rea.labelB_side.append(MaxPathRec)

    elif stype == "B_side":
        for rec in product.labelB_side:
            CheckDownRec = product.CheckDownRec(rec)
            if CheckDownRec != []:
                for path in CheckDownRec:
                    path["pathnode"].append(product)
                product.path = CheckDownRec
                candidate.add((product, rec.name))
                break
                # total_path = []
                # for upRec in product.getUpedge():
                #     CheckUpRec = CheckRec(upRec, CheckDownRec.copy())
                #     if CheckUpRec != []:
                #         # total_path += CheckUpRec.copy()
                #         # candidate.add((product, rec.name))
                #         for path in CheckUpRec:
                #             path["pathnode"].append(product)
                #         product.path = CheckUpRec
                #         candidate.add((product, rec.name))
                #         break
                # if total_path != []:
                #     for path in total_path:
                #         path["pathnode"].append(product)
                #     product.path = total_path
                #     candidate.add((product, rec.name))
                #     break

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

        # for node in self.output_species:
        # for node in self.output_species:
        self.output_species.label.add("output")
        SearchUp(self.output_species, self.up_limit)
        for node in self.input_species:
            print (node.name , node.layer)

        candidate = []
        self.BuildStart(input_species)
        input_species[0].CopyToPath("A")
        candidate0 = SearchDown(input_species[0], "A", self.up_limit, 1, self.output_species)
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord("A")
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0)) 

        self.BuildStart(input_species)
        input_species[1].CopyToPath("A")
        candidate1 = SearchDown(input_species[1], "A", self.up_limit, 1, self.output_species)
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
            candidate[i].CopyToPath("A")
            candidate0 = SearchDown(candidate[i], "B", self.up_limit, 1, self.output_species)
            print ("candidate0 : ", len(candidate0))
            print ("candidate name : ", candidate[i].name)
            for species in candidate0:
                species[0].CopyToRecord("A1")
            self.ClearPath()
            self.ClearLevel()

            candidate_B = []
            for j in range(len(candidate)):
                if i != j:
                    candidate[j].CopyToPath("A")
                    candidate1 = SearchDown(candidate[j], "B_side", self.up_limit, 1, self.output_species)
                    print ("start merging B")
                    for species in candidate1:
                        targetRec = self.mapToEdge[species[1]]
                        targetB = targetRec.pathA1
                        targetB.path = species[0].path
                        # for upRec in species[0].getUpedge():
                            

                        merge = self.MergeBside(targetB.path, targetB.pathA1)
                        if merge != []:
                            print("merge success")
                            targetB.MergeAll("B_side")
                            candidate_B.append((targetB, targetRec, species[0]))
                    print (len(candidate_B))
                        # species[0].CopyToRecord("A2")
                    self.ClearPath()
                    self.ClearLevel()

                    for target in candidate_B:
                        target[0].CopyToPath("A2")
                        # print (target.show())
                        output = SearchDown(target[0], "output", self.up_limit, 3, self.output_species)

                        for can in output:
                            assert (can[0] == self.output_species)
                        print ("target : ", target[0].show(), target[0].layer)
                        c = 0
                        if len(output) != 0:
                            print (target[1].show(c))
                            c = c+1
                            for path in self.output_species.path:
                                # for upRec in target[2]:
                                    
                                for rec in path["pathlist"]:
                                    print (rec.show(c))
                                    c = c+1 
                                for node in path["pathnode"]:
                                    print (node.show())
                                print ("******************other reactions****************")
                                allnodes = set()
                                for rec in path["pathlist"]:
                                    for node in rec.getrea():
                                        allnodes.add(node)
                                    for node in rec.getpro():
                                        allnodes.add(node)
                                    allnodes.add(rec.getenz())
                                allrec = self.CollectAll(allnodes, path["pathlist"])
                                for rec in allrec:
                                    print (rec.show(c))
                                    c = c+1
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

    def Check2Path(self, path, sidepath):
        path_tmp = {}
        if self.CheckMerge(path, sidepath):
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
            path_tmp = self.Merge(path, sidepath)
            if self.CheckAll(path_tmp):
                return path_tmp
            for rec in path["pathlist"]:
                rec.getenz().mark = 0
            for rec in sidepath["pathlist"]:
                rec.getenz().mark = 0
            for node in path["pathnode"]:
                node.mark = 0
            for node in sidepath["pathnode"]:
                node.mark = 0
        return {}

    def MergeBside(self, label0, label1):
        path_tmp = []
        for sidepath in label0:
            for i in sidepath["pathlist"][-1].getpro():
                sidepath["related"].add(i)
        for path in label1:
            for i in path["pathlist"][-1].getpro():
                path["related"].add(i)

        for sidepath in label0:
            for path in label1:
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
                # path0 = self.Merge(path, sidepath)
                # downRec  = self.mapToEdge[label1[1]]
                # pathnode = []
                # path1 = self.CreatePath(pathnode, downRec)

                if self.CheckMerge(path, sidepath):
                    tmp = self.Merge(path, sidepath)
                    if self.CheckAll(tmp):
                        path_tmp.append(tmp)
                for rec in path["pathlist"]:
                    rec.getenz().mark = 0
                for rec in sidepath["pathlist"]:
                    rec.getenz().mark = 0
                for node in path["pathnode"]:
                    node.mark = 0
                for node in sidepath["pathnode"]:
                    node.mark = 0
        return path_tmp

    # def MergeBside(self, label0, label1):
    #     path_tmp = []
    #     for sidepath in label0[0].sidepath:
    #         for i in sidepath["pathlist"][-1].getpro():
    #             if i != label0[0]:
    #                 sidepath["related"].add(i)
    #     for path in label1[0].path:
    #         for i in path["pathlist"][-1].getpro():
    #                 if i != label1[0]:
    #                     path["related"].add(i)

    #     for sidepath in label0[0].sidepath:
    #         for path in label1[0].path:
    #             # flag = False
    #             for downRec0 in label0[0].getDownedge():
    #                 if downRec0 not in sidepath["pathlist"] and downRec0 not in path["pathlist"] and downRec0.getenz().name != "spontaneous_reaction":
    #                     sidepath["pathlist"].append(downRec0)
    #                     break

    #             if self.CheckMerge(path, sidepath): 
    #                 label0[0].mark = 1
    #                 label1[0].mark = 2
    #                 for rec in path["pathlist"]:
    #                     if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
    #                         rec.getenz().mark = 1
    #                 for rec in sidepath["pathlist"]:
    #                     if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
    #                         rec.getenz().mark = 2
    #                 for node in path["pathnode"]:
    #                     if (node not in self.input_species):
    #                         node.mark = 1
    #                 for node in sidepath["pathnode"]:
    #                     if (node not in self.input_species):
    #                         node.mark = 2
    #                 path0 = self.Merge(path, sidepath)
    #                 downRec  = self.mapToEdge[label1[1]]
    #                 pathnode = []
    #                 pathnode.append(label1[0])
    #                 pathnode.append(label0[0])
    #                 path1 = self.CreatePath(pathnode, downRec)

    #                 if self.CheckMerge(path0, path1):
    #                     tmp = self.Merge(path0, path1)
    #                     if self.CheckAll(tmp):
    #                         path_tmp.append(tmp)
    #                 label0[0].mark = 0
    #                 label1[0].mark = 0
    #                 for rec in path["pathlist"]:
    #                     rec.getenz().mark = 0
    #                 for rec in sidepath["pathlist"]:
    #                     rec.getenz().mark = 0
    #                 for node in path["pathnode"]:
    #                     node.mark = 0
    #                 for node in sidepath["pathnode"]:
    #                     node.mark = 0
    #     return path_tmp

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
