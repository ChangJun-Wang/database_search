from node import node

from edge import edge

import copy

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
 
def SearchDown(startNode, Target_type, up_limit, down_limit, inNodes, outNode):
    search_limit = down_limit
    bfs          = [startNode]
    candidate    = set()

    # startNode.CopyToPath()
    # start to traverse
    while (bfs != []):
        startPro = bfs.pop(0)
        if (Target_type in startPro.label):
            found(Target_type, startPro, candidate, inNodes)

        if startPro.level+1 <= search_limit:
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                CheckDownRec = startPro.CheckDownRec(downRec)
                if CheckDownRec != []:
                    for product in downRec.getpro():
                        if product.layer != 0 and product not in downRec.getrea():
                            CheckProduct = product.CheckProduct(CheckDownRec, up_limit, outNode)
                            if CheckProduct != []:
                                # if Target_type == "output":
                                # product.SetPath(startPro, downRec, CheckProduct)
                                # else:
                                product.AddPath(startPro, downRec, CheckProduct)
                                product.level = startPro.level + 1
                                if (product not in bfs and product.level <= search_limit):
                                    bfs.append(product)
                                # assert product.path != CheckProduct
                                # assert len(product.path) != 0
    return candidate

def found(stype, product, candidate, inNodes):
    if stype == "A" and (product.getpin() > 2):
        if product.path != []:
            candidate.add((product, product.labelA))
    elif stype == "C_side":
        for label in product.labelC_side:
            candidate.add((product, label))
    elif stype == "B" and (product.getpout() > 1):
        MaxPathRec = product.getDownedge()[0]
        MaxCount = 0
        for downRec in product.getDownedge():
            CheckDownRec = product.CheckDownRec(downRec)
            if len(CheckDownRec) > MaxCount:
                MaxCount = len(CheckDownRec)
                MaxPathRec = downRec

        candidate.add((product, MaxPathRec.name))
        if MaxPathRec.rec_node == None :
            MaxPathRec.rec_node = product
        else:
            if product.layer != 0 and product.layer < MaxPathRec.rec_node.layer:
                MaxPathRec.rec_node = product

        for rea in MaxPathRec.getrea():
            if rea != product:
                rea.label.add("B_side")
                rea.labelB_side.append(MaxPathRec)

    elif stype == "B_side":
        for rec in product.labelB_side:
            if ~rec.is_reverse():
                CheckDownRec = product.CheckDownRec(rec)
                if CheckDownRec != []:
                    candidate.add((product, rec.name))
                    break

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
        self.typeF          = 0

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
            if (len(nodes.getCatedge()) > 1) and (nodes.getpin() > 1 or nodes.getpout() > 1):
                self.typeF += 1

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
                self.reaction = self.reaction + 1
                self.mapToEdge[self.reaction] = edge(self.reaction)
                self.edgeList.append(self.mapToEdge[self.reaction])

                self.mapToEdge[self.reaction].setenz(self.mapToNode[enzyme])
                self.mapToNode[enzyme].addCatedge(self.mapToEdge[self.reaction])

                for species in tmpList[1:]:
                    HashKey = HashKey + species[0]
                    if (species == "="):
                        is_product = 1
                    else:
                        if (species == "NO"):
                            species = "nitric_oxide"
                        if (species not in self.mapToNode):
                            self.mapToNode[species] = node(species)
                            self.nodeList.append(self.mapToNode[species])
                        if (is_product == 0):
                            self.mapToEdge[self.reaction].addrea(self.mapToNode[species])
                            self.mapToNode[species].addDownedge(self.mapToEdge[self.reaction])
                            self.mapToNode[species].addpout()
                        else:
                            self.mapToEdge[self.reaction].addpro(self.mapToNode[species])
                            self.mapToNode[species].addUpedge(self.mapToEdge[self.reaction])
                            self.mapToNode[species].addpin()

    def parsing_label(self, inputfile):
        # print("in record label parsing")
        # tmpList = []
        # typeC = 0
        # forbid = False
        for line in inputfile:
            tmpList = self.split(line)
            if len(tmpList) > 0:
                if tmpList[0] == "TypeC:":
                    self.typeC += 1
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
        # self.output_species = output_species

        # for node in self.output_species:
        for node in output_species:
            if node in self.mapToNode:
                self.output_species.append(self.mapToNode[node])
        for node in self.output_species:
            node.label.add("output")
            SearchUp(node, self.up_limit)
        for node in self.input_species:
            print (node.name , node.layer)

        candidate = []
        self.BuildStart(input_species)
        # input_species[0].CopyToPath("A")
        candidate0 = SearchDown(input_species[0], "A", self.up_limit, 1, self.input_species, self.output_species)
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord("A")
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0))

        # print ("///////////////////GMP////////////////////")
        # for record in self.mapToNode["GMP"].recordA:
        #     for rec in record["pathlist"]:
        #         print (rec.show(0))
        #     print ("***************record******************")
        # print ("///////////////////AMP////////////////////")
        # for record in self.mapToNode["AMP"].recordA:
        #     for rec in record["pathlist"]:
        #         print (rec.show(0))
        #     print ("***************record******************")

        # for inNode in input_species:
        self.BuildStart(input_species)
        # input_species[1].CopyToPath("A")
        candidate1 = SearchDown(input_species[1], "A", self.up_limit, 1, self.input_species, self.output_species)
        # check if there is any common label found between input0 and input1 path
        print("start merging type A")
        for label in candidate0:
            if label in candidate1:
                label[0].MergeAll("A")
                candidate.append(label[0])
            else:
                label[0].ClearRecordA()
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate1))
        print(len(candidate))

        # print ("///////////////////GMP////////////////////")
        # for record in self.mapToNode["GMP"].recordA:
        #     for rec in record["pathlist"]:
        #         print (rec.show(0))
        #     print ("***************record******************")
        # print ("///////////////////AMP////////////////////")
        # for record in self.mapToNode["AMP"].recordA:
        #     for rec in record["pathlist"]:
        #         print (rec.show(0))
        #     print ("***************record******************")

        ans = 100
        record_ans = []
        record_oth = []
        record_node = []
        # cand = [self.mapToNode["NADH"], self.mapToNode["GMP"]]
                
        for i in range(len(candidate)):
            self.ClearPathA1()
            self.ClearPathA2()
            self.ClearPath()
            self.ClearLevel()
            print ("candidate name : ", candidate[i].name)
            # if self.input_species == [self.mapToNode["acetone"], self.mapToNode["D-glucose"]]:
            #     if candidate[i].name != "GMP":
            #         continue
            candidate[i].CopyToPath("A")
            candidateB = SearchDown(candidate[i], "B", self.up_limit, 1, self.input_species, self.output_species)
            print ("candidateB : ", len(candidateB))
            # print ("candidate name : ", candidate[i].name)
            for species in candidateB:
                species[0].CopyToRecord("A1")
            
            for j in range(len(candidate)):
                if candidate[i] != candidate[j]:
                    self.ClearPath()
                    self.ClearLevel()
                    candidate_B = []
                    candidate[j].CopyToPath("A")
                    candidate_Bside = SearchDown(candidate[j], "B_side", self.up_limit, 1, self.input_species, self.output_species)
                    print ("start merging B")
                    for species in candidate_Bside:
                        targetRec = self.mapToEdge[species[1]]
                        targetB = targetRec.rec_node
                        targetB.path = species[0].path
                        # output = SearchDown(targetB, "output", self.up_limit, 2, self.input_species, self.output_species)
                        merge = self.MergeBside(targetB.path, targetB.pathA1, targetB, species[0])
                        if merge != []:
                            print("merge success")
                            targetB.pathA2 = merge
                            for path in targetB.pathA2:
                                path["pathlist"].append(targetRec)
                                for pro in targetRec.getpro():
                                    path["related"].add(pro)
                                for rea in targetRec.getrea():
                                    if rea not in path["pathnode"]:
                                        path["related"].add(rea)
                            candidate_B.append((targetB, species[0], targetRec))
                    print (len(candidate_B))
                        # species[0].CopyToRecord("A2")

                    for target in candidate_B:
                        self.ClearPath()
                        self.ClearLevel()
                        target[0].CopyToPath("A2")
                        target[1].mark = 1
                        target[2].getenz().mark = 1
                        # print (target.show())
                        output = SearchDown(target[0], "output", self.up_limit, 3, self.input_species, self.output_species)
                        target[1].mark = 0
                        target[2].getenz().mark = 0
                        # for can in output:
                        #     assert (can[0] == self.output_species)
                        print ("target : ", target[0].show(), target[0].layer)
                        # tmp = []
                        for outNode in self.output_species:
                            tmp = []
                            if len(outNode.path) == 0:
                                continue
                            for path in outNode.path:
                                allnodes = set()
                                for node in self.input_species:
                                    allnodes.add(node)
                                for rec in path["pathlist"]:
                                    for rea in rec.getrea():
                                        if rea not in path["pathnode"]:
                                            allnodes.add(rea)
                                for node in path["related"]:
                                    allnodes.add(node)
                                for rec in path["pathlist"]:
                                    allnodes.add(rec.getenz())
                                flag = True
                                for rec in outNode.getUpedge():
                                    if rec.activated(allnodes):
                                        flag = False
                                for node in path["pathnode"]:
                                    allnodes.add(node)
                                allnodes.add(outNode)
                                for rec in outNode.getDownedge():
                                    if rec.activated(allnodes):
                                        flag = False
                                if flag:
                                    tmp.append(path)
                            outNode.path = tmp
                            if len(outNode.path) != 0:
                                print ("target : ", target[0].show(), target[0].layer)
                                for path in outNode.path:
                                    thres = []
                                    for upRec in target[1].getUpedge():
                                        # thres.append(upRec)
                                        flag = True
                                        for rea in upRec.getrea():
                                            if rea in path["pathnode"]:
                                                flag = False
                                        for pro in upRec.getpro():
                                            if pro != target[1] and pro in path["pathnode"]:
                                                flag = False
                                            if pro == outNode:
                                                flag = False
                                        if flag and upRec.getenz().name != "spontaneous_reaction":
                                            thres.append(upRec)
                                    # ans = 100
                                    for x in thres:
                                        allnodes = set()
                                        for rec in path["pathlist"]:
                                            for node in rec.getrea():
                                                allnodes.add(node)
                                            for node in rec.getpro():
                                                allnodes.add(node)
                                            allnodes.add(rec.getenz())
                                        for node in x.getrea():
                                            allnodes.add(node)
                                        for node in x.getpro():
                                            allnodes.add(node)
                                        allnodes.add(x.getenz())
                                        allrec = self.CollectAll(allnodes, path["pathlist"])
                                        c = 2
                                        for rec in path["pathlist"]:
                                            c = c+1 
                                        for rec in allrec:
                                            c = c+1
                                        # return True
                                        # print ("c : ", c)
                                        self.count += 1
                                        if c < ans:
                                            ans = c
                                            record_node = path["pathnode"].copy()
                                            record_ans = path["pathlist"].copy()
                                            record_ans.append(x)
                                            record_oth = allrec.copy()

        print ("ans : ", ans, len(record_ans), len(record_oth))
        c = 0
        collect_rea = set()
        for rec in record_ans:
            for rea in rec.getrea():
                collect_rea.add(rea)
            print (rec.show(c))
            c = c+1
        print ("******************pathnode****************")
        for node in record_node:
            print( node.show() )
        for node in collect_rea:
            print("present(", node.show(), ", 1).")
        print ("******************other reactions****************")
        for rec in record_oth:
            print (rec.show(c))
            c = c+1

        print ("total solution number : ", self.count)
        # return False

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

    # def Check2Path(self, path, sidepath):
    #     path_tmp = {}
    #     if self.CheckMerge(path, sidepath):
    #         for rec in path["pathlist"]:
    #             if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
    #                 rec.getenz().mark = 1
    #         for rec in sidepath["pathlist"]:
    #             if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
    #                 rec.getenz().mark = 2
    #         for node in path["pathnode"]:
    #             if (node not in self.input_species):
    #                 node.mark = 1
    #         for node in sidepath["pathnode"]:
    #             if (node not in self.input_species):
    #                 node.mark = 2
    #         path_tmp = self.Merge(path, sidepath)
    #         if self.CheckAll(path_tmp):
    #             return path_tmp
    #         for rec in path["pathlist"]:
    #             rec.getenz().mark = 0
    #         for rec in sidepath["pathlist"]:
    #             rec.getenz().mark = 0
    #         for node in path["pathnode"]:
    #             node.mark = 0
    #         for node in sidepath["pathnode"]:
    #             node.mark = 0
    #     return {}

    def MergeBside(self, label0, label1, node0, node1):
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
                node0.mark = 3
                node1.mark = 4
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
                node0.mark = 0
                node1.mark = 0
        return path_tmp

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
            cr3 = 0
            cr4 = 0
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
                elif rea.mark == 3:
                    cr3 = 1
                elif rea.mark == 4:
                    cr4 = 1
            if cr1+cr2 == 0:
                if enzMark == 0:
                    return False
            if cr1 != 0 and cp2 != 0:
                return False
            if cp1 != 0 and cr2 != 0:
                return False
            if cr4 == 1 and cr3 == 0:
                return False
        return True

    def CheckPathnode(self, path0, path1):
        # assert (len(path0) > 2)
        for species in path1["pathnode"]:
            if (species not in self.input_species) and (species in path0["pathnode"]):
                return False
        return True

    def CheckMerge(self, path0, path1):
        if (path0["pathlist"][0].getenz() == path1["pathlist"][0].getenz()):
            return False
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

    def ClearPath(self):
        for node in self.nodeList:
            node.path       = []

    def ClearPathA1(self):
        for node in self.nodeList:
            node.pathA1       = []

    def ClearPathA2(self):
        for node in self.nodeList:
            node.pathA2       = []
            if "B_side" in node.label:
                node.label.remove("B_side")
            node.labelB_side = []
        for edge in self.edgeList:
            edge.rec_node   = None

    def ClearRecordA(self):
        for node in self.nodeList:
            node.recordA    = []

    def ClearLevel(self):
        for node in self.nodeList:
            node.level = 0
