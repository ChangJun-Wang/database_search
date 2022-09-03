from node import node

from edge import edge

import copy

import time

def SearchUp(startNode, search_limit):
    startNode.layer = 1
    bfs          = [startNode]

    while (bfs != []):
        startNode = bfs.pop(0)

        if startNode.layer <= search_limit:
            for upRec in startNode.getUpedge():
                for rea in upRec.getrea():
                    if (rea.layer == 0):
                        rea.layer = startNode.layer + 1
                        bfs.append(rea)
 
def SearchDown(startNode, Target_type, struc, down_limit, inNodes, outNode, out, ans):
    search_limit = down_limit
    bfs          = [startNode]
    candidate    = set()

    while (bfs != []):
        startPro = bfs.pop(0)
        if (Target_type in startPro.label):
            found(Target_type, startPro, candidate, inNodes, struc, out, ans)

        if startPro.level+1 <= search_limit:
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                CheckDownRec = startPro.CheckDownRec(downRec, out, ans)
                if CheckDownRec != []:
                    for product in downRec.getpro():
                        if product not in downRec.getrea(): #product.layer != 0 and 
                            CheckProduct = product.CheckProduct(CheckDownRec, 0, outNode)
                            if CheckProduct != []:
                                product.AddPath(startPro, downRec, CheckProduct)
                                product.level = startPro.level + 1
                                if (product not in bfs and product.level <= search_limit):
                                    bfs.append(product)
    return candidate

def found(stype, product, candidate, inNodes, struc, out, ans):
    if stype == "A" and (product.getpin() > 2):
        if product.path != []:
            candidate.add((product, product.labelA))
    elif stype == "B" and (product.getpout() > 1):
        MaxPathRec = product.getDownedge()[0]
        MaxCount = 0
        for downRec in product.getDownedge():
            CheckDownRec = product.CheckDownRec(downRec, out, ans)
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
                CheckDownRec = product.CheckDownRec(rec, out, ans)
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
        self.time           = 0
        self.bound          = 0
        self.out            = False

        self.ans = 100
        self.distribute = {}
        self.record_ans = []
        self.record_oth = []
        self.record_node = []
        self.record_struc = ()

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
            nodes.H2O = self.mapToNode["H2O"]
            nodes.spontaneous_reaction = self.mapToNode["spontaneous_reaction"]
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
        print("start initializing")
        self.initialize()

    def test(self):
        count = 0
        for node in self.nodeList:
            for downRec in node.getDownedge():
                count = count + len(downRec.getpro())
        print ("average next node count : ", count/len(self.nodeList))

    def main(self, mode, input_species, output_species):
        if (mode == 0):
            self.main_default(input_species, output_species)
        elif (mode == 1):
            self.main_bd(input_species, output_species)
        elif (mode == 2):
            self.main_XOR(input_species, output_species)

    def main_bd(self, input_species, output_species):
        self.input_species  = input_species

        for node in output_species:
            if node in self.mapToNode:
                self.output_species.append(self.mapToNode[node])
        for node in self.output_species:
            node.label.add("output")
            SearchUp(node, 6)
        for node in self.input_species:
            print (node.name , node.layer)

        candidate = []
        self.BuildStart(input_species)
        candidate0 = SearchDown(input_species[0], "A", self.up_limit, 1, self.input_species, self.output_species, self.out, self.ans)
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord("A")
        self.ClearPath()
        self.ClearLevel()
        print(len(candidate0))
        self.BuildStart(input_species)
        candidate1 = SearchDown(input_species[1], "A", self.up_limit, 1, self.input_species, self.output_species, self.out, self.ans)

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
                
        for i in range(len(candidate)):
            self.ClearPathA1()
            self.ClearPathA2()
            self.ClearPath()
            self.ClearLevel()
            candidate[i].CopyToPath("A")
            candidateB = SearchDown(candidate[i], "B", self.up_limit, 1, self.input_species, self.output_species, False, self.ans)
            print ("candidateB : ", len(candidateB))
            print ("candidate name : ", candidate[i].name)
            for candB in candidateB:
                candB[0].CopyToRecord("A1")
            
            for j in range(len(candidate)):
                if candidate[i] == candidate[j]:
                    continue
                self.ClearPath()
                self.ClearLevel()
                candidate_B = []
                candidate[j].CopyToPath("A")
                candidate_Bside = SearchDown(candidate[j], "B_side", self.up_limit, 1, self.input_species, self.output_species, False, self.ans)
                print ("start merging B")

                for species in candidate_Bside:
                    targetRec = self.mapToEdge[species[1]]
                    targetB = targetRec.rec_node
                    targetB.path = species[0].path
                    startTime = time.process_time()
                    merge = self.MergeBside(targetB.path, targetB.pathA1, targetB, species[0])
                    endTime = time.process_time()
                    self.time = self.time + endTime - startTime
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

                for target in candidate_B:
                    self.ClearPath()
                    self.ClearLevel()
                    target[0].CopyToPath("A2")
                    target[1].mark = 1
                    target[2].getenz().mark = 1
                    output = SearchDown(target[0], "output", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
                    target[1].mark = 0
                    target[2].getenz().mark = 0

                    print ("target : ", target[0].show(), target[0].layer)
                    for outNode in self.output_species:
                        tmp = []
                        if len(outNode.path) == 0:
                            continue
                        print ("outNode : ", outNode.show())
                        for path in outNode.path:
                            startTime1 = time.process_time()
                            Check_Vio = self.CheckVio(path["struc"][0], path["struc"][1], target[0], target[1], path, None)
                            endTime1 = time.process_time()
                            self.time = self.time + endTime1 - startTime1
                            if not Check_Vio:
                                continue
                            
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
                                        startTime1 = time.process_time()
                                        Check_Vio = self.CheckVio(path["struc"][0], path["struc"][1], target[0], target[1], path, upRec)
                                        endTime1 = time.process_time()
                                        self.time = self.time + endTime1 - startTime1
                                        if Check_Vio:
                                            thres.append(upRec)
                                            self.out = True
                                for x in thres:
                                    allnodes = set()
                                    pathlist = path["pathlist"].copy()
                                    for rec in pathlist:
                                        for node in rec.getrea():
                                            allnodes.add(node)
                                        for node in rec.getpro():
                                            allnodes.add(node)
                                        allnodes.add(rec.getenz())

                                    pathlist.append(x)
                                    for node in x.getrea():
                                        allnodes.add(node)
                                    for node in x.getpro():
                                        allnodes.add(node)
                                    allnodes.add(x.getenz())
                                    allrec = self.CollectAll(allnodes, pathlist)
                                    c = len(pathlist)+len(allrec)

                                    self.count += 1
                                    if c in self.distribute:
                                        self.distribute[c] = self.distribute[c] + 1
                                    else:
                                        self.distribute[c] = 1
                                    if c < self.ans:
                                        self.ans = c
                                        self.record_struc = path["struc"]
                                        self.record_node = path["pathnode"].copy()
                                        self.record_ans = pathlist
                                        self.record_oth = allrec.copy()
        self.Summary()


    def main_default(self, input_species, output_species):
        self.input_species  = input_species

        for node in output_species:
            if node in self.mapToNode:
                self.output_species.append(self.mapToNode[node])
        for node in self.output_species:
            node.label.add("output")
            SearchUp(node, 6)
        for node in self.input_species:
            print (node.name , node.layer)

        candidate = []
        self.BuildStart(input_species)
        # input_species[0].CopyToPath("A")
        candidate0 = SearchDown(input_species[0], "A", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord("A")

        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0))

        # for inNode in input_species:
        self.BuildStart(input_species)
        # input_species[1].CopyToPath("A")
        candidate1 = SearchDown(input_species[1], "A", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
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
                
        for i in range(len(candidate)):
            self.ClearPathA1()
            self.ClearPathA2()
            self.ClearPath()
            self.ClearLevel()
            candidate[i].CopyToPath("A")
            candidateB = SearchDown(candidate[i], "B", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
            print ("candidateB : ", len(candidateB))
            print ("candidate name : ", candidate[i].name)
            for candB in candidateB:
                candB[0].CopyToRecord("A1")
            
            for j in range(len(candidate)):
                if candidate[i] == candidate[j]:
                    continue
                self.ClearPath()
                self.ClearLevel()
                candidate_B = []
                candidate[j].CopyToPath("A")
                candidate_Bside = SearchDown(candidate[j], "B_side", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
                print ("start merging B")
                # for species in candidate_Bside:
                #     targetRec = self.mapToEdge[species[1]]
                #     targetB = targetRec.rec_node
                #     targetB.pathA2 = []
                #     for spath in species[0].path:
                #         tmp    = {}
                #         tmp["related"]  = spath["related"].copy()
                #         tmp["pathlist"] = spath["pathlist"].copy()
                #         tmp["pathenz"]  = spath["pathenz"].copy()
                #         tmp["pathnode"] = spath["pathnode"].copy()
                #         targetB.pathA2.append(tmp)

                # for species in candidate_Bside:
                #     targetRec = self.mapToEdge[species[1]]
                #     targetB = targetRec.rec_node

                #     # print ("target : ", targetB.show(), targetB.layer)
                #     # print ("species[0] : ", species[0].show(), species[0].layer)
                #     self.ClearOut()
                #     # if targetB.name != "phosphate":
                #     #     continue
                #     # for pathA2 in targetB.pathA2:
                #     #     print ("***********************a new path2******************************")
                #     #     for reaction in pathA2["pathlist"]:
                #     #         print (reaction.show(0))
                #     # for pathA2 in targetB.pathA1:
                #     #     print ("***********************a new path1******************************")
                #     #     for reaction in pathA2["pathlist"]:
                #     #         print (reaction.show(0))
                #     self.CheckVio(targetB.pathA2, targetB.pathA1, targetB, species[0], targetRec)
                #     if time.process_time() - startTime > 3600:
                #         break
                    # print ("target : ", targetB.show(), targetB.layer)
                for species in candidate_Bside:
                    targetRec = self.mapToEdge[species[1]]
                    targetB = targetRec.rec_node
                    targetB.path = species[0].path
                    startTime = time.process_time()
                    merge = self.MergeBside(targetB.path, targetB.pathA1, targetB, species[0])
                    endTime = time.process_time()
                    self.time = self.time + endTime - startTime
                    if merge != []:
                        print("merge success")
                        targetB.pathA2 = merge
                        for path in targetB.pathA2:
                            # assert(len(path["struc"]) == 2)
                            path["pathlist"].append(targetRec)
                            for pro in targetRec.getpro():
                                path["related"].add(pro)
                            for rea in targetRec.getrea():
                                if rea not in path["pathnode"]:
                                    path["related"].add(rea)
                        candidate_B.append((targetB, species[0], targetRec))
                print (len(candidate_B))

                for target in candidate_B:
                    self.ClearPath()
                    self.ClearLevel()
                    target[0].CopyToPath("A2")
                    target[1].mark = 1
                    target[2].getenz().mark = 1
                    output = SearchDown(target[0], "output", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
                    target[1].mark = 0
                    target[2].getenz().mark = 0

                    print ("target : ", target[0].show(), target[0].layer)
                    for outNode in self.output_species:
                        tmp = []
                        if len(outNode.path) == 0:
                            continue
                        print ("outNode : ", outNode.show())
                        for path in outNode.path:
                            startTime1 = time.process_time()
                            Check_Vio = self.CheckVio(path["struc"][0], path["struc"][1], target[0], target[1], path, None)
                            endTime1 = time.process_time()
                            self.time = self.time + endTime1 - startTime1
                            if not Check_Vio:
                                continue
                            
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
                                        startTime1 = time.process_time()
                                        Check_Vio = self.CheckVio(path["struc"][0], path["struc"][1], target[0], target[1], path, upRec)
                                        endTime1 = time.process_time()
                                        self.time = self.time + endTime1 - startTime1
                                        if Check_Vio:
                                            thres.append(upRec)
                                # ans = 100
                                for x in thres:
                                    allnodes = set()
                                    pathlist = path["pathlist"].copy()
                                    for rec in pathlist:
                                        for node in rec.getrea():
                                            allnodes.add(node)
                                        for node in rec.getpro():
                                            allnodes.add(node)
                                        allnodes.add(rec.getenz())
                                    # allrec_tmp = len(self.CollectAll(allnodes, pathlist))

                                    pathlist.append(x)
                                    for node in x.getrea():
                                        allnodes.add(node)
                                    for node in x.getpro():
                                        allnodes.add(node)
                                    allnodes.add(x.getenz())
                                    allrec = self.CollectAll(allnodes, pathlist)
                                    # if len(allrec) == allrec_tmp:
                                    c = 0
                                    # for rec in path["pathlist"]:
                                    c = c+len(pathlist) 
                                    # for rec in allrec:
                                    c = c+len(allrec)

                                    self.count += 1
                                    if c in self.distribute:
                                        self.distribute[c] = self.distribute[c] + 1
                                    else:
                                        self.distribute[c] = 1
                                    if c < self.ans:
                                        self.ans = c
                                        self.record_struc = path["struc"]
                                        self.record_node = path["pathnode"].copy()
                                        self.record_ans = pathlist
                                        self.record_oth = allrec.copy()
        self.Summary()

    def main_XOR(self, input_species, output_species):
        self.input_species  = input_species

        for node in output_species:
            if node in self.mapToNode:
                self.output_species.append(self.mapToNode[node])
        for node in self.output_species:
            node.label.add("output")
            SearchUp(node, 6)
        for node in self.input_species:
            print (node.name , node.layer)

        candidate = []
        self.BuildStart(input_species)
        # input_species[0].CopyToPath("A")
        candidate0 = SearchDown(input_species[0], "A", self.up_limit, 1, self.input_species, self.output_species, self.out, self.ans)
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord("A")

        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0))

        self.BuildStart(input_species)
        candidate1 = SearchDown(input_species[1], "A", self.up_limit, 1, self.input_species, self.output_species, self.out, self.ans)
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
                
        for i in range(len(candidate)):
            self.ClearPathA1()
            self.ClearPathA2()
            self.ClearPath()
            self.ClearLevel()
            candidate[i].CopyToPath("A")
            candidateB = SearchDown(candidate[i], "B", self.up_limit, 1, self.input_species, self.output_species, self.out, self.ans)
            print ("candidateB : ", len(candidateB))
            print ("candidate name : ", candidate[i].name)
            for candB in candidateB:
                candB[0].CopyToRecord("A1")
            
            for j in range(len(candidateB)):
                self.ClearPath()
                self.ClearLevel()
                candidate_B = []
                candidateB[j].CopyToPath("A1")
                candidate_Bside = SearchDown(candidateB[j], "B_side", self.up_limit, 1, self.input_species, self.output_species, self.out, self.ans)
                print ("start merging B")

                for species in candidate_Bside:
                    targetRec = self.mapToEdge[species[1]]
                    targetB = targetRec.rec_node
                    targetB.path = species[0].path
                    startTime = time.process_time()
                    merge = self.MergeBside(targetB.path, targetB.pathA1, targetB, species[0])
                    endTime = time.process_time()
                    self.time = self.time + endTime - startTime
                    if merge != []:
                        print("merge success")
                        targetB.pathA2 = merge
                        for path in targetB.pathA2:
                            # assert(len(path["struc"]) == 2)
                            path["pathlist"].append(targetRec)
                            for pro in targetRec.getpro():
                                path["related"].add(pro)
                            for rea in targetRec.getrea():
                                if rea not in path["pathnode"]:
                                    path["related"].add(rea)
                        candidate_B.append((targetB, species[0], targetRec))
                print (len(candidate_B))

                for target in candidate_B:
                    self.ClearPath()
                    self.ClearLevel()
                    target[0].CopyToPath("A2")
                    target[1].mark = 1
                    target[2].getenz().mark = 1
                    output = SearchDown(target[0], "output", self.up_limit, 3, self.input_species, self.output_species, self.out, self.ans)
                    target[1].mark = 0
                    target[2].getenz().mark = 0

                    print ("target : ", target[0].show(), target[0].layer)
                    for outNode in self.output_species:
                        tmp = []
                        if len(outNode.path) == 0:
                            continue
                        print ("outNode : ", outNode.show())
                        for path in outNode.path:
                            startTime1 = time.process_time()
                            Check_Vio = self.CheckVio(path["struc"][0], path["struc"][1], target[0], target[1], path, None)
                            endTime1 = time.process_time()
                            self.time = self.time + endTime1 - startTime1
                            if not Check_Vio:
                                continue
                            
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
                                        startTime1 = time.process_time()
                                        Check_Vio = self.CheckVio(path["struc"][0], path["struc"][1], target[0], target[1], path, upRec)
                                        endTime1 = time.process_time()
                                        self.time = self.time + endTime1 - startTime1
                                        if Check_Vio:
                                            thres.append(upRec)
                                # ans = 100
                                for x in thres:
                                    allnodes = set()
                                    pathlist = path["pathlist"].copy()
                                    for rec in pathlist:
                                        for node in rec.getrea():
                                            allnodes.add(node)
                                        for node in rec.getpro():
                                            allnodes.add(node)
                                        allnodes.add(rec.getenz())
                                    # allrec_tmp = len(self.CollectAll(allnodes, pathlist))

                                    pathlist.append(x)
                                    for node in x.getrea():
                                        allnodes.add(node)
                                    for node in x.getpro():
                                        allnodes.add(node)
                                    allnodes.add(x.getenz())
                                    allrec = self.CollectAll(allnodes, pathlist)
                                    # if len(allrec) == allrec_tmp:
                                    c = 0
                                    # for rec in path["pathlist"]:
                                    c = c+len(pathlist) 
                                    # for rec in allrec:
                                    c = c+len(allrec)

                                    self.count += 1
                                    if c in self.distribute:
                                        self.distribute[c] = self.distribute[c] + 1
                                    else:
                                        self.distribute[c] = 1
                                    if c < self.ans:
                                        self.ans = c
                                        self.record_struc = path["struc"]
                                        self.record_node = path["pathnode"].copy()
                                        self.record_ans = pathlist
                                        self.record_oth = allrec.copy()
        self.Summary()


    def Summary(self):
        print ("ans : ", self.ans, len(self.record_ans), len(self.record_oth))
        c = 0
        collect_rea = set()
        for rec in self.record_ans:
            for rea in rec.getrea():
                collect_rea.add(rea)
            print (rec.show(c))
            c = c+1
        print ("******************pathnode****************")
        for node in self.record_node:
            print( node.show() )
        for node in collect_rea:
            print("present(", node.show(), ", 1).")
        print ("******************other reactions****************")
        for rec in self.record_oth:
            print (rec.show(c))
            c = c+1
        print ("total solution number : ", self.count)
        print ("total solution distribution : ", self.distribute)

    def Merge(self, path, record, targetrec):
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

        tmp["pathlist"].append(targetrec)
        for pro in targetrec.getpro():
            tmp["related"].add(pro)
        for rea in targetrec.getrea():
            if rea not in tmp["pathnode"]:
                tmp["related"].add(rea)
        return tmp

    def MarkPath(self, pw_p, nw_p, out_p, vp, vn):
        for node in out_p["pathnode"]:
            if node.mark == 0:
                node.mark = 6
        for rec in out_p["pathlist"]:
            for rea in rec.getrea():
                if rea.mark == 0:
                    rea.mark = 5
        for rec in pw_p["pathlist"]:
            if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
                rec.getenz().mark = 1
            for rea in rec.getrea():
                if ((rea != self.input_species[0]) and (rea != self.input_species[1])):
                    rea.mark = 5
        for rec in nw_p["pathlist"]:
            if (self.input_species[0] in rec.getrea()) or (self.input_species[1] in rec.getrea()):
                rec.getenz().mark = 2
            for rea in rec.getrea():
                if ((rea != self.input_species[0]) and (rea != self.input_species[1])):
                    rea.mark = 5
        for node in pw_p["pathnode"]:
            if (node not in self.input_species):
                node.mark = 1
        for node in nw_p["pathnode"]:
            if (node not in self.input_species):
                node.mark = 2
        vp.mark = 3
        vn.mark = 4

    def ClearMP(self, pw_p, nw_p, out_p, vp, vn):
        for node in self.nodeList:
            node.mark    = 0

    # def CheckVio(self, pw_p, nw_p, out_p, vp, vn):
    # def CheckVio(self, label0, label1, node0, node1, targetRec):
    #     merge = self.MergeBside(label0, label1, node0, node1)
    #     if merge == []:
    #         return False
        
    #     print("merge success")
    #     # path_tmp = []
    #     for cand in merge:
    #         path = cand[0]
    #         sidepath = cand[1]
    #         self.ClearPath()
    #         self.ClearLevel()
    #         tmp = self.Merge(path, sidepath, targetRec)
    #         # for rec in tmp["pathlist"]:
    #         #     print (rec.show(0))
    #         node0.path.append(tmp)
    #         output = SearchDown(node0, "output", self.up_limit, 3, self.input_species, self.output_species)
    #         # print ("**********************output : ",len(output))
    #         for outNode in output:
    #             for outpath in outNode[0].path:                
    #                 self.MarkPath(path, sidepath, outpath, node0, node1)
    #                 outNode[0].mark = 6
    #                 for reaction in outpath["pathlist"]:
    #                     if node0 in reaction.getrea() and node1 not in reaction.getrea():
    #                         reaction.getenz().mark = 6
    #                 # for rec in outpath["pathlist"]:
    #                 #     print (rec.show(0))
    #                 if self.CheckAll(outpath):
    #                     print ("find one outpath")
    #                     outNode[0].pathOut.append(outpath)
    #                     self.ClearMP(path, sidepath, outpath, node0, node1)
    #                     break
    #                 self.ClearMP(path, sidepath, outpath, node0, node1)
    #                 outNode[0].mark = 0
    #     return True

    def CheckVio(self, path, sidepath, node0, node1, outpath, upRec):
        self.ClearMP(path, sidepath, outpath, node0, node1)
        self.MarkPath(path, sidepath, outpath, node0, node1)
        if upRec != None:
            upRec.getenz().mark = 2
        for reaction in outpath["pathlist"]:
            if node0 in reaction.getrea() and node1 not in reaction.getrea():
                reaction.getenz().mark = 6
        # for rec in outpath["pathlist"]:
        #     print (rec.show(0))
        if self.CheckAll(outpath,upRec):
            self.ClearMP(path, sidepath, outpath, node0, node1)
            return True
        self.ClearMP(path, sidepath, outpath, node0, node1)
        return False

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
                    tmp = self.Merge_tmp(path, sidepath)
                    if self.CheckAll_tmp(tmp):
                        pw_p = {}
                        pw_p["pathlist"] = path["pathlist"].copy()
                        pw_p["pathnode"] = path["pathnode"].copy()
                        pw_p["related"] = path["related"].copy()
                        pw_p["pathenz"] = path["pathenz"].copy()
                        nw_p = {}
                        nw_p["pathlist"] = sidepath["pathlist"].copy()
                        nw_p["pathnode"] = sidepath["pathnode"].copy()
                        nw_p["related"] = sidepath["related"].copy()
                        nw_p["pathenz"] = sidepath["pathenz"].copy()
                        tmp["struc"] = ((pw_p, nw_p))
                        # assert (type(tmp["struc"][0]) == type({}))
                        # assert (type(tmp["struc"][1]) == type({}))
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

    def Merge_tmp(self, path, record):
        tmp    = {}
        tmp["struc"]  = path["struc"]
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

    def CheckAll_tmp(self, path):
        allnodes = set()
        allnodes.add(self.mapToNode["H2O"])
        allnodes.add(self.mapToNode["spontaneous_reaction"])
        for rec in path["pathlist"]:
            for node in rec.getrea():
                allnodes.add(node)
            for node in rec.getpro():
                allnodes.add(node)
            allnodes.add(rec.getenz())

        allrec = self.CollectAll_tmp(allnodes, path["pathlist"])

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

    def CollectAll_tmp(self, allnodes, pathlist):
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

    # def CreatePath(self, pathnode, downRec):
    #     related  = set()
    #     enz      = []
    #     enz.append(downRec.getenz())
    #     for rea in downRec.getrea():
    #         if rea not in pathnode:
    #             related.add(rea)
    #     for pro in downRec.getpro():
    #         related.add(pro)
    #     path = {}
    #     path["pathnode"] = pathnode
    #     path["related"]  = related
    #     path["pathlist"] = []
    #     path["pathlist"].append(downRec)
    #     path["pathenz"]  = enz

    #     return path

    def BuildStart(self, input_species):
        for node in input_species:
            tmp             = {}
            tmp["struc"]    = ()
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

    def CheckAll(self, path, upRec):
        # return True
        allnodes = set()
        allnodes.add(self.mapToNode["H2O"])
        allnodes.add(self.mapToNode["spontaneous_reaction"])
        pathlist = path["pathlist"].copy()
        for rec in path["pathlist"]:
            for node in rec.getrea():
                allnodes.add(node)
            for node in rec.getpro():
                allnodes.add(node)
            allnodes.add(rec.getenz())
        if upRec != None:
            for node in upRec.getrea():
                allnodes.add(node)
            for node in upRec.getpro():
                allnodes.add(node)
            allnodes.add(upRec.getenz())
            pathlist.append(upRec)

        allrec = self.CollectAll(allnodes, pathlist)

        for rec in allrec:
            cr = False
            c_c = 0
            cr3 = 0
            cr4 = 0
            for rea in rec.getrea():
                if rea.mark == 3:
                    cr3 = 1
                if rea.mark == 4:
                    cr4 = 1
                elif rea.mark != 5:
                    cr = True
                elif rea.mark != 0:
                    c_c += 1
            if rec.getenz().mark != 0:
                cr = True
            if cr == False:
                return False
            if cr4 == 1 and cr3 != 1:
                return False
            if cr4 == 1 and cr3 == 1:
                c_c += 1
            if c_c < 2:
                for pro in rec.getpro():
                    if pro.mark == 1 or pro.mark == 3:
                        if pro.pro_m[2]:
                            return False
                    elif pro.mark == 2:
                        if pro.pro_m[1] or pro.pro_m[3]:
                            return False
                    elif pro.mark == 6:
                        if pro.pro_m[1] or pro.pro_m[2]:
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
    
    def CollectAll(self, allnodes, pathlist):
        self.ClearMark()
        temp = ""
        count = 1
        allnode = allnodes.copy()
        allnode.add(self.mapToNode["H2O"])
        allnode.add(self.mapToNode["spontaneous_reaction"])
        allrec  = set()
        for reaction in pathlist:
            c = [False]*7
            c_c = 0
            for rea in reaction.getrea():
                if rea.mark == 0:
                    c = rea.pro_m
                else:
                    c[rea.mark] = True
                    if rea.mark != 5 or rea.mark != 4:
                        c_c += 1
            if reaction.getenz().mark != 0 and reaction.getenz().mark != 6:
                c[reaction.getenz().mark] = True
                c_c += 1

            if c_c < 2:
                for pro in reaction.getpro():
                    if pro.mark == 0:
                        for pm in range(7):
                            pro.pro_m[pm] = pro.pro_m[pm] | c[pm]
        while count != 0:
            tmp = len(allrec)
            for enz in allnodes:
                for reaction in enz.getCatedge():
                    if (reaction not in pathlist) and (reaction not in allrec) and (reaction.activated(allnode)):
                        c = [False]*7
                        c_c = 0
                        for rea in reaction.getrea():
                            if rea.mark == 0:
                                c = rea.pro_m
                            else:
                                c[rea.mark] = True
                                if rea.mark != 5 or rea.mark != 4:
                                    c_c += 1
                        if reaction.getenz().mark != 0 and reaction.getenz().mark != 6:
                            c[reaction.getenz().mark] = True
                            if reaction.getenz().mark != 0:
                                c_c +=1
                        if c_c < 2:
                            for pro in reaction.getpro():
                                for pm in range(7):
                                    if pro.mark == 0:
                                        pro.pro_m[pm] = pro.pro_m[pm] | c[pm]
                                allnode.add(pro)
                        allrec.add(reaction)
            count = len(allrec) - tmp
        return allrec

    def ClearMark(self):
        for node in self.nodeList:
            node.pro_m    = [False]*7

    def ClearPath(self):
        for node in self.nodeList:
            node.path.clear()

    def ClearPathA1(self):
        for node in self.nodeList:
            node.pathA1.clear()

    def ClearPathA2(self):
        for node in self.nodeList:
            node.pathA2.clear()
            if "B_side" in node.label:
                node.label.remove("B_side")
            node.labelB_side.clear()
        for edge in self.edgeList:
            edge.rec_node   = None

    def ClearRecordA(self):
        for node in self.nodeList:
            node.recordA.clear()

    def ClearOut(self):
        for node in self.nodeList:
            node.pathOut.clear()

    def ClearLevel(self):
        for node in self.nodeList:
            node.level = 0
