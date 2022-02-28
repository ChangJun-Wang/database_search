from node import node

from edge import edge

from bfs import BFS

class search:
    def __init__(self):
        self.c = 0
        self.mapTotar    = {}
        self.core_node   = []
        self.tmp_path    = {}
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
        bfs = BFS(input_species[0], "A", self.mapToCnode, self.mapToClist)
        bfs.BuildStart(input_species)
        candidate0 = bfs.search()
        # record all the returned type A
        for species in candidate0:
            species[0].CopyToRecord()
        self.ClearPath()
        self.ClearLevel()

        print(len(candidate0))

        bfs = BFS(input_species[1], "A", self.mapToCnode, self.mapToClist)
        candidate1 = bfs.search()
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

        print("start merging type C_side")
        for node in candidate: 
            bfs0 = BFS(node, "C_side", self.mapToCnode, self.mapToClist)
            candidate0 = bfs0.search()
            for node0 in candidate0:
                allnodes = set()
                for path0 in node0[0].path:
                    pathnode = self.mapToCnode[str(node0[1])].copy()
                    related  = []
                    enz      = []
                    for rec in self.mapToClist[str(node0[1])]:
                        enz.append(rec.getenz())
                        for rea in rec.getrea():
                            if rea not in pathnode:
                                related.append(rea)
                    path1 = {}
                    path1["pathnode"] = pathnode
                    path1["related"]  = related
                    path1["pathlist"] = self.mapToClist[str(node0[1])]
                    path1["pathenz"]  = enz
                    for species in path0["pathnode"]:
                        species.mark = 1
                    for species in path1["pathnode"]:
                        species.mark = 2
                    if (self.CheckMerge(path0, path1) 
                        and self.CheckLabelC(path0, path1)
                        and node0[0].name != "H+"):
                        path_tmp = self.Merge(path0, path1)
                        for rec in path_tmp["pathlist"]:
                            print(rec.show())
                            for species in rec.getrea():
                                allnodes.add(species)
                            for species in rec.getpro():
                                allnodes.add(species)
                            allnodes.add(rec.getenz())
                        print("/////////////////threshold reaction///////////////")
                        for cnode in self.mapToCnode[str(node0[1])][2:]:
                            if cnode != node0[0]:
                                print (cnode.show())
                                print (len(cnode.getUpedge()))
                                for upEdge in cnode.getUpedge():
                                    if (upEdge.getenz().name != "alcohol_dehydrogenase" 
                                        and upEdge.getenz().name != "spontaneous_reaction"):
                                        print(upEdge.show())
                                        for species in upEdge.getrea():
                                            allnodes.add(species)
                                        for species in upEdge.getpro():
                                            allnodes.add(species)
                                        allnodes.add(upEdge.getenz())
                                        break
                        print("/////////////////all related reaction///////////////")
                        print(self.check(allnodes, path_tmp["pathlist"]))
                        for nodes in allnodes:
                            print ("present(" + nodes.show() + ").")
                        return True
                    for species in path0["pathnode"]:
                        species.mark = 0
                    for species in path1["pathnode"]:
                        species.mark = 0


        # for species in candidate:
        #     assert species.recordA != []

        # for i in range(len(candidate)):
        #     for j in range(i+1, len(candidate)):
        #         bfs1 = BFS(candidate[j], "C_side", self.mapToCnode, self.mapToClist)
        #         candidate1 = bfs1.search()
        #         print("start merging type C_side")
        #         for label0 in candidate0:
        #             for label1 in candidate1:
        #                 if int(label0[1]) == int(label1[1]) and label0[0] != label1[0]:
        #                     accept = self.MergeCside(label0, label1)
        #                     if accept != {}:
        #                         allnodes = set()
        #                         result.append(accept)
        #                         for rec in accept["pathlist"]:
        #                             print(rec.show())
        #                             for species in rec.getrea():
        #                                 allnodes.add(species)
        #                             for species in rec.getpro():
        #                                 allnodes.add(species)
        #                             allnodes.add(rec.getenz())

        #                         print("/////////////////all related reaction///////////////")
        #                         print(self.check(allnodes, accept["pathlist"]))
        #                         for node in allnodes:
        #                             print ("present(" + node.show() + ").")
        #                         return True
        #         self.ClearPath()
        #         self.ClearLevel()

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

    def MergeCside(self, label0, label1):
        path_tmp = {}
        for sidepath in label0[0].sidepath:
            for path in label1[0].path:
                for species in sidepath["pathnode"]:
                    species.mark = 1
                for species in path["pathnode"]:
                    species.mark = 2
                if self.CheckMerge(path, sidepath): # and self.CheckPathnode(path, sidepath):
                    path0 = self.Merge(path, sidepath)
                    labelnum = int(label0[1])
                    pathnode = self.mapToCnode[str(labelnum)].copy()
                    related  = []
                    enz      = []
                    for rec in self.mapToClist[str(labelnum)]:
                        enz.append(rec.getenz())
                        for rea in rec.getrea():
                            if rea not in pathnode:
                                related.append(rea)
                    path1 = {}
                    path1["pathnode"] = pathnode
                    path1["related"]  = related
                    path1["pathlist"] = self.mapToClist[str(labelnum)]
                    path1["pathenz"]  = enz
                    if self.CheckMerge(path0, path1) and self.CheckLabelC(path0, path1):
                        path_tmp = self.Merge(path0, path1)
                    for species in sidepath["pathnode"]:
                        species.mark = 0
                    for species in path["pathnode"]:
                        species.mark = 0
        return path_tmp

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

    def CheckLabelC(self, path0, path1):
        AllSpecies = set()

        for rec in path0["pathlist"]:
            for species in rec.getrea():
                AllSpecies.add(species)
            AllSpecies.add(rec.getenz())
            for species in rec.getpro():
                AllSpecies.add(species)
        for rec in path1["pathlist"]:
            for species in rec.getrea():
                AllSpecies.add(species)
            AllSpecies.add(rec.getenz())
            for species in rec.getpro():
                AllSpecies.add(species)

        # for species in path0["related"]:
        #     AllSpecies.add(species)
        # for species in path1["related"]:
        #     AllSpecies.add(species)
        # for species in path0["pathnode"]:
        #     AllSpecies.add(species)
        # for species in path1["pathnode"]:
        #     AllSpecies.add(species)

        AllSpecies.add(self.mapToNode["H2O"])
        AllSpecies.add(self.mapToNode["spontaneous_reaction"])

        for enz in path0["pathenz"]:
            for rec in enz.getCatedge():
                if rec.activated(AllSpecies):
                    if not self.CheckCycle(path0, rec):
                        return False
        for enz in path1["pathenz"]:
            for rec in enz.getCatedge():
                if rec.activated(AllSpecies):
                    if not self.CheckCycle(path1, rec):
                        return False

        # n = 1
        # while n != 0:
        #     AllRec     = set()
        #     tmp        = len(AllSpecies)
        #     for enz in path0["pathenz"]:
        #         for rec in enz.getCatedge():
        #             if rec.activated(AllSpecies):
        #                 AllRec.add(rec)
        #                 if not self.CheckCycle(path0, rec):
        #                     return False
        #     for enz in path1["pathenz"]:
        #         for rec in enz.getCatedge():
        #             if rec.activated(AllSpecies):
        #                 AllRec.add(rec)
        #                 if not self.CheckCycle(path1, rec):
        #                     return False
        #     for rec in AllRec:
        #         for species in rec.getrea():
        #             AllSpecies.add(species)
        #         AllSpecies.add(rec.getenz())
        #         for species in rec.getpro():
        #             AllSpecies.add(species)
        #     n = len(AllSpecies) - tmp

        return True

    def CheckCycle(self, path, rec):
        if rec in path["pathlist"]:
            return True
        c1 = 0
        c2 = 0
        cr = 0
        for rea in rec.getrea():
            if rea.mark == 1:
                c1 = c1 + 1
            elif rea.mark == 2:
                c2 = c2 + 1
        for pro in rec.getpro():
            if pro.mark == 1:
                c1 = c1 + 1
            elif pro.mark == 2:
                c2 = c2 + 1
        if c1 == 0 or c2 == 0:
            return True
        return False

    def check_rea(self, enz, product):
        for catrec in enz.getCatedge():
            if product in catrec.getrea():
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
            tmp             = {}
            tmp["related"]  = set()
            tmp["pathlist"] = []
            tmp["pathenz"]  = []
            tmp["pathnode"] = []
            node.path.append(tmp)

        for node in self.nodeList:
            node.path[0]["related"].add(self.mapToNode["H2O"])
            node.path[0]["pathenz"].append(self.mapToNode["spontaneous_reaction"])

    def ClearRecordA(self):
        for node in self.nodeList:
            node.recordA    = []
            tmp             = {}
            tmp["related"]  = set()
            tmp["pathlist"] = []
            tmp["pathenz"]  = []
            tmp["pathnode"] = []
            node.recordA.append(tmp)
        for node in self.nodeList:
            node.recordA[0]["related"].add(self.mapToNode["H2O"])
            node.recordA[0]["pathenz"].append(self.mapToNode["spontaneous_reaction"])

    def ClearRecordC(self):
        for node in self.nodeList:
            node.recordC = {}

    def ClearLevel(self):
        for node in self.nodeList:
            node.level = 0

    