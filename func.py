import sys
 
# the setrecursionlimit function is
# used to modify the default recursion
# limit set by python. Using this,
# we can increase the recursion limit
# to satisfy our needs
 
sys.setrecursionlimit(10**6)

from node import node

from edge import edge

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
        self.forbid_node.add(self.mapToNode["H2O"])
        self.mapToNode["H2O"].visited = 1
        for nodes in self.nodeList:
            if (nodes.getpin() > 1):
                enz = set()
                for rec in nodes.getUpedge():
                    if rec.getenz()[0] != "spontaneous_reaction":
                        enz.add(rec.getenz()[0])
                self.typeA += 1
                nodes.setpin(len(enz))
                nodes.label.add("A")
                nodes.labelA = str(self.typeA)
            if (nodes.getpout() > 1):
                enz = set()
                for rec in nodes.getDownedge():
                    if rec.getenz()[0] != "spontaneous_reaction":
                        enz.add(rec.getenz()[0])
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
            # if nodes.visited == 0:
            #     # print(nodes.name)
            #     nodes.visited = 1
            #     self.find_sec_part(nodes)
                # if self.typeC == 10:
                #     break
                    # self.typeC += 1
                # self.Cycle_Find(nodes, [], [])
        # for nodes in self.nodeList:
        #     assert len(nodes.labelA) == 0 or len(nodes.labelA) == 1
        #     assert len(nodes.labelB) == 0 or len(nodes.labelB) == 1
        #     assert len(nodes.labelD) == 0 or len(nodes.labelD) == 1
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
            else:
                is_product = 0
                self.reaction = tmpList[0]
                self.mapToEdge[self.reaction] = edge(self.reaction)
                self.edgeList.append(self.mapToEdge[self.reaction])

                self.mapToEdge[self.reaction].addenz(self.mapToNode[enzyme])
                self.mapToNode[enzyme].addCatedge(self.mapToEdge[self.reaction])

                for species in tmpList[1:]:
                    HashKey = HashKey + species[0]
                    if (species == "="):
                        is_product = 1
                    # elif(species == "{r}"):
                    #     self.mapToEdge[self.reaction].reverse = True
                    # elif(species == "{ir}" or species == "?" or species == "more"):
                    #     pass
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
                            # if (enzyme != 'spontaneous_reaction'): # and enzyme not in self.mapToNode[species].enzyme):
                            #     #self.mapToNode[species].enzyme.add(enzyme)
                            self.mapToNode[species].addpout()
                        else:
                            self.mapToEdge[self.reaction].addpro(self.mapToNode[species])
                            self.mapToNode[species].addUpedge(self.mapToEdge[self.reaction])
                            (self.mapToNode[species]).upEdgePar[self.mapToEdge[self.reaction]] = par
                            # if (enzyme != 'spontaneous_reaction'): # and enzyme not in self.mapToNode[species].enzyme):
                            #     #self.mapToNode[species].enzyme.add(enzyme)
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



    def check_rec(self, edge, nodelist):
        temp = []
        for rea in edge.getrea():
            if (rea not in nodelist) and (rea != self.mapToNode["H2O"]):
                temp.append(rea)
        return temp

    def check_activated(self, edgelist):
        assert(len(edgelist) == 3)
        enzyme    = set()
        species   = set()
        temp      = []

        for rec in edgelist:
            enzyme.add(rec.getenz()[0])
            for i in rec.getrea():
                species.add(i)
            for i in rec.getpro():
                species.add(i)
        for enz in enzyme:
            for rec in enz.getCatedge():
                if rec not in edgelist:
                    check = self.check_rec(rec, species)
                    temp.append(check)
                    if check == []:
                        return False
        self.mapToCforbid[self.typeC] = temp
        return True


    def check_legal(self, node1, edge1, node2, edge2):
        for label in edge1.labelC:
            if label in edge2.labelC:
                return False
        for species1 in edge1.getrea():
            if (species1 != node1 and species1 != self.mapToNode["H2O"] 
                and "B" in species1.label):
                for species2 in edge2.getrea():
                    if (species2 != node2 and species2 != species1 
                        and species2 != self.mapToNode["H2O"]
                        and "B" in species2.label):
                        for downRec in species1.getDownedge():
                            if (downRec in species2.getDownedge() 
                                and species1 not in downRec.getpro() and species2 not in downRec.getpro()):
                                edgelist = []
                                edgelist.append(edge1)
                                edgelist.append(edge2)
                                edgelist.append(downRec)
                                if self.check_activated(edgelist):
                                    print("find typeC", self.typeC)
                                    self.mapToClist[self.typeC] = []
                                    self.mapToClist[self.typeC].append(edge1)
                                    self.mapToClist[self.typeC].append(edge2)
                                    self.mapToClist[self.typeC].append(downRec)
                                    self.mapToCnode[self.typeC] = []
                                    self.mapToCnode[self.typeC].append(node1)
                                    self.mapToCnode[self.typeC].append(node2)
                                    self.mapToCnode[self.typeC].append(species1)
                                    self.mapToCnode[self.typeC].append(species2)
                                    assert(len(self.mapToClist[self.typeC]) == 3)
                                    assert(len(self.mapToCnode[self.typeC]) == 4)

                                    node1.label.add("C")
                                    node2.label.add("C")
                                    edge1.label.add("C")
                                    edge2.label.add("C")
                                    node1.labelC.append(self.typeC)
                                    node2.labelC.append(self.typeC)
                                    edge1.labelC.append(self.typeC)
                                    edge2.labelC.append(self.typeC)
                                    for rea in edge1.getrea():
                                        if rea != node1:
                                            rea.label.add("C_side")
                                            rea.labelC_side.append(self.typeC)
                                    for rea in edge2.getrea():
                                        if rea != node2:
                                            rea.label.add("C_side")
                                            rea.labelC_side.append(self.typeC)
                                    self.typeC += 1
                                    return True
        return False

    def find_sec_part(self, node):
        for downRec in node.getDownedge():
            for product in downRec.getpro():
                if product.visited == 0:
                    for d_downRec in product.getDownedge():
                        if node in d_downRec.getpro():
                            return self.check_legal(node, downRec, product, d_downRec)
        return False        

    def ClearVis(self):
        for nodes in self.nodeList:
            nodes.visited  = 0
            nodes.recStack = False
        for rec in self.edgeList:
            assert(len(rec.getenz()) == 1)
            rec.visited  = 0
            rec.recStack = False

    def ClearTmp(self):
        for nodes in self.nodeList:
            nodes.tmp_vis  = 0
            nodes.recStack = False
        for rec in self.edgeList:
            assert(len(rec.getenz()) == 1)
            rec.tmp_vis  = 0
            rec.recStack = False

    def ClearPath(self):
        pass


    def BFS_all(self, input_species, notes):
        BFS    = [input_species[0]]
        tmp    = {}
        tmp["related"]  = set()
        tmp["pathlist"] = []
        tmp["pathenz"]  = []
        tmp["pathnode"] = []
        tmp["pathnode"].append(input_species[1])
        input_species[0].path.append(tmp)
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                check_downRec = startPro.check_downRec(downRec)
                if check_downRec != []:
                    for product in downRec.getpro():
                        check_pro = product.check_pro(check_downRec)
                        if check_pro != []:  #((product not in self.forbid_node) and (product not in startPro.related)):
                            BFS.append(product)
                            product.add_path(startPro, downRec, check_pro)
                            assert product.path != check_pro
                            assert len(product.path) != 0
                            
                            if ("A" in product.label) and (product.getpin() > 2):
                                labelA = product.labelA
                                # if self.BFS_tmp(input_species, target):
                                (result, findC, labelC) = self.BFS_findC(product)
                                if result:
                                    result1 = False
                                    find_tmp = None
                                    for i in self.mapToCnode[labelC][2:]:
                                        if i != findC:
                                            (result1, find_tmp) = self.BFS_tmp(input_species, i, findC)
                                    if result1:
                                        allnodes = set()
                                        allnodes.add(self.mapToNode["spontaneous_reaction"])
                                        print("*********** From i1 : ************ ")
                                        for rec in (findC.path[0])["pathlist"]:
                                            print(rec.show())
                                            for node in rec.getrea():
                                                allnodes.add(node)
                                            for node in rec.getenz():
                                                allnodes.add(node)
                                            for node in rec.getpro():
                                                allnodes.add(node)
                                        print("*********** From i1 : ************ ")
                                        for node in (findC.path[0])["pathnode"]:
                                            print(node.show())
                                        print("*********** From i2 : ************ ")
                                        for rec in (find_tmp.path_tmp[0])["pathlist"]:
                                            print(rec.show())
                                            for node in rec.getrea():
                                                allnodes.add(node)
                                            for node in rec.getenz():
                                                allnodes.add(node)
                                            for node in rec.getpro():
                                                allnodes.add(node)
                                        # print("*********** Find A : ************ ")
                                        # c = 0
                                        # for rec in product.getUpedge():
                                        #     if c == 4:
                                        #         break
                                        #     c += 1
                                        #     print(rec.show())
                                        print("*********** Find C : ************ ")
                                        for rec in self.mapToClist[labelC]:
                                            print(rec.show())
                                            for node in rec.getrea():
                                                allnodes.add(node)
                                            for node in rec.getenz():
                                                allnodes.add(node)
                                            for node in rec.getpro():
                                                allnodes.add(node)

                                        print("*********** Threshold : ************ ")
                                        for i in self.mapToCnode[labelC][2:]:
                                            if i != findC:
                                                print(i.show())
                                                for rec in i.getUpedge():
                                                    if (rec.getenz()[0] != self.mapToNode["spontaneous_reaction"]):
                                                        if self. check_rea(rec.getenz()[0], i):
                                                            print(rec.show())
                                                            for node in rec.getrea():
                                                                allnodes.add(node)
                                                            for node in rec.getenz():
                                                                allnodes.add(node)
                                                            for node in rec.getpro():
                                                                allnodes.add(node)
                                                            break

                                        print("/////////////////all related reaction///////////////")
                                        print(self.check(allnodes))
                                        return True
        return False

    def check_rea(self, enz, product):
        for catrec in enz.getCatedge():
            if product in catrec.getrea() or self.mapToNode["NAD+"] in catrec.getrea() or self.mapToNode["H2O2"] in catrec.getpro():
                return False
        return True

    def check_cycle(self, path, c_side):
        for node in self.mapToCnode[c_side]:
            if node in path["related"]:
                return False
        for nodelist in self.mapToCforbid[c_side]:
            c = 0
            for node in nodelist:
                if node in path["related"]:
                    c += 1
            if c == len(nodelist):
                return False
        return True

    def BFS_findC(self, input_node):
        BFS    = [input_node]
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                check_downRec = startPro.check_downRec(downRec)
                if check_downRec != []:
                    for product in downRec.getpro():
                        check_pro = product.check_pro(check_downRec)
                        if check_pro != []:
                            BFS.append(product)
                            if ("C_side" in product.label):
                                for c_side in product.labelC_side:
                                    for path in product.path:
                                        if self.check_cycle(path, c_side):
                                            product.add_path(startPro, downRec, check_pro)
                                            return (True, product, c_side)
                            product.add_path(startPro, downRec, check_pro)
        return (False, node(None), None)


    def BFS_tmp(self, input_species, target, relate):
        # targetlist   = target.path[0]["pathnode"][2:]
        targetlist = []
        targetlist.append(target)
        tar_relate   = relate.path[0]["related"]
        tar_relate.add(relate)
        BFS    = [input_species[1]]
        tmp    = {}
        tmp["related"]  = set()
        tmp["pathlist"] = []
        tmp["pathenz"]  = []
        tmp["pathnode"] = []
        tmp["pathnode"].append(input_species[0])
        input_species[1].path_tmp.append(tmp)
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                #check whether this reaction has conflict to previous reactions and input
                check_downRec = startPro.check_downRectmp(downRec, relate)
                if check_downRec != []:
                    for product in downRec.getpro():
                        check_pro = product.check_pro(check_downRec)
                        if check_pro != []:
                            if product in targetlist:
                                product.add_pathtmp(startPro, downRec, check_pro)
                                return (True, product)
                            if product not in tar_relate:
                                BFS.append(product)
                                product.add_pathtmp(startPro, downRec, check_pro)
                                assert product.path_tmp != startPro.path_tmp
                                assert product.path_tmp != check_downRec
                                assert len(product.path_tmp) != 0
        return (False, node(None))

    def check(self, allnodes):
        temp = ""
        for enz in allnodes:
            for reaction in enz.getCatedge():
                count = 0
                for reactant in reaction.getrea():
                    if reactant in allnodes:
                        count += 1
                if count == len(reaction.getrea()):
                    temp += (reaction.show()) + str(count)
        return temp
        