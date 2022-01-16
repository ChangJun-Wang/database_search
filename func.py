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
                nodes.labelA.append("A" + str(self.typeA))
            if (nodes.getpout() > 1):
                enz = set()
                for rec in nodes.getDownedge():
                    if rec.getenz()[0] != "spontaneous_reaction":
                        enz.add(rec.getenz()[0])
                self.typeB += 1
                nodes.setpout(len(enz))
                nodes.label.add("B")
                nodes.labelB.append("B" + str(self.typeB))
            if (nodes.getpout() > 0) and (nodes.getpin() > 0):
                self.typeD += 1
                nodes.label.add("D")
                nodes.labelD.append("D" + str(self.typeD))
            if (len(nodes.getCatedge()) > 1):
                self.typeE += 1
                nodes.label.add("E")
                nodes.labelE.append("E" + str(self.typeE))
            # if nodes.visited == 0:
            #     # print(nodes.name)
            #     nodes.visited = 1
            #     self.find_sec_part(nodes)
                # if self.typeC == 10:
                #     break
                    # self.typeC += 1
                # self.Cycle_Find(nodes, [], [])

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

    def check_species(self, pathnode, downRec):
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

    def check_downRec(self, startPro, downRec):
        tmp = []
        for path in startPro.path:
            if self.check_species(path["pathnode"], downRec):
                tmp.append(path)
        return tmp

    def add_path(self, startPro, product, downRec, downPath):
        for path in downPath:
            tmp_path = {}
            tmp_path["related"]  = set()
            tmp_path["pathnode"] = []
            tmp_path["pathlist"] = []
            for i in path["related"]:
                tmp_path["related"].add(i)
            for i in downRec.getrea():
                tmp_path["related"].add(i)
            for i in path["pathlist"]:
                tmp_path["pathlist"].append(i)
            tmp_path["pathlist"].append(downRec)
            for i in path["pathnode"]:
                tmp_path["pathnode"].append(i)
            tmp_path["pathnode"].append(startPro)
            product.path.append(tmp_path)

    def BFS_all(self, input_species, notes):
        BFS    = [input_species[0]]
        tmp    = {}
        tmp["related"]  = set()
        tmp["pathlist"] = []
        tmp["pathnode"] = []
        tmp["pathnode"].append(input_species[1])
        input_species[0].path.append(tmp)
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                #check whether this reaction has conflict to previous reactions and input
                check_downRec = self.check_downRec(startPro, downRec)
                if check_downRec != []:
                    for product in downRec.getpro():
                        if ((product not in self.forbid_node) and (product not in startPro.related)):
                            BFS.append(product)
                            self.add_path(startPro, product, downRec, check_downRec)
                            assert product.path != startPro.path
                            assert product.path != check_downRec
                            assert len(product.path) != 0
                            if ("A" in product.label) and (product.getpin() > 2):
                                labelA = product.labelA[0]
                                (result, findC, labelC) = self.BFS_findC(product)
                                if result:
                                    for rec in (findC.path[0])["pathlist"]:
                                        print(rec.show())
                                    print("*********** Find A : ************ \n")
                                    for rec in product.getUpedge():
                                        print(rec.show())
                                    print("*********** Find C : ************ \n")
                                    for rec in self.mapToClist[labelC]:
                                        print(rec.show())
                                    return True
        return False

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
                check_downRec = self.check_downRec(startPro, downRec)
                if check_downRec != []:
                    for product in downRec.getpro():
                        if ((product not in self.forbid_node) and (product not in startPro.related)):
                            BFS.append(product)
                            self.add_path(startPro, product, downRec, check_downRec)
                            if ("C_side" in product.label):
                                for c_side in product.labelC_side:
                                    for path in product.path:
                                        if self.check_cycle(path, c_side):
                                            return (True, product, c_side)
        return (False, node(None), None)


    def sec_input(self, input_species, notes):
        count = 2
        while count > 0:
            (result0, product0, target0) = self.BFS_tmp(input_species, [], [])
            pathlist = product0.tmp_list
            pathnode = product0.tmp_node
            if result0:
                print("///////////find first one//////////////")
                for rec in product0.pathlist:
                    print(rec.show())
                for rec in pathlist:
                    print(rec.show())
                count -= 1
                print("/////////////////end////////////////////")
                target0.show()
            (result1, product1) = self.BFS_findA(target0, input_species, notes)
            if result1:
                print("///////////find second one//////////////")
                for rec in product1.tmp_list:
                    print(rec.show())
                print("/////////////////end////////////////////")
                count -= 1
                c = 0
                print("///////////selection for threshold//////////////")
                for rec in self.mapToNode["acceptor"].getUpedge():
                    if (self.mapToNode["spontaneous_reaction"] not in rec.getenz() 
                        and self.mapToNode["H+"] not in rec.getrea()
                        and self.mapToNode["reduced_acceptor"] not in rec.getrea()
                        ):
                        print(rec.show())
                        c += 1
                    if c == 3:
                        break
                print("/////////////////end////////////////////")
            self.ClearTmp()



    def BFS_tmp(self, input_species, pathlist, pathnode):
        startPro = input_species[1]
        BFS    = [startPro]
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                if downRec.tmp_vis == 0:
                    downRec.tmp_vis = 1
                    mark = True

                    #check whether this reaction has conflict to previous reactions and input
                    for reactant in downRec.getrea():
                        if reactant == input_species[0] or reactant.name == "NAD+"  or reactant.name == "H+":
                            mark = False
                    for product in downRec.getpro():
                        if product in startPro.tmp_node or product == startPro or product == input_species[0]:
                            mark = False
                        if product in pathnode:
                            mark = False
                    for enzyme in downRec.getenz():
                        if enzyme in startPro.tmp_node:
                            mark = False
                        if enzyme in pathnode:
                            mark = False
                    if mark:
                        for product in downRec.getpro():
                            if (product.tmp_vis == 0 and product not in self.forbid_node):
                                product.tmp_vis = 1
                                BFS.append(product)
                                for i in startPro.forbid:
                                    product.tmp_forbid.add(i)
                                for i in downRec.getrea():
                                    product.tmp_forbid.add(i)
                                for i in startPro.tmp_list:
                                    (product.tmp_list).append(i)
                                (product.tmp_list).append(downRec)
                                assert product.tmp_list != startPro.tmp_list
                                for i in startPro.tmp_node:
                                    product.tmp_node.append(i)
                                product.tmp_node.append(startPro)
                                assert product.tmp_node != startPro.tmp_node
                                #if ("A" in product.label):
                                if len(product.tmp_label) == 3:
                                    product.tmp_node.append(product)
                                    tmp_str = ""
                                    for i in product.tmp_label:
                                        tmp_str += i
                                    
                                    target = self.mapTotar[tmp_str]
                                    return (True, product, target)
        return (False, node(None), node(None))


    def BFS_findA(self, target, input_species, notes):
        BFS    = [input_species[0]]
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                if downRec.visited == 0:
                    downRec.visited = 1
                    mark = True
                    #check whether this reaction has conflict to previous reactions and input
                    for reactant in downRec.getrea():
                        if reactant == input_species[1]:
                            mark = False
                    for product in downRec.getpro():
                        if product in startPro.tmp_node or product == input_species[0]:
                            mark = False
                    for enzyme in downRec.getenz():
                        if enzyme in startPro.tmp_node:
                            mark = False
                    if mark:
                        for product in downRec.getpro():
                            if (product.visited == 0
                                and product not in self.forbid_node):
                                product.visited = 1
                                BFS.append(product)
                                for i in startPro.tmp_forbid:
                                    product.tmp_forbid.add(i)
                                for i in downRec.getrea():
                                    product.tmp_forbid.add(i)
                                for i in startPro.tmp_list:
                                    product.tmp_list.append(i)
                                product.tmp_list.append(downRec)
                                assert product.tmp_list != startPro.tmp_list
                                for i in startPro.tmp_node:
                                    product.tmp_node.append(i)
                                product.tmp_node.append(startPro)
                                assert product.tmp_node != startPro.tmp_node
                                #if ("A" in product.label):
                                if product == target:
                                    return (True, product)
        return (False, node(None))

    def check(self, allnodes):
        temp = ""
        for enz in allnodes:
            for reaction in enz.getCatedge():
                count = 0
                for reactant in reaction.getrea():
                    if reactant in allnodes:
                        count += 1
                for product in reaction.getpro():
                    if reactant in allnodes and count == len(reaction.getrea()):
                        temp += (reaction.show())
                        break
        return temp
        


    # def backTrace(self, startRec, pathlist, pathnode):
    #     #print("in backtrae")
    #     for index in range(len(pathlist)):
    #         if startRec == pathlist[index]:
    #             self.mapToClist[self.typeC] = pathlist[index:]
    #             self.mapToCnode[self.typeC] = pathnode[index:]
    #             for rec in pathlist[index:]:
    #                 rec.label.add("C")
    #                 for enz in rec.getenz():
    #                     enz.labelC.append("C" + str(self.typeC))
    #                     enz.label.add("C")
    #                 for species in rec.getrea():
    #                     if species not in pathnode[index:]:
    #                         species.labelC_side.append(self.typeC)
    #                         species.label.add("C_side")
    #             for species in pathnode[index:]:
    #                 species.labelC_side.append(self.typeC)
    #                 species.label.add("C")
    #             return True
    #     return False

    # def Cycle_Find(self, startRea, pathlist, pathnode):
    #     #print("in DFS")
    #     for downRec in startRea.getDownedge():
    #         if downRec.visited == 0 and downRec.enable:
    #             downRec.visited = 1
    #             pathlist.append(downRec)
    #             pathnode.append(startRea)
    #             downRec.recStack = True
    #             for product in downRec.getpro():
    #                 if product not in self.forbid_node:
    #                     #product.visited = 1
    #                     if self.Cycle_Find(product, pathlist, pathnode):
    #                         return True
    #         else:
    #             if (downRec.recStack == True):
    #                 (self.typeC) += 1
    #                 pathlist.append(downRec)
    #                 return self.backTrace(downRec, pathlist, pathnode)
    #         if downRec.recStack == True:
    #             downRec.recStack = False
    #     return False
     

    # def BFS_findBC(self, label, input_node, notes):
    #     BFS    = [input_node]
    #     while BFS != []:
    #         startPro = BFS.pop(0)
    #         for downRec in startPro.getDownedge():
    #             if downRec.visited == 0:
    #                 downRec.visited = 1
    #                 mark = True

    #                 # check whether this reaction has conflict to previous reactions and input
    #                 # for reactant in downRec.getrea():
    #                 #     if reactant == input_species[1]:
    #                 #         mark = False
    #                 for product in downRec.getpro():
    #                     if product in startPro.pathnode or product in startPro.forbid:
    #                         mark = False
    #                 for enzyme in downRec.getenz():
    #                     if enzyme in startPro.pathnode:
    #                         mark = False
    #                 if mark:
    #                     for product in downRec.getpro():
    #                         if (product.visited == 0 
    #                             and product not in self.forbid_node):
    #                             product.visited = 1
    #                             BFS.append(product)
    #                             for i in startPro.forbid:
    #                                 product.forbid.add(i)
    #                             for i in downRec.getrea():
    #                                 product.forbid.add(i)
    #                             for i in startPro.pathlist:
    #                                 product.pathlist.append(i)
    #                             product.pathlist.append(downRec)
    #                             assert product.pathlist != startPro.pathlist
    #                             for i in startPro.pathnode:
    #                                 product.pathnode.append(i)
    #                             product.pathnode.append(startPro)
    #                             assert product.pathnode != startPro.pathnode
    #                             #if ("A" in product.label):
    #                             if (("B" in product.label) and product not in self.forbid_node
    #                                 and ("C_side" in product.label) and ("C" not in product.label)
    #                                 and (len(self.mapToClist[product.labelC[0]]) == 3)):

    #                                 reaction = edge(None)
    #                                 reactant = node(None)
    #                                 vbuf = []
    #                                 mark_tmp = True
    #                                 for rec in self.mapToClist[product.labelC[0]]:
    #                                     for rea in rec.getrea():
    #                                         if rea in product.forbid:
    #                                             mark_tmp = False
    #                                 for rec in self.mapToClist[product.labelC[0]]:
    #                                     if rec not in product.getDownedge():
    #                                         reaction = rec
    #                                 for rea in reaction.getrea():
    #                                     for rec in product.getDownedge():
    #                                         if ("C" not in rea.label and rea.name != "ethanol"
    #                                             and "C" not in rec.label 
    #                                             and rea in rec.getrea()):
    #                                             mark0 = True
    #                                             for i in rec.getrea():
    #                                                 if i in product.pathnode:
    #                                                     mark0 = False
    #                                             for i in rec.getpro():
    #                                                 if i in product.pathnode:
    #                                                     mark0 = False
    #                                             if mark0:
    #                                                 vbuf.append(rec)
    #                                             reactant = rea
    #                                             rea.tmp_label = ([product.labelB[0],product.labelC[0]])
    #                                 if vbuf != [] and mark_tmp:
    #                                     for i in self.mapToClist[product.labelC[0]]:
    #                                         product.pathlist.append(i)
    #                                     product.pathlist.append(vbuf[0])
    #                                     product.pathnode.append(product)
    #                                     return (True, product, reactant, [label, product.labelB[0],product.labelC[0]])
    #     return (False, node(None), node(None), [])   