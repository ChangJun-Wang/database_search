import sys
 
# the setrecursionlimit function is
# used to modify the default recursion
# limit set by python. Using this,
# we can increase the recursion limit
# to satisfy our needs
 
sys.setrecursionlimit(10**6)

class node :
    def __init__(self, name):
        self.name     = name
        self.upEdge   = []
        self.downEdge = []
        self.catEdge  = []
        self.pin      = 0
        self.pout     = 0
        self.visited  = 0
        self.recStack = False
        self.typeC    = 0
        self.label    = set()
        self.enzyme   = set()
        self.tmp_label= ""
        self.tmp_path = []
        self.labelA   = []
        self.labelB   = []
        self.labelC   = []
        self.labelD   = []
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
        print(self.name)

class edge :
    def __init__(self, name):
        self.name     = name
        self.tmp_path = []
        self.toRea    = []  #
        self.toPro    = []  #
        self.toEnz    = []
        self.label    = set()
        self.visited  = 0
        self.recStack = False
    def getpro(self):
        return self.toPro
    def getrea(self):
        return self.toRea
    def getenz(self):
        return self.toEnz

    def addpro(self, pro):
        (self.toPro).append(pro)
    def addrea(self, rea):
        (self.toRea).append(rea)
    def addenz(self, enz):
        (self.toEnz).append(enz)

    def show(self):
        temp = "\n" + "reaction : "
        temp += (str(self.name))
        temp += "\n" + ("Enzyme : ") + "\n"
        for i in self.getenz():
            temp +=  (str(i.name)) + "  "
        temp += "\n" + ("Reactant : ") + "\n"
        for i in self.getrea():
            temp +=  (str(i.name)) + "  "
        temp += "\n" + ("Product : ") + "\n"
        for i in self.getpro():
            temp += (str(i.name)) + "  "
        return temp


class search:
    def __init__(self):
        self.tmp_path = {}
        self.mapToNode   = {}
        self.mapToEdge   = {}
        self.mapToCycle  = {}
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
        self.store       = []

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
            if (self.mapToNode[nodes].getpin() > 1):
                self.typeA += 1
                self.mapToNode[nodes].label.add("A")
                self.mapToNode[nodes].labelA.append("A" + str(self.typeA))
            if (self.mapToNode[nodes].getpout() > 1):
                self.typeB += 1
                self.mapToNode[nodes].label.add("B")
                self.mapToNode[nodes].labelB.append("B" + str(self.typeB))
            if (self.mapToNode[nodes].getpout() > 0) and (self.mapToNode[nodes].getpin() > 0):
                self.typeD += 1
                self.mapToNode[nodes].label.add("D")
                self.mapToNode[nodes].labelD.append("D" + str(self.typeD))
            if (len(self.mapToNode[nodes].getCatedge()) > 1):
                self.typeE += 1
                self.mapToNode[nodes].label.add("E")
                self.mapToNode[nodes].labelD.append("E" + str(self.typeE))
            if self.mapToNode[nodes].visited == 0:
                self.mapToNode[nodes].visited = 1
                self.DFS(self.mapToNode[nodes], [], [])
        self.ClearVis()

    def parsing(self, inputfile):
        temp = ""
        tmpSpecies = ""
        tmpList = []
        enzyme = ""
        is_product = 0
        HashKey = "" 

        print("in class parsing")
        for line in inputfile:
            tmpList = self.split(line)

            if (tmpList[0] == "RN"):
                enzyme = tmpList[1]
                if (enzyme not in self.mapToNode):
                    self.mapToNode[enzyme] = node(enzyme)
                    self.nodeList.append(enzyme)
            else:
                is_product = 0
                self.reaction += 1
                self.mapToEdge[self.reaction] = edge(self.reaction)
                #self.edgeList.append(self.reaction)

                self.mapToEdge[self.reaction].addenz(self.mapToNode[enzyme])
                self.mapToNode[enzyme].addCatedge(self.mapToEdge[self.reaction])

                for species in tmpList[1:]:
                    HashKey = HashKey + species[0]
                    if (species == "="):
                        is_product = 1
                    elif(species == "{r}"):
                        pass
                    elif(species == "{ir}" or species == "?" or species == "more"):
                        pass
                    else:
                        if (species not in self.mapToNode):
                            self.mapToNode[species] = node(species)
                            self.nodeList.append(species)
                        if (is_product == 0):
                            self.mapToEdge[self.reaction].addrea(self.mapToNode[species])
                            self.mapToNode[species].addDownedge(self.mapToEdge[self.reaction])
                            self.mapToNode[species].addpout()
                        else:
                            self.mapToEdge[self.reaction].addpro(self.mapToNode[species])
                            self.mapToNode[species].addUpedge(self.mapToEdge[self.reaction])
                            if (enzyme != 'spontaneous_reaction' and enzyme not in self.mapToNode[species].enzyme):
                                self.mapToNode[species].enzyme.add(enzyme)
                                self.mapToNode[species].addpin()


                if (HashKey in self.HashTable):
                    typeD_tmp = 0
                    for rec in self.HashTable[HashKey]:
                        if ((self.mapToEdge[rec].getrea() == self.mapToEdge[self.reaction].getrea()) 
                            and (self.mapToEdge[rec].getpro() == self.mapToEdge[self.reaction].getpro())
                            and (self.mapToNode[enzyme] not in self.mapToEdge[rec].getenz())):
                            typeD_tmp += 1
                            ((self.mapToEdge[rec]).label).add("D")
                            for i in self.mapToEdge[rec].getpro():
                                i.label.add("D")
                            if (self.reaction not in self.HashCollect):
                                templist = []
                                templist.append(rec)
                                self.HashCollect[self.reaction] = templist
                            else:
                                self.HashCollect[self.reaction].append(rec)

                    if (typeD_tmp == 0):
                        self.HashTable[HashKey].append(self.reaction)
                    elif (typeD_tmp == 1):
                        self.typeD += 1
                    else:
                        pass
                        #print("bug")
                        #print("duplicate")
                    #HashTable[HashKey] = reaction
                else:
                    self.HashTable[HashKey] = []
                    self.HashTable[HashKey].append(self.reaction)
                HashKey = ""
        self.initialize()


    def backTrace(self, startRec, pathlist, pathnode):
        #print("in backtrae")
        for index in range(len(pathlist)):
            if startRec == pathlist[index]:
                #print("add label")
                #(self.typeC) += 1
                #startRec.label.add("C" + str(self.typeC))
                #startRec.label.add("C")
                #for enz in startRec.getenz():
                #    (enz).label.add("C" + str(self.typeC))
                #    (enz).label.add("C")
                self.mapToCycle["C" + str(self.typeC)] = pathlist[index:]
                for rec in pathlist[index:]:
                    rec.label.add("C")
                    for species in rec.getpro():
                        if species not in pathnode[index:]:
                            species.labelC.append("C" + str(self.typeC))
                            species.label.add("C_side")
                for species in pathnode[index:]:
                    species.labelC.append("C" + str(self.typeC))
                    species.label.add("C")
                return True
        return False

    def DFS(self, startRea, pathlist, pathnode):
        #print("in DFS")
        tmp_node = pathnode
        tmp_node.append(startRea)
        for downRec in startRea.getDownedge():
            if downRec.visited == 0:
                downRec.visited = 1
                tmp_rec = pathlist
                tmp_rec.append(downRec)
                downRec.recStack = True
                for product in downRec.getpro():
                    if product.name != "H2O":
                        if self.DFS(product, tmp_rec, tmp_node):
                            return True
            else:
                if (downRec.recStack == True):# and (self.count > 2):
                    (self.typeC) += 1
                    pathlist.append(downRec)
                    return self.backTrace(downRec, pathlist, pathnode)
            if downRec.recStack == True:
                downRec.recStack = False
        return False
        
    def reveal(self):
        for rec in range(self.reaction):
            if self.mapToEdge[rec+1].recStack:
                print(self.mapToEdge[rec+1].show())


    def DFS_all(self, startNode, notes, pathlist, pathnode):
        #print("in DFS")
        for downRec in startNode.getDownedge():
            mark = True
            for reactant in downRec.getrea():
                if reactant.name == "D-glucose":
                    mark = False
            for product in downRec.getpro():
                if product in pathnode:
                    mark = False
            if mark:
                for product in downRec.getpro():
                    if product.visited == 0 and product.name != "H2O":
                        product.visited = 1
                        tmp_rec = pathlist
                        tmp_rec.append(downRec)
                        tmp_node = pathnode
                        tmp_node.append(startNode)
                        if (notes in product.label):
                            label = product.labelA[0]
                            product.tmp_label = label
                            self.tmp_path[label]   = tmp_rec
                        self.DFS_all(product, notes, tmp_rec, tmp_node)
        #return False

    def DFS_down(self, startNode, notes, pathlist, pathnode):
        #print("DFS down")
        for downRec in startNode.getDownedge():
            #print("DFS down")
            mark = True
            for reactant in downRec.getrea():
                if reactant.name == "acetone":
                    mark = False
            for product in downRec.getpro():
                if product in pathnode:
                    mark = False
            if mark:
                for product in downRec.getpro():
                    if(product.visited == 0):
                        #print("find", notes)
                        product.visited = 1
                        #product.recStack = True
                        #downRec.recStack = True
                        pathlist.append(downRec)
                        pathnode.append(startNode)
                        if ((product.tmp_label != "") and (product.name != "D-glucono-1,5-lactone") and (product.name != "H2O") and (product.name != "H+") and (product.tmp_label[0] == notes)): #and (len(self.tmp_path[product.tmp_label]) < 10)):
                            print("find0", notes)
                            print("\n********** path0 *********\n")
                            path = product.tmp_path[0]
                            while path != self.mapToNode["acetone"]:
                                print(path.show())
                                path = path.tmp_path[0]
                            print("\n********** path0 *********")

                            print("find1", notes)
                            print("\n********** path1 *********\n")                            
                            for path in pathlist:
                                print(path.show())
                            #for node in pathnode:
                            print(product.name)
                            print("\n********** path1 *********")
                            if len(product.getUpedge()) > 2 and notes == "A":
                                count = 0
                                print("////////other option///////////")
                                for x in product.getUpedge():
                                    if count > 3:
                                        break
                                    else:
                                        if self.mapToNode["alcohol_dehydrogenase"] not in x.getenz():
                                            print(x.show())
                                            count += 1
                                print("////////other option end///////////")
                            pathnode = []
                            pathnode.append(product)
                            return (True, pathlist, pathnode)
                        (r, pl, pn)  = self.DFS_down(product, notes, pathlist, pathnode)
                        if r:
                            return (r, pl, pn)
                    #product.recStack = False
                    #downRec.recStack = False
        return (False, [], [])


    def DFS_up(self, startPro, notes, pathlist, pathnode):
        print("in DFS up")
        mark = True
        for upRec in startPro.getUpedge():
            for reactant in upRec.getrea():
                if(reactant.visited == 0 and reactant.name != "H2O" and mark):
                    reactant.visited = 1
                    reactant.recStack = True
                    upRec.recStack = True
                    pathlist.append(upRec)
                    pathnode.append(startPro)
                    if (notes in (reactant).label) and mark:
                        print("find", notes)
                        (reactant).label.remove(notes)
                        for path in pathlist:
                            print(path.show())
                        return (True, pathnode, pathlist)
                    return self.DFS_up(reactant, notes, pathlist, pathnode)
                reactant.recStack = False
                upRec.recStack = False
        return (False, [], [])

    def ClearVis(self):
        for nodes in self.nodeList:
            self.mapToNode[nodes].visited  = 0
            self.mapToNode[nodes].recStack = False
        for rec in range(self.reaction):
            self.mapToEdge[rec+1].visited  = 0
            self.mapToEdge[rec+1].recStack = False

    def show(self, node):
        if "C" in node.label:
            for label in node.labelC:
                if ((len(self.mapToCycle[label]) > 2) 
                    and (len(self.mapToCycle[label]) < 5)):
                    if self.mapToCycle[label][0].getenz() != self.mapToCycle[label][1].getenz() :
                        print("There are the related reaction: \n")
                        for i in self.mapToCycle[label]:
                            print(label)
                            print(i.show())
                        return True
        return False
    
    def showRec(self, node):
        if "C" in node.label:
            for label in node.labelC:
                if ((len(self.mapToCycle[label]) == 3)):
                    #print("There are the related reaction: \n")
                    pathlist = self.mapToCycle[label]
                    for index in range(len(pathlist) - 1):
                        for product in pathlist[index].getpro():
                            if (product in pathlist[index + 1].getrea()):
                                product.label.add(label)
                                product.label.add("C")
                                break
                    return pathlist
        return []

    def BFS_down(self, startPro, notes, pathlist, pathnode):
        #record = [startPro]
        BFS    = [startPro]
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                if downRec.visited == 0:
                    downRec.visited = 1
                    downRec.tmp_path = []
                    downRec.tmp_path.append(startPro)
                    mark = True
                    for reactant in downRec.getrea():
                        if reactant.name == "D-glucose" or reactant.name == "glucose":
                            mark = False
                    for product in downRec.getpro():
                        if product in pathnode:
                            mark = False
                    if mark:
                        for product in downRec.getpro():
                            if (product.visited == 0 and product.name != "H2O"):
                                product.visited = 1
                                product.tmp_path = []
                                product.tmp_path.append(downRec)
                                BFS.append(product)
                                pathlist.append(downRec)
                                pathnode.append(startPro)
                                if (notes in product.label):
                                    label = product.labelA[0]
                                    product.tmp_label       = label

    def BFS_downC(self, startPro, notes, pathlist, pathnode):
        #record = [startPro]
        #print("find C_side")
        #print(startPro)
        BFS    = [startPro]
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                if downRec.visited == 0:
                    downRec.visited = 1
                    downRec.tmp_path = []
                    downRec.tmp_path.append(startPro)
                    mark = True
                    for reactant in downRec.getrea():
                        if reactant.name == "D-glucose" or reactant.name == "glucose":
                            mark = False
                    for product in downRec.getpro():
                        if product in pathnode:
                            mark = False
                    if mark:
                        for product in downRec.getpro():
                            if (product.visited == 0 and product.name != "H2O"):
                                product.visited = 1
                                product.tmp_path = []
                                product.tmp_path.append(downRec)
                                BFS.append(product)
                                pathlist.append(downRec)
                                pathnode.append(startPro)
                                if (notes in product.label and len(product.getDownedge()) > 2):
                                    #print("find C_side")
                                    label = product.labelC[0]
                                    #print(label)
                                    if product not in self.mapToCycle[label][0].getpro():
                                        for x in self.mapToCycle[label][0].getpro():
                                            if "C" not in x.label:
                                                x.tmp_label = label
                                                x.tmp_path = []
                                                x.tmp_path.append(self.mapToCycle[label][0])
                                                
                                    else:
                                        for x in self.mapToCycle[label][0].getrea():
                                            if "C" not in x.label:
                                                x.tmp_label = label
                                                x.tmp_path = []
                                                x.tmp_path.append(self.mapToCycle[label][0])
                                                


    def BFS_downC1(self, startPro, notes, pathlist, pathnode):
        BFS        = [startPro]
        while BFS != []:
            startPro = BFS.pop(0)
            for downRec in startPro.getDownedge():
                if downRec.visited == 0:
                    downRec.visited = 1
                    downRec.tmp_path = []
                    downRec.tmp_path.append(startPro)
                    mark = True
                    for reactant in downRec.getrea():
                        if reactant.name == "D-glucose" or reactant.name == "glucose":
                            mark = False
                    for product in downRec.getpro():
                        if product in pathnode:
                            mark = False
                    if mark:
                        for product in downRec.getpro():
                            if(product.visited == 0):
                                #print("find", notes)
                                product.visited = 1
                                #product.recStack = True
                                #downRec.recStack = True
                                pathlist.append(downRec)
                                #pathnode.append(startNode)
                                if ((product.tmp_label != "") and (product.name != "D-glucono-1,5-lactone") and (product.name != "H2O") and (product.name != "H+") and (product.tmp_label[0] == notes)): #and (len(self.tmp_path[product.tmp_label]) < 10)):
                                    print("find0", notes)
                                    print("\n********** pathC0 *********\n")
                                    path = product.tmp_path[0]
                                    #while path != self.mapToNode["NADH"]:
                                    print(path.show())
                                    #    path = path.tmp_path[0]
                                    print("\n********** pathC0 *********")

                                    print("find1", notes)
                                    print("\n********** pathC1 *********\n")                            
                                    #for path in pathlist:
                                    print(downRec.show())
                                    #for node in pathnode:
                                    print(product.name)
                                    print("\n********** pathC1 *********")
                                    pathnode = []
                                    pathnode.append(product)
                                    return True
        return False


    def DFS_downC(self, startNode, notes, pathlist, pathnode):
        #print("DFS downC")
        #startNode.show()
        for downRec in startNode.getDownedge():
            #print("DFS down")
            mark = True
            for reactant in downRec.getrea():
                if reactant.name == "acetone" and reactant.name == "D-glucose":
                    mark = False
            for product in downRec.getpro():
                if product in pathnode:
                    mark = False
            if mark:
                #print("DFS down")
                for product in downRec.getpro():
                    if(product.visited == 0):
                        #print("find", notes)
                        product.visited = 1
                        #product.recStack = True
                        #downRec.recStack = True
                        pathlist.append(downRec)
                        #pathnode.append(startNode)
                        if ((product.tmp_label != "") and (product.name != "D-glucono-1,5-lactone") and (product.name != "H2O") and (product.name != "H+") and (product.tmp_label[0] == notes)): #and (len(self.tmp_path[product.tmp_label]) < 10)):
                            print("find0", notes)
                            print("\n********** pathC0 *********\n")
                            path = product.tmp_path[0]
                            #while path != self.mapToNode["NADH"]:
                            print(path.show())
                            #    path = path.tmp_path[0]
                            print("\n********** pathC0 *********")

                            print("find1", notes)
                            print("\n********** pathC1 *********\n")                            
                            #for path in pathlist:
                            print(downRec.show())
                            #for node in pathnode:
                            print(product.name)
                            print("\n********** pathC1 *********")
                            pathnode = []
                            pathnode.append(product)
                            return True
                        if self.DFS_downC(product, notes, pathlist, pathnode):
                            return True
                    #product.recStack = False
                    #downRec.recStack = False
        return False


if __name__ == '__main__':
    #traversal(store)

    searcher = search()
    f = open("parseddata.txt", encoding="utf-8", mode='r')
    searcher.parsing(f)
    f.close()
    print("num of reaction : ", searcher.reaction)
    print("num of species : ", len(searcher.nodeList))
    
    temp = 0
    fw = open("outputNode.txt", encoding="utf-8", mode='w')


    print("typeA : ", searcher.typeA)
    print("typeB : ", searcher.typeB)
    print("typeC : ", searcher.typeC)
    print("typeD : ", searcher.typeD)

    for nodes in searcher.nodeList:
        fw.write(searcher.mapToNode[nodes].name + " : ")
        for x in searcher.mapToNode[nodes].label:
            fw.write(x + ",")
        fw.write("\n")
    fw.close()


    count = 0
    notes = ["C","C","C"]
    print ("\n", "*************************************************")

    """for nodes in searcher.nodeList:
        if len(notes) == 0:
            break
        if ("C" in searcher.mapToNode[nodes].label and searcher.mapToNode[nodes].name != "H2O"):
            temp = searcher.showRec(searcher.mapToNode[nodes])
            if temp != []:
                for rea0 in temp[0].getrea():
                    for rea1 in temp[1].getrea():
                        if ("C" not in rea0.label) and ("C" not in rea1.label):
                            for downRec in rea0.getDownedge():
                                if len(notes) == 0:
                                    break
                                if (downRec in rea1.getDownedge()):
                                    notes.pop()
                                    print(downRec.show())
                                    print("\n", "find C", "\n")
                                    for rec in temp:
                                        print(rec.show())
                    
                                    count = 0
                                    for rec in rea1.getUpedge():
                                        if count == 2:
                                            break
                                        if (rea0 not in rec.getpro()) and (rea0 not in rec.getrea()) and (searcher.mapToNode["alcohol_dehydrogenase"] not in rec.getenz()):
                                            print(rec.show())
                                            print("find A rea1", "\n")
                                            count += 1
                                    count = 0
                                    enzyme = "alcohol_dehydrogenase"
                                    for rec in rea0.getUpedge():
                                        if count == 7:
                                            break
                                        if (rea1 not in rec.getpro()) and (rea1 not in rec.getrea()) and (searcher.mapToNode[enzyme] not in rec.getenz()):
                                            enzyme = rec.getenz()[0].name
                                            print(rec.show())
                                            print("find A rea0", "\n")
                                            count += 1
                                    print("\n", "*************************************************")"""

    #specified input species :
    input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    notes         = ["A","A","E","C"]

    #searcher.DFS_all(input_species[0], "A", [], [])
    searcher.BFS_down(input_species[0], "A", [], [])
    searcher.ClearVis()
    (result0, pathlist0, pathnode0) = searcher.DFS_down(input_species[1], "A", [], [])
    print(result0)

    print(pathlist0)
    
    (result1, pathlist1, pathnode1) = searcher.DFS_down(input_species[1], "A", [], [])
    print(result1)

    print(pathlist1)

    pathnode0[0].show()
    pathnode1[0].show()

    searcher.BFS_downC(pathnode0[0], "C_side", [], [])
    searcher.ClearVis()
    result3 = searcher.BFS_downC1(pathnode1[0], "C", [], [])

    print("result3 : ", result3)

    print ("count : ", count)

