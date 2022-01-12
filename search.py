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
        self.pin      = 0
        self.pout     = 0
        self.visited  = 0
        self.recStack = False
        self.typeC = 0
        self.label = []
    def getUpedge(self):
        return self.upEdge
    def getDownedge(self):
        return self.downEdge
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

class edge :
    def __init__(self, name):
        self.name     = name
        self.toRea    = []  #
        self.toPro    = []  #
        self.toEnz    = []
        self.label = []
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
        self.mapToNode   = {}
        self.mapToEdge   = {}
        self.mapToCycle  = {}
        self.nodeList    = []
        self.edgeList    = []
        self.HashTable   = {}
        self.HashCollect = {}
        self.reaction    = 0
        self.typeC       = 0
        self.typeD       = 0
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
                self.mapToNode[enzyme].addUpedge(self.mapToEdge[self.reaction])
                self.mapToNode[enzyme].addDownedge(self.mapToEdge[self.reaction])

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
                            self.mapToNode[species].addpin()

                if (HashKey in self.HashTable):
                    typeD_tmp = 0
                    for rec in self.HashTable[HashKey]:
                        if ((self.mapToEdge[rec].getrea() == self.mapToEdge[self.reaction].getrea()) 
                            and (self.mapToEdge[rec].getpro() == self.mapToEdge[self.reaction].getpro())
                            and (self.mapToNode[enzyme] not in self.mapToEdge[rec].getenz())):
                            typeD_tmp += 1
                            ((self.mapToEdge[rec]).label).append("D")
                            for i in self.mapToEdge[rec].getpro():
                                i.label.append("D")
                            #for j in self.mapToEdge[rec].getrea():
                                #j.label.append("D")
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

    def backTrace(self, startRec, pathlist):
        #print("in backtrae")
        for index in range(len(pathlist)):
            if startRec == pathlist[index]:
                #print("add label")
                #(self.typeC) += 1
                startRec.label.append("C" + str(self.typeC))
                startRec.label.append("C")
                for enz in startRec.getenz():
                    (enz).label.append("C" + str(self.typeC))
                    (enz).label.append("C")
                self.mapToCycle["C" + str(self.typeC)] = pathlist[index:]
                return True
        return False

    def DFS(self, startRea, pathlist):
        #print("in DFS")
        for downRec in startRea.getDownedge():
            if downRec.visited == 0:
                downRec.visited = 1
                tmp = pathlist
                tmp.append(downRec)
                downRec.recStack = True
                #self.count += 1
                for product in downRec.getpro():
                    if product.name != "H2O":
                        if self.DFS(product, tmp):
                            return True
            else:
                if (downRec.recStack == True):# and (self.count > 2):
                    (self.typeC) += 1
                    #print("A cycle")
                    pathlist.append(downRec)
                    return self.backTrace(downRec, pathlist)
                    #self.count = 0
            if downRec.recStack == True:
                downRec.recStack = False
                #self.count -= 1
        return False
        
    def reveal(self):
        for rec in range(self.reaction):
            if self.mapToEdge[rec+1].recStack:
                print(self.mapToEdge[rec+1].show())

    def DFSsearch(self, startPro, notes, pathlist):
        for downRec in startPro.getDownedge():
            for product in downRec.getpro():
                if(product.visited == 0 and product.name != "H2O"):
                    product.visited = 1
                    product.recStack = True
                    downRec.recStack = True
                    pathlist.append(downRec)
                    for enz in downRec.getenz():
                        if ("C" in (enz).label) and self.show(enz):
                            print("find C")
                            for path in pathlist:
                                print(path.show())
                            #self.reveal()
                            #for rec in product.getUpedge():
                            #    if "C" in rec.label:
                            #        print(rec.show())
                            #        print ("Find type ", "note")
                            return True
                    return self.DFSsearch(product, notes, pathlist)
                product.recStack = False
                downRec.recStack = False
        return False

    def ClearVis(self):
        for nodes in self.nodeList:
            self.mapToNode[nodes].visited  = 0
            self.mapToNode[nodes].recStack = False
        for rec in range(self.reaction):
            self.mapToEdge[rec+1].recStack = False

    def show(self, node):
        for label in node.label:
            if len(label) > 1:
                if (len(self.mapToCycle[label]) > 2) and (len(self.mapToCycle[label]) < 5):
                    if self.mapToCycle[label][0].getenz() != self.mapToCycle[label][1].getenz() :
                        print("There are the related reaction: \n")
                        for i in self.mapToCycle[label]:
                            print(label)
                            print(i.show())
                        return True
        return False
    
    def showRec(self, edge):
        for label in edge.label:
            if len(label) > 1:
                #if (len(self.mapToCycle[label]) > 2) and (len(self.mapToCycle[label]) < 10):
                #if self.mapToCycle[label][0].getenz() != self.mapToCycle[label][1].getenz() :
                print("There are the related reaction: \n")
                for i in self.mapToCycle[label]:
                    print(i.show())
                return True


if __name__ == '__main__':
    #traversal(store)
    typeA = 0
    typeB = 0
    typeC = 0

    searcher = search()
    f = open("parseddata.txt", encoding="utf-8", mode='r')
    searcher.parsing(f)
    f.close()
    print("num of reaction : ", searcher.reaction)
    print("num of species : ", len(searcher.nodeList))
    
    temp = 0
    fw = open("outputNode.txt", encoding="utf-8", mode='w')
    for nodes in searcher.nodeList:
        if (searcher.mapToNode[nodes].getpin() > 1):
            typeA += 1
            searcher.mapToNode[nodes].label.append("A")
        if (searcher.mapToNode[nodes].getpout() > 1):
            typeB += 1
            searcher.mapToNode[nodes].label.append("B")
        if searcher.mapToNode[nodes].visited == 0:
            searcher.mapToNode[nodes].visited = 1
            searcher.DFS(searcher.mapToNode[nodes], [])

    #for edges in range(searcher.reaction):
    #    if (searcher.mapToEdge[edges+1].visited == 0):
    #        searcher.mapToEdge[edges+1].visited = 1
    #        searcher.mapToEdge[edges+1].recStack = True
    #        for product in searcher.mapToEdge[edges+1].getpro():
     #           if product.name != "H2O":
     #               searcher.DFS(searcher.mapToNode[nodes], [])
     #       searcher.mapToEdge[edges+1].recStack = False
            #for edges in range(searcher.reaction):
            #    if searcher.mapToEdge[edges+1].recStack == True:
            #        print("still true stack exists")
    #for nodes in searcher.nodeList:
    #    typeC += searcher.mapToNode[nodes].typeC


    print("typeA : ", typeA)
    print("typeB : ", typeB)
    print("typeC : ", searcher.typeC)
    print("typeD : ", searcher.typeD)

    for nodes in searcher.nodeList:
        fw.write(searcher.mapToNode[nodes].name + " : ")
        for x in searcher.mapToNode[nodes].label:
            fw.write(x + ",")
        fw.write("\n")
    fw.close()

    searcher.ClearVis()

    #fc = open("cycle.txt", encoding="utf-8", mode='w')
    #for i in range(searcher.typeC):
    #    fc.write(str(i+1))
    #    fc.write(searcher.mapToCycle["C" + str(i+1)])
    #    fc.write("\n")
    #fc.close()

    count = 0
    notes = ["C"]
    print ("\n", "*************************************************")
    '''for edges in range(searcher.reaction):
        if "C" in searcher.mapToEdge[edges+1].label: # and searcher.mapToNode[nodes].name != "H2O":
            if searcher.show(searcher.mapToEdge[edges+1]):
                print("find C")
                break'''
    for nodes in searcher.nodeList:
        if "A" in searcher.mapToNode[nodes].label and searcher.mapToNode[nodes].name != "H2O":
            result = searcher.DFSsearch(searcher.mapToNode[nodes], notes, [])
            if result:
                for rec in searcher.mapToNode[nodes].getUpedge():
                    print(rec.show())
                    print("find A")
                notes.pop()
                print ("*************************************************")
                count += 1

            if not notes:
                break

    print ("count : ", count)