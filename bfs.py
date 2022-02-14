import copy

from node import node

from edge import edge

class BFS:
    def __init__(self, startNode, Target_type, mapToCnode, mapToClist):
        self.search_limit = 1
        self.bfs          = [startNode]
        self.type         = Target_type
        self.candidate    = set()
        self.mapToCnode   = mapToCnode
        self.mapToClist   = mapToClist

    # def search(self):
    #     self.BuildStart(input_species)
    #     input_species[0].CopyToRecord("A", 0)
    #     input_species[1].CopyToRecord("A", 0)
    #     self.BFS()
    #     return self.candidate
	
    def BuildStart(self, input_species):
        (input_species[0].path[0])["pathnode"].append(input_species[0])
        (input_species[0].path[0])["pathnode"].append(input_species[1])
        (input_species[1].path[0])["pathnode"].append(input_species[0])
        (input_species[1].path[0])["pathnode"].append(input_species[1])
        # for startPro in input_species:
        #     for downRec in startPro.getDownedge():
        #         # check whether this reaction has conflict to previous reactions and input
        #         CheckDownRec = startPro.CheckDownRec(downRec)
        #         if CheckDownRec != []:
        #             for product in downRec.getpro():
        #                 CheckProduct = product.CheckProduct(CheckDownRec)
        #                 if CheckProduct != []:  # ((product not in self.forbid_node) and (product not in startPro.related)):
        #                     product.AddPath(startPro, downRec, CheckProduct)
        #                     product.level = startPro.level + 1
        #                     if (product not in self.bfs and product.level < self.search_limit):
        #                         self.bfs.append(product)
        #                     assert product.path != CheckProduct
        #                     assert len(product.path) != 0
        input_species[0].CopyToRecord("A", 0)
        input_species[1].CopyToRecord("A", 0)

    def search(self):
		# initializtion for starting node
        self.bfs[0].CopyToPath()

        # start to traverse
        count = 0
        while (self.bfs != []): # and count < self.search_limit):
            # print ("level :", count)
            # count = count + 1
            startPro = self.bfs.pop(0)
            # print ("bfslen :", len(self.bfs))
            if (self.type in startPro.label): # and (startPro.getpin() > 2):
                assert self.type == "A" or self.type == "C_side"
                if self.type == "A" and (startPro.getpin() > 2):
                    # print ("find A")
                    self.foundA(startPro)
                elif self.type == "C_side":
                    self.foundC_side(startPro)
                # elif self.type == "B":
                #   pass
                # elif self.type == "C":
                #   pass
                # elif self.type == "D":
                #   pass
            if startPro.level + 1 <= self.search_limit:
                for downRec in startPro.getDownedge():
                    # check whether this reaction has conflict to previous reactions and input
                    CheckDownRec = startPro.CheckDownRec(downRec)
                    if CheckDownRec != []:
                        for product in downRec.getpro():
                            CheckProduct = product.CheckProduct(CheckDownRec)
                            if CheckProduct != []:  # ((product not in self.forbid_node) and (product not in startPro.related)):
                                product.AddPath(startPro, downRec, CheckProduct)
                                product.level = startPro.level + 1
                                if (product not in self.bfs and product.level <= self.search_limit):
                                    self.bfs.append(product)
                                assert product.path != CheckProduct
                                assert len(product.path) != 0
        return self.candidate

    def foundA(self, product):
    	if product.path == []:
    		return False
    	else:
    		self.candidate.add((product, product.labelA))
    		return True

    def foundC_side(self, product):
    	if product.path == []:
    		return False
    	for label in product.labelC_side:
    		self.candidate.add((product, label))
    	return True
 

    def CheckLabelC(self, input_node):
        for c_side in product.labelC_side:
            for path in product.path:
                if self.CheckCycle(path, c_side):
                    product.AddPath(startPro, downRec, CheckProduct)
                    return (True, product, c_side)
        return (False, node(None), None)

    def CheckCycle(self, path, c_side):
        for node in self.mapToCnode[c_side]:
            if node in path["related"]:
                return False
        return True

