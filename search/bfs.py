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
	
    def BuildStart(self, input_species):
        (input_species[0].path[0])["pathnode"].append(input_species[1])
        (input_species[1].path[0])["pathnode"].append(input_species[0])

    def search(self):
		# initializtion for starting node
        self.bfs[0].CopyToPath()

        # start to traverse
        count = 0
        while (self.bfs != []):
            startPro = self.bfs.pop(0)
            if (self.type in startPro.label):
                assert self.type == "A" or self.type == "C_side"
                if self.type == "A" and (startPro.getpin() > 2):
                    # print ("find A")
                    self.foundA(startPro)
                elif self.type == "C_side":
                    self.foundC_side(startPro)

            if startPro.level + 1 <= self.search_limit:
                for downRec in startPro.getDownedge():
                    # check whether this reaction has conflict to previous reactions and input
                    CheckDownRec = startPro.CheckDownRec(downRec)
                    if CheckDownRec != []:
                        for product in downRec.getpro():
                            CheckProduct = product.CheckProduct(CheckDownRec)
                            if CheckProduct != []:
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
