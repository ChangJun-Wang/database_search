import copy

from node import node

from edge import edge

class BFS:
    def __init__(self, startNode, Target_type):
        self.search_limit = 1
        self.bfs          = [startNode]
        self.type         = Target_type
        self.candidate    = set()

    # def initial(self, input_species):
    #     candidate = []
    #     for downRec in input_species.getDownedge():
    #         # check whether this reaction has conflict to previous reactions and input
    #         CheckDownRec = input_species.CheckDownRec(downRec)
    #         if CheckDownRec != []:
    #             for product in downRec.getpro():
    #                 CheckProduct = product.CheckProduct(CheckDownRec)
    #                 if CheckProduct != []:
    #                     product.AddPath(startPro, downRec, CheckProduct)
    #                     product.level = startPro.level + 1
    #                     if (product not in self.bfs and product.level <= self.search_limit):
    #                         self.bfs.append(product)

    def search(self):
		# initializtion for starting node
        self.bfs[0].CopyToPath()

        # start to traverse
        count = 0
        while (self.bfs != []):
            startPro = self.bfs.pop(0)
            if (self.type in startPro.label):
                if self.type == "A" and (startPro.getpin() > 2):
                    # print ("find A")
                    self.foundA(startPro)
                elif self.type == "C_side":
                    self.foundC_side(startPro)
                elif self.type == "B":
                    self.foundB(startPro)
                elif self.type == "B_side":
                    self.foundB_side(startPro)
            if startPro.level + 1 <= self.search_limit:
                # for path in startPro.path:
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
            # startPro.level = 0
            # startPro.path = []
        return self.candidate


    def foundA(self, product):
        if product.path == []:
            return False
        else:
            self.candidate.add((product, product.labelA))
            return True

    def foundC_side(self, product):
        for label in product.labelC_side:
            self.candidate.add((product, label))

    def foundB(self, product):
        assert(len(product.getDownedge()) > 1)
        for downRec in product.getDownedge():
            self.candidate.add((product, downRec.name))
            downRec.labelB.append(product.labelB)
            for rea in downRec.getrea():
                if rea != product:
                    rea.label.add("B_side")
                    # rea.labelB_side.append(product.labelB)
                    # rea.labelBmap[product.labelB] = downRec.name

    def foundB_side(self, product):
        for downRec in product.getDownedge():
            if downRec.labelB!=[]:
                self.candidate.add((product, downRec.name))
        # for label in product.labelB_side:
        #     self.candidate.add((product, label))

