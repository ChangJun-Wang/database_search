import copy

from node import node

from edge import edge

class BFS:
	def __init__(self, startNode, Target_type, mapToCnode, mapToClist):
		self.search_limit = 2
		self.bfs 	  	  = [startNode]
		self.type 	  	  = Target_type
		self.candidate 	  = set()
        self.mapToCnode   = mapToCnode
        self.mapToClist   = mapToClist

	def search(self):
		self.BFS()
        return self.candidate
		
	def BFS(self):
		# initializtion for starting node
        self.bfs[0].CopyToPath()

        # start to traverse
        count = 0
        while (self.bfs != [] or count >= self.search_limit):
        	count = count + 1
            startPro = self.bfs.pop(0)
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                CheckDownRec = startPro.CheckDownRec(downRec)
                if CheckDownRec != []:
                    for product in downRec.getpro():
                        CheckProduct = product.CheckProduct(CheckDownRec)
                        if CheckProduct != []:  #((product not in self.forbid_node) and (product not in startPro.related)):
                            self.bfs.append(product)
                            product.AddPath(startPro, downRec, CheckProduct)
                            assert product.path != CheckProduct
                            assert len(product.path) != 0
                            if (self.type in product.label): # and (product.getpin() > 2):
                            	assert self.type == "A" or self.type == "C_side"
                            	if self.type == "A":
                            		self.foundA(product)
                            	elif self.type == "C_side":
                            		self.foundC_side(product)
                            	# elif self.type == "B":
                            	# 	pass
                            	# elif self.type == "C":
                            	# 	pass
                            	# elif self.type == "D":
                            	# 	pass

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

