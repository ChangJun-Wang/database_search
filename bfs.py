import copy

from node import node

from edge import edge

class BFS:
	def __init__(self, startNode, Target_type):
		self.search_limit = 2
		self.bfs 	  	  = [startNode]
		self.type 	  	  = Target_type
		self.candidate 	  = set()

	def search(self):
		self.BFS()
		for species in self.candidate:
			species.record = copy.deepcopy(species.path)
		

	def BFS(self):
		# initializtion for starting node
        initial_node    		 = {}
        initial_node["related"]  = set()
        initial_node["pathlist"] = []
        initial_node["pathenz"]  = []
        initial_node["pathnode"] = []
        self.bfs[0].path.append(initial_node)

        # start to traverse
        count = 0
        while (self.bfs != [] or count >= self.search_limit):
        	count = count + 1
            startPro = self.bfs.pop(0)
            for downRec in startPro.getDownedge():
                # check whether this reaction has conflict to previous reactions and input
                check_downRec = startPro.check_downRec(downRec)
                if check_downRec != []:
                    for product in downRec.getpro():
                        check_pro = product.check_pro(check_downRec)
                        if check_pro != []:  #((product not in self.forbid_node) and (product not in startPro.related)):
                            self.bfs.append(product)
                            product.add_path(startPro, downRec, check_pro)
                            assert product.path != check_pro
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
    		self.candidate.add(product.labelA)
    		return True

    def foundC_side(self, product):
    	if product.path == []:
    		return False
    	for label in product.labelC_side:
    		self.candidate.add(label)
    	return True
 


