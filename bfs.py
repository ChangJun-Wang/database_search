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
		# for species in self.candidate:
		# 	species.record = copy.deepcopy(species.path)
        return self.candidate
		

	def BFS(self):
		# initializtion for starting node
            # initial_node    		 = {}
            # initial_node["related"]  = set()
            # initial_node["pathlist"] = []
            # initial_node["pathenz"]  = []
            # initial_node["pathnode"] = []
        # if self.type == "A":
        self.bfs[0].CopyToPath()
        # elif self.type == "C_side":
        #     self.bfs[0].CopyToPath("A", self.bfs[0].labelA)
            # self.bfs[0].path.append(initial_node)

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
    		self.candidate.add((product, product.labelA))
    		return True

    def foundC_side(self, product):
    	if product.path == []:
    		return False
    	for label in product.labelC_side:
    		self.candidate.add((product, label))
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

