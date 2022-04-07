class path :
    def __init__(self):
        self.pathnode     = []
        self.pathlist     = []
        self.related      = []
        self.pathenz      = []

    def length(self):
        return len(self.pathlist)
        
    def getNode(self):
        return self.pathnode
    def getList(self):
        return self.pathlist
    def getRelate(self):
        return self.related
    def getEnz(self):
        return self.pathenz