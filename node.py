class node :
    def __init__(self, name):
        self.name     = name
        self.enable   = True
        self.upEdge   = []
        self.downEdge = []
        self.catEdge  = []
        self.upEdgePar   = {}
        self.downEdgePar = {}
        self.pin      = 0
        self.pout     = 0
        self.visited  = 0
        self.tmp_vis  = 0
        self.recStack = False
        self.label    = set()
        self.path     = []
        self.pathnode = []
        self.pathlist = []
        self.forbid   = set()
        self.tmp_node = []
        self.tmp_list = []
        self.related  = set()
        self.tmp_forbid   = set()

        self.tmp_label= []
        self.tmp_path = []
        self.labelA   = []
        self.labelB   = []
        self.labelC   = []
        self.labelC_side = []
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

    def setpin(self, pin):
        self.pin = pin
    def setpout(self, pout):
        self.pout = pout

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
        return (self.name)
