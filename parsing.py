def split(inputline):
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
#exodeoxyribonuclease I

inputfile = open("parseddata.txt", encoding="utf-8", mode='r')
fw = open("outputNode.txt", encoding="utf-8", mode='w')
tmpline  = ""
tmpList  = []
HashTable= {}
HashKey  = ""
tempstr  = ""
stack1   = []
stack2   = []
temp     = ""
print("parsing.......")
for line in inputfile:
    tempstr = ""
    HashKey = ""
    tmpList = split(line)
    mark = 0
    #print(tmpList)
    fw.write(tmpList[0] + " ")
    for x in tmpList[1:]:
        if x == "+":
            fw.write(tempstr[:-1] + " ")
            tempstr = ""
        elif x == "=":
            fw.write(tempstr[:-1] + " = ")
            tempstr = ""
        elif x == "{r}":
            fw.write(tempstr[:-1] + " {r}")
            tempstr = ""
        elif x == "{ir}":
            fw.write(tempstr[:-1] + " {ir}")
            tempstr = ""
        else:
            tempstr += x + "_"
    fw.write(tempstr[:-1] + "\n")
    tempstr = ""

inputfile.close()
fw.close()