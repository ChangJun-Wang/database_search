def split(inputline):
    temp    = ""
    tmpList = []
    for w in inputline:
        if (w == " " or w == "\n"):
            if temp != "":
                if len(temp) > 3:
                    if temp[2] == "_":
                        if temp[0:2].isdigit():
                            temp = temp[3:]
                        elif temp[0:2] == "an":
                            temp = temp[3:]
                if len(temp) > 2:
                    if temp[1] == "_":
                        if temp[0].isdigit():
                            temp = temp[2:]
                        elif temp[0] == "a":
                            temp = temp[2:]
                tmpList.append(temp)
            temp = ""
        else:
            temp += w
    return tmpList
#exodeoxyribonuclease I

def checkRea(list0, list1):
    if len(list0) != len(list1):
        return False
    c = 0
    for i in list0:
        if i in list1:
            c += 1
    if c == len(list0):
        return True
    return False

# def checkPro():
#     pass

def check(table, species):
    for i in table:
        if checkRea(species[0], i[0]) and checkRea(species[1], i[1]):
            return False
    return True

inputfile = open("parseddata_2.txt", encoding="utf-8", mode='r')
fw = open("outputNode.txt", encoding="utf-8", mode='w')
print("parsing.......")
table = []
count = 0
for line in inputfile:
    reactant = []
    product  = []
    is_product = False

    if line[0:2] == "RN":
        table = []
        fw.write(line)
    else:
        tmpList = split(line)
        tempstr = ""
        for x in tmpList:
            if x == "=":
                is_product = True
            if is_product:
                product.append(x)
            else:
                reactant.append(x)
        if check(table, (reactant, product)):
            table.append((reactant, product))
            fw.write(str(count) + " ")
            fw.write(line)
            count += 1

####################################################################
# def split(inputline):
#     temp    = ""
#     tmpList = []
#     for w in inputline:
#         if (w == " " or w == "\n"):
#             if temp != "":
#                 if len(temp) > 3:
#                     if temp[2] == "_":
#                         if temp[0:2].isdigit():
#                             temp = temp[3:]
#                         elif temp[0:2] == "an":
#                             temp = temp[3:]
#                 if len(temp) > 2:
#                     if temp[1] == "_":
#                         if temp[0].isdigit():
#                             temp = temp[2:]
#                         elif temp[0] == "a":
#                             temp = temp[2:]
#                 tmpList.append(temp)
#             temp = ""
#         else:
#             temp += w
#     return tmpList

# inputfile = open("parseddata_2.txt", encoding="utf-8", mode='r')
# fw = open("outputNode.txt", encoding="utf-8", mode='w')
# print("parsing.......")
# for line in inputfile:
#     find = False
#     if line[0:2] == "RN":
#         fw.write(line)
#     else:
#         tmpList = split(line)
#         tempstr_NADP = ""
#         tempstr = ""
#         is_found = False
#         for x in tmpList:
#             if x[0:6] == "NAD(P)":
#                 is_found = True
#                 tempstr += "NAD" + x[6:] + " "
#                 tempstr_NADP += "NADP" + x[6:] + " "
#             else:
#                 tempstr += x + " "
#                 tempstr_NADP += x + " "
#         tempstr += "\n"
#         tempstr_NADP += "\n"
#         if is_found:
#             fw.write(tempstr_NADP)
#         fw.write(tempstr)


########################################################################
# inputfile = open("parseddata_origin.txt", encoding="utf-8", mode='r')
# fw = open("outputNode.txt", encoding="utf-8", mode='w')
# tmpline  = ""
# tmpList  = []
# HashTable= {}
# HashKey  = ""
# tempstr  = ""
# stack1   = []
# stack2   = []
# temp     = ""
# print("parsing.......")
# for line in inputfile:
#     reactant = []
#     product  = []
#     is_product = False

#     if line[0:2] == "RN":
#         fw.write(line)
#     else:
#         tmpList = split(line)
#         for x in tmpList[1:]:
#             if x == "=":
#                 is_product = True
#             elif x == "{ir}":
#                 pass
#             elif x == "{r}" :
#                 tempstr = ""
#                 if len(product) != 0:
#                     for i in product:
#                         tempstr += i + " "
#                     tempstr += " = "
#                     for i in reactant:
#                         tempstr += i + " "
#                     tempstr += "\n"
#                     fw.write(tempstr)

#             elif x != "more" and x != "More" and x != "?"  and x[0] != "{":
#                 if is_product:
#                     product.append(x)
#                 else:
#                     reactant.append(x)
#     tempstr = ""
#     if len(reactant) != 0:
#         for i in reactant:
#             tempstr += i + " "
#         tempstr += " = "
#         for i in product:
#             tempstr += i + " "
#         tempstr += "\n"
#         fw.write(tempstr)







    # tempstr = ""
    # HashKey = ""
    # tmpList = split(line)
    # mark = 0
    # #print(tmpList)
    # fw.write(tmpList[0] + " ")
    # for x in tmpList[1:]:
    #     if x == "+":
    #         fw.write(tempstr[:-1] + " ")
    #         tempstr = ""
    #     elif x == "=":
    #         fw.write(tempstr[:-1] + " = ")
    #         tempstr = ""
    #     elif x == "{r}":
    #         fw.write(tempstr[:-1] + " {r}")
    #         tempstr = ""
    #     elif x == "{ir}":
    #         fw.write(tempstr[:-1] + " {ir}")
    #         tempstr = ""
    #     else:
    #         tempstr += x + "_"
    # fw.write(tempstr[:-1] + "\n")
    # tempstr = ""

inputfile.close()
fw.close()