from node import node

from edge import edge

from func import search

if __name__ == '__main__':
    searcher = search()
    f = open("parseddata.txt", encoding="utf-8", mode='r')
    searcher.parsing(f)
    f.close()
    fl = open("label.txt", encoding="utf-8", mode='r')
    searcher.parsing_label(fl)
    fl.close()
    print("num of reaction : ", searcher.reaction)
    print("num of species : ", len(searcher.nodeList))
    
    temp = 0

    print("typeA : ", searcher.typeA)
    print("typeB : ", searcher.typeB)
    print("typeC : ", searcher.typeC)
    print("typeD : ", searcher.typeD)
    print("typeE : ", searcher.typeE)

    # count = 0
    print ("\n", "*************************************************")

    #specified input species :
    input_species = [searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    # input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    # input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["lactate"]]

    searcher.main(input_species)
    # print (searcher.c)

    