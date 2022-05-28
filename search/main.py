#update

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
    # input_species = [searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    # output_species = searcher.mapToNode["NADH"]
    # output_species = searcher.mapToNode["2,2'-azino-bis(3-ethylbenzthiazoline-6-sulfonic_acid)"]
    # input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    # output_species = searcher.mapToNode["NADH"]
    # output_species = searcher.mapToNode["resorufin"]
    input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    output_species = searcher.mapToNode["NADH"]
    # output_species = [searcher.mapToNode["resorufin"], searcher.mapToNode["NADH"]] 

    # searcher.test()
    searcher.main(input_species, output_species)
    # print (searcher.c)

    