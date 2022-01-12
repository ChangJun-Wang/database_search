import sys
 
# the setrecursionlimit function is
# used to modify the default recursion
# limit set by python. Using this,
# we can increase the recursion limit
# to satisfy our needs
 
sys.setrecursionlimit(10**6)

from node import node

from edge import edge

from func import search

if __name__ == '__main__':
    #traversal(store)

    searcher = search()
    f = open("parseddata.txt", encoding="utf-8", mode='r')
    searcher.parsing(f)
    f.close()
    print("num of reaction : ", searcher.reaction)
    print("num of species : ", len(searcher.nodeList))
    
    temp = 0
    fw = open("outputNode.txt", encoding="utf-8", mode='w')


    print("typeA : ", searcher.typeA)
    print("typeB : ", searcher.typeB)
    print("typeC : ", searcher.typeC)
    print("typeD : ", searcher.typeD)
    print("typeE : ", searcher.typeE)

    # for node in searcher.nodeList:
    #     fw.write(nodes.name + " : ")
    #     fw.write("\n")
    fw.close()


    # count = 0
    # notes = ["C","C","C"]
    # print ("\n", "*************************************************")

    # #specified input species :
    # input_species = [searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    # notes         = ["A","A","B","C"]

    # searcher.BFS_all(input_species, "A")
    # print (searcher.c)
    # searcher.ClearVis()
    # searcher.sec_input(input_species, notes)
    # searcher.ClearTmp()




    # allnodes = set()
    # allnodes.add(searcher.mapToNode["lactate"])
    # allnodes.add(searcher.mapToNode["oxidized_acceptor"])
    # allnodes.add(searcher.mapToNode["2-oxo-acid_reductase"])
    # allnodes.add(searcher.mapToNode["pyruvate"])
    # allnodes.add(searcher.mapToNode["reduced_acceptor"])
    # allnodes.add(searcher.mapToNode["[arylsulfatase]-L-serine"])
    # allnodes.add(searcher.mapToNode["S-adenosyl-L-methionine"])
    # allnodes.add(searcher.mapToNode["reduced_acceptor"])
    # allnodes.add(searcher.mapToNode["Ser-type_anaerobic_sulfatase-maturating_enzyme"])
    # allnodes.add(searcher.mapToNode["[arylsulfatase]-3-oxo-L-alanine"])
    # allnodes.add(searcher.mapToNode["L-methionine"])
    # allnodes.add(searcher.mapToNode["5'-deoxyadenosine"])
    # allnodes.add(searcher.mapToNode["acceptor"])
    # allnodes.add(searcher.mapToNode["H2O"])
    # allnodes.add(searcher.mapToNode["4-pyridoxate"])
    # allnodes.add(searcher.mapToNode["4-pyridoxic_acid_dehydrogenase"])
    # allnodes.add(searcher.mapToNode["5-formyl-3-hydroxy-2-methylpyridine-4-carboxylate"])
    # allnodes.add(searcher.mapToNode["NADH"])
    # allnodes.add(searcher.mapToNode["H+"])
    # allnodes.add(searcher.mapToNode["NAD+"])
    # allnodes.add(searcher.mapToNode["(9E)-tetradecenoate"])
    # allnodes.add(searcher.mapToNode["O2"])
    # allnodes.add(searcher.mapToNode["acyl-CoA_11-(E)-desaturase"])
    # allnodes.add(searcher.mapToNode["(9E,11E)-tetradecadienoate"])
    # allnodes.add(searcher.mapToNode["oxidized_acceptor"])
    # allnodes.add(searcher.mapToNode["2-oxo-acid_reductase"])
    # allnodes.add(searcher.mapToNode["ethanol"])
    # allnodes.add(searcher.mapToNode["beta-NAD+"])
    # allnodes.add(searcher.mapToNode["alcohol_dehydrogenase"])
    # allnodes.add(searcher.mapToNode["acetaldehyde"])
    # allnodes.add(searcher.mapToNode["alcohol_oxidase"])
    # allnodes.add(searcher.mapToNode["decanal"])
    # allnodes.add(searcher.mapToNode["succinate-semialdehyde_dehydrogenase_[NAD(P)+]"])
    # allnodes.add(searcher.mapToNode["decanoate"])
    # allnodes.add(searcher.mapToNode["1-octanol"])
    # allnodes.add(searcher.mapToNode["NADP+"])
    # allnodes.add(searcher.mapToNode["octanoic_acid"])
    # allnodes.add(searcher.mapToNode["NADPH"])
    # allnodes.add(searcher.mapToNode["stearoyl-CoA"])
    # allnodes.add(searcher.mapToNode["electron_donor"])
    # allnodes.add(searcher.mapToNode["stearoyl-[acyl-carrier-protein]_9-desaturase"])
    # allnodes.add(searcher.mapToNode["9-octadecenoyl-CoA"])


    # related_rec = searcher.check(allnodes)
    # fw.write(related_rec)
    fw.close()

