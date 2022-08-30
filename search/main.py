<<<<<<< HEAD
#update
import sys

import time

from func import search

if __name__ == '__main__':
    start = time.process_time()
    searcher = search()
    f = open("parseddata.txt", encoding="utf-8", mode='r')
    # f = open("testsys.txt", encoding="utf-8", mode='r')
    searcher.parsing(f)
    f.close()
    fl = open("label.txt", encoding="utf-8", mode='r')
    searcher.parsing_label(fl)
    fl.close()
    print("num of reaction : ", searcher.reaction)
    print("num of species : ", len(searcher.nodeList))
    
    temp = 0

    print("typeA : ", searcher.typeB)
    print("typeB : ", searcher.typeA)
    # print("typeC : ", searcher.typeE)
    # print("typeD : ", searcher.typeD)
    # print("typeE : ", searcher.typeC)
    # print("typeF : ", searcher.typeF)
    # print("typeG : ", searcher.reaction)

    # count = 0
    print ("\n", "*************************************************")
    input_species = []
    if sys.argv[1] == 'le':
        input_species = [searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    elif sys.argv[1] == 'ag':
        input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    elif sys.argv[1] == 'gn':
        input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    elif sys.argv[1] == 'all':
        input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"], searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    # elif sys.argv[1] == 'gn':
    #     input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    # elif sys.argv[1] == 'gn':
    #     input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    else:
        print ("no input")
    output_species = ["NADH", "resorufin", "2,2'-azino-bis(3-ethylbenzothiazoline-6-sulphonic_acid)", " 4,5-Diaminofluorescein"]
    # output_species = searcher.mapToNode["resorufin"]

    if input_species != []:
        searcher.main(input_species, output_species)
        # searcher.XORmain(input_species, output_species)

    end = time.process_time()
    print("總共執行時間：%f 秒" % (end - start))
    print("檢查violation：%f 秒" % (searcher.time))
    # print ("return true")


    #specified input species :
    # output_species = searcher.mapToNode["NADH"]
    # output_species = searcher.mapToNode["2,2'-azino-bis(3-ethylbenzthiazoline-6-sulfonic_acid)"]
    # input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    # output_species = searcher.mapToNode["NADH"]
    # output_species = searcher.mapToNode["resorufin"]
    # input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    # output_species = [searcher.mapToNode["resorufin"], searcher.mapToNode["NADH"]] 

    # searcher.test()
    # print (searcher.c)
    # a = [searcher.mapToNode["NAD+"], searcher.mapToNode["phenol"], searcher.mapToNode["FADH2"], searcher.mapToNode["NADH"], searcher.mapToNode["trans,octacis-decaprenylphospho-beta-D-ribofuranose"], searcher.mapToNode["2-bromobenzoate"], searcher.mapToNode["H2O2"], searcher.mapToNode["FAD"], searcher.mapToNode["resazurin"], searcher.mapToNode["O2"], searcher.mapToNode["bromide"], searcher.mapToNode["reduced_ferredoxin_[iron-sulfur]_cluster"], searcher.mapToNode["1H-pyrrole-2-carbonyl-[PltL_peptidyl-carrier_protein]"], searcher.mapToNode["ethanol"], searcher.mapToNode["phenyl_phosphate"], searcher.mapToNode["4-phenoxyphenol"], searcher.mapToNode["lactate"], searcher.mapToNode["cytidine"], searcher.mapToNode["NADPH"]]
    # b = [searcher.mapToNode["alcohol_dehydrogenase"], searcher.mapToNode["lactate_dehydrogenase"], searcher.mapToNode["2-halobenzoate_1,2-dioxygenase"], searcher.mapToNode["alcohol_oxidase"], searcher.mapToNode["lignin_peroxidase"], searcher.mapToNode["unspecific_peroxygenase"], searcher.mapToNode["1H-pyrrole-2-carbonyl-[peptidyl-carrier_protein]_chlorinase"], searcher.mapToNode["decaprenylphospho-beta-D-ribofuranose_2-dehydrogenase"], searcher.mapToNode["nucleoside_phosphotransferase"]]
    # allnodes = set()
    # for x in a:
    #     allnodes.add(x)
    # for x in b:
    #     allnodes.add(x)
    # allrecs = searcher.CollectAll(allnodes, [])
    # c = 0
    # for i in allrecs:
    #     print (i.show(c))
    #     c =c +1



=======
#update
import sys

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

    print("typeA : ", searcher.typeB)
    print("typeB : ", searcher.typeA)
    print("typeC : ", searcher.typeE)
    print("typeD : ", searcher.typeD)
    print("typeE : ", searcher.typeC)
    print("typeF : ", searcher.typeF)
    print("typeG : ", searcher.reaction)

    # count = 0
    print ("\n", "*************************************************")
    input_species = []
    if sys.argv[1] == 'le':
        input_species = [searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    elif sys.argv[1] == 'ag':
        input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    elif sys.argv[1] == 'gn':
        input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    elif sys.argv[1] == 'all':
        input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"], searcher.mapToNode["lactate"], searcher.mapToNode["ethanol"]]
    # elif sys.argv[1] == 'gn':
    #     input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    # elif sys.argv[1] == 'gn':
    #     input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    else:
        print ("no input")
    output_species = ["NADH", "resorufin", "2,2'-azino-bis(3-ethylbenzothiazoline-6-sulphonic_acid)", " 4,5-Diaminofluorescein"]
    # output_species = searcher.mapToNode["resorufin"]

    if input_species != []:
        searcher.main(input_species, output_species)

    #specified input species :
    # output_species = searcher.mapToNode["NADH"]
    # output_species = searcher.mapToNode["2,2'-azino-bis(3-ethylbenzthiazoline-6-sulfonic_acid)"]
    # input_species = [searcher.mapToNode["acetone"], searcher.mapToNode["D-glucose"]]
    # output_species = searcher.mapToNode["NADH"]
    # output_species = searcher.mapToNode["resorufin"]
    # input_species = [searcher.mapToNode["D-glucose"], searcher.mapToNode["nitric_oxide"]]
    # output_species = [searcher.mapToNode["resorufin"], searcher.mapToNode["NADH"]] 

    # searcher.test()
    # print (searcher.c)
    # a = [searcher.mapToNode["NAD+"], searcher.mapToNode["phenol"], searcher.mapToNode["FADH2"], searcher.mapToNode["NADH"], searcher.mapToNode["trans,octacis-decaprenylphospho-beta-D-ribofuranose"], searcher.mapToNode["2-bromobenzoate"], searcher.mapToNode["H2O2"], searcher.mapToNode["FAD"], searcher.mapToNode["resazurin"], searcher.mapToNode["O2"], searcher.mapToNode["bromide"], searcher.mapToNode["reduced_ferredoxin_[iron-sulfur]_cluster"], searcher.mapToNode["1H-pyrrole-2-carbonyl-[PltL_peptidyl-carrier_protein]"], searcher.mapToNode["ethanol"], searcher.mapToNode["phenyl_phosphate"], searcher.mapToNode["4-phenoxyphenol"], searcher.mapToNode["lactate"], searcher.mapToNode["cytidine"], searcher.mapToNode["NADPH"]]
    # b = [searcher.mapToNode["alcohol_dehydrogenase"], searcher.mapToNode["lactate_dehydrogenase"], searcher.mapToNode["2-halobenzoate_1,2-dioxygenase"], searcher.mapToNode["alcohol_oxidase"], searcher.mapToNode["lignin_peroxidase"], searcher.mapToNode["unspecific_peroxygenase"], searcher.mapToNode["1H-pyrrole-2-carbonyl-[peptidyl-carrier_protein]_chlorinase"], searcher.mapToNode["decaprenylphospho-beta-D-ribofuranose_2-dehydrogenase"], searcher.mapToNode["nucleoside_phosphotransferase"]]
    # allnodes = set()
    # for x in a:
    #     allnodes.add(x)
    # for x in b:
    #     allnodes.add(x)
    # allrecs = searcher.CollectAll(allnodes, [])
    # c = 0
    # for i in allrecs:
    #     print (i.show(c))
    #     c =c +1



>>>>>>> d9fa1734b190534bfa52384f299936f809ef9a8f
    