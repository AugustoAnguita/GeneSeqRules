#GNU GENERAL PUBLIC LICENSE - Version 3
#
# main.py Copyright (C) 2019
#
# This program is free software: you can redistribute it and/or modify
#
# it under the terms of the GNU General Public License as published by
#
# the Free Software Foundation, either version 3 of the License, or
#
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#
# but WITHOUT ANY WARRANTY; without even the implied warranty of
#
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#
# GNU General Public License for more details.
#
# 
#
# You should have received a copy of the GNU General Public License
#
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information: augustoanguitaruiz@gmail.com
#
#Anguita-Ruiz A, Segura-Delgado A, Alcala R & Alcala-Fdez J


import sys, os
import xml.etree.ElementTree
from xml.dom import minidom

def readPMMLFile(path):
    root = xml.etree.ElementTree.parse(path).getroot()
    ns = {
        "d": root.tag[root.tag.find("{")+1:root.tag.find("}")]
    }

    items = {}
    itemsets = {}
    rules = {}

    for item in root.findall("d:AssociationModel/d:Item", ns):
        value = item.get("value")
        probe = value.split("/")[0]
        gen = value.split("/")[1]
        expression = gen.split("=")[1]
        gen = gen.split("=")[0]

        items[item.get("id")] = {
            "probe": probe,
            "gen": gen,
            "expression": expression
        }
    #print(items)

    for item in root.findall("d:AssociationModel/d:Itemset", ns):
        litems = []
        for i in item:
            litems.append(i.get("itemRef"))

        itemsets[item.get("id")] = litems
    #print(itemsets["98"])

    for item in root.findall("d:AssociationModel/d:AssociationRule", ns):
        id = "{}:{}".format(item.get("antecedent"), item.get("consequent"))
        rules[id] = [
            item.get("antecedent"),
            item.get("consequent")
        ]
    #print(rules)

    return (items, itemsets, rules)

def readGOCSV(path):
    # Format MUST be:
    # "";"PROBEID";"SYMBOL";"GO";"ONTOLOGY";"GOALL";"ONTOLOGYALL";"PATH"
    # including the header in the first line

    probes_go = {}
    gos = {}

    f = open(path, "r")
    for l in f.readlines()[1:]:
        row = l.replace("\"", "").replace("\n", "").split(";")
        probeid = row[1]
        GO = row[3]
        GOALL = row[5]
        ONCOLOGY = row[4]
        ONCOLOGYALL = row[6]

        if probeid not in probes_go:
            probes_go[probeid] = {
                "BP": {},
                "MF": {},
                "CC": {},
                "path":{
                    row[7]: True
                }
            }
            if ONCOLOGY in probes_go[probeid]:
                probes_go[probeid][ONCOLOGY][GO] = True
                probes_go[probeid][ONCOLOGY][GOALL] = True
            if ONCOLOGYALL in probes_go[probeid]:
                probes_go[probeid][ONCOLOGYALL][GO] = True
                probes_go[probeid][ONCOLOGYALL][GOALL] = True
        else:
            if ONCOLOGY in probes_go[probeid]:
                probes_go[probeid][ONCOLOGY][GO] = True
                probes_go[probeid][ONCOLOGY][GOALL] = True
            if ONCOLOGYALL in probes_go[probeid]:
                probes_go[probeid][ONCOLOGYALL][GO] = True
                probes_go[probeid][ONCOLOGYALL][GOALL] = True
            probes_go[probeid]["path"][row[7]] = True

        gos[GO] = True
        gos[GOALL] = True

    f.close()

    return (probes_go, gos)

def readTRRUSTCSV(path):
    # Format MUST be:
    # FACTOR DE TRANSCRIPCION;GENES DIANA;MODO DE REGULACION;ARTICULOS QUE LO SOPORTAN
    # including the header in the first line

    # { GEN_ID: { DIANA_GEN:{ GEN_ID:REGULATION } } }
    transcription_factor = {}

    f = open(path, "r")
    for l in f.readlines()[1:]:
        row = l.replace("\"", "").replace("\n", "").split(";")
        genid = row[0]
        dianagen = row[1]
        regulation = row[2]

        if genid not in transcription_factor:
            transcription_factor[genid] = {
                "diana_gen": {
                    dianagen: regulation
                }
            }
        else:
            transcription_factor[genid]["diana_gen"][dianagen] = regulation

    f.close()

    return transcription_factor

def calculateMeasureWithMatches(ant_items, con_items, rule_matches):
    ant = 0.0
    con = 0.0
    total = 0.0
    ant_matches = 0.0
    con_matches = 0.0
    if len(rule_matches) < 1:
        return (ant, con, total)

    for go in rule_matches:
        ant += len(rule_matches[go]["ant"])
        con += len(rule_matches[go]["con"])
        if len(rule_matches[go]["ant"]) > 0:
            ant_matches += 1
        if len(rule_matches[go]["con"]) > 0:
            con_matches += 1

    total = ant + con
    total_items = len(ant_items) + len(con_items)

    # make it [0, 1]
    if ant_matches > 0:
        ant = float(ant) / float(len(ant_items)*ant_matches)
    if con_matches > 0:
        con = float(con) / float(len(con_items)*con_matches)
    if len(rule_matches) > 0:
        total = float(total) / float(total_items*len(rule_matches))

    return (ant, con, total)

def calculateProbeMeasures(items, itemsets, rules, probes_go, gos):
    rules_measures = {} # RULE_ID -> {MEASURE: {RAW: Value, LEVEL: level}}
    rules_matches = {} # RULE_ID -> {ONCOLOGY:{ GO: {ANT: [ITEM_ID], CON: [ITEM_ID]} }}
    rules_matches_sp = {} # RULE_ID -> {PATH: {ANT: [ITEM_ID], CON: [ITEM_ID]} }

    # var to save max and min value for measure
    maxmin = {
        "BP": [0, 999999999],
        "MF": [0, 999999999],
        "CC": [0, 999999999],
        "path": [0, 999999999]
    }

    # Check matches for each rule and item
    for r in rules:
        ant = rules[r][0]
        con = rules[r][1]
        ant_itemset = itemsets[ant]
        con_itemset = itemsets[con]

        rules_matches[r] = {
            "BP":{},
            "MF":{},
            "CC":{},
            "path": {}
        }
        rules_matches_sp[r] = {}

        rules_measures[r] = {
            "BP": {},
            "MF": {},
            "CC": {},
            "path": {}
        }
        rules_measures[r]["BP"] = {
            "raw": 0,
            "level": 0
        }
        rules_measures[r]["MF"] = {
            "raw": 0,
            "level": 0
        }
        rules_measures[r]["CC"] = {
            "raw": 0,
            "level": 0
        }
        rules_measures[r]["path"] = {
            "raw": 0,
            "level": 0
        }

        for o in ["BP", "MF", "CC", "path"]:
            gos_matches = {
                r: {}
            }
            for g in gos:
                # for each oncology, lets check each probe
                ant_matches = []
                for i in ant_itemset:
                    probe = items[i]["probe"]
                    if probe in probes_go and o in probes_go[probe] and g in probes_go[probe][o]:
                        ant_matches.append(i)
                con_matches = []
                for i in con_itemset:
                    probe = items[i]["probe"]
                    if probe in probes_go and o in probes_go[probe] and g in probes_go[probe][o]:
                        con_matches.append(i)

                # check level 1
                if len(ant_matches) == len(ant_itemset) and len(con_matches) == len(con_itemset):
                    # all have been matched, it is level 1
                    rules_measures[r][o]["raw"] += 4
                    rules_measures[r][o]["level"] = 1
                else:
                    # save matches for level combinations
                    gos_matches[r][g] = {
                        "ant": ant_matches,
                        "con": con_matches
                    }

                # save matches
                rules_matches[r][o][g] = {
                    "ant": ant_matches,
                    "con": con_matches
                }

            # all gos matches saved, lets check for levels 2,3,4
            # check level 2
            for g in gos:
                if g not in gos_matches or len(gos_matches[g]["con"]) <= 0:
                    # no match in consequent -> no level 2
                    continue

                combined_gos = [g] # list of combined gos to achieve level 2
                not_matched_ant = [ x for x in ant_itemset if x not in gos_matches[g]["ant"] ] # items ant no matcheados
                not_matched_con = [x for x in con_itemset if x not in gos_matches[g]["con"]]  # items con no matcheados

                for gg in gos_matches:
                    # compruebo, de los no matcheados por el GO anterior, los que si matchea este
                    not_matched_ant2 = [x for x in not_matched_ant if x not in gos_matches[g]["ant"]]  # items ant no matcheados
                    not_matched_con2 = [x for x in not_matched_con if x not in gos_matches[g]["con"]]  # items con no matcheados

                    if len(not_matched_ant2) <= 0 and len(not_matched_con2) < len(con_itemset):
                        # todos matcheados con al menos 1 consecuente! es nivel 2. Elimino los GO del diccionario para no volverlos a usar
                        rules_measures[r][o]["raw"] += 3
                        if rules_measures[r][o]["level"] == 0:
                            rules_measures[r][o]["level"] = 2

                        for ggg in combined_gos:
                            # delete sombined gos to do the following combinations
                            del gos_matches[ggg]
                        break # stop this loop and continue searching matches with the remaining gos
                    elif len(not_matched_ant2) < len(not_matched_ant) and len(not_matched_con2) < len(not_matched_con):
                        # he conseguido mas matches, pero aun no he combinado suficientes..
                        # anado el go a la lista de combinaciones, y actualizo la lista de matches
                        combined_gos.append(gg)
                        not_matched_ant = not_matched_ant2
                        not_matched_con = not_matched_con2

            # once the level 2 checks end, we have in gos_matches the remaining
            # gos matches to continue checking level 3 and 4

            # check level 3
            for g in gos:
                if g not in gos_matches:
                    continue
                if len(gos_matches[g]["ant"]) > 1 and len(gos_matches[g]["con"]) > 1:
                    rules_measures[r][o]["raw"] += 2
                    if rules_measures[r][o]["level"] == 0:
                        rules_measures[r][o]["level"] = 3
                    del gos_matches[g] # delete go, it is level 3

            # check level 4
            for g in gos:
                if g not in gos_matches:
                    continue
                if len(gos_matches[g]["ant"]) > 2 or len(gos_matches[g]["con"]) > 2:
                    rules_measures[r][o]["raw"] += 1
                    if rules_measures[r][o]["level"] == 0:
                        rules_measures[r][o]["level"] = 4
                    del gos_matches[g]

            # update max min values for each measure
            if maxmin[o][0] < rules_measures[r][o]["raw"]:
                maxmin[o][0] = rules_measures[r][o]["raw"]
            if maxmin[o][1] > rules_measures[r][o]["raw"]:
                maxmin[o][1] = rules_measures[r][o]["raw"]

    for r in rules:
        for o in ["BP", "MF", "CC", "path"]: 
            if rules_measures[r][o]["level"] == 0:
                rules_measures[r][o]["level"] = 5
        
        rules_measures[r] = {
            "BP": rules_measures[r]["BP"]["level"] + (1.0 - float(rules_measures[r]["BP"]["raw"] / float(maxmin["BP"][0]+1))),

            "MF": rules_measures[r]["MF"]["level"] + (1.0 - float(rules_measures[r]["MF"]["raw"] / float(maxmin["MF"][0]+1))),

            "CC": rules_measures[r]["CC"]["level"] + (1.0 - float(rules_measures[r]["CC"]["raw"] / float(maxmin["CC"][0]+1))),

            "SP": rules_measures[r]["path"]["level"] + (1.0 - float(rules_measures[r]["path"]["raw"] / float(maxmin["path"][0]+1)))
        }

    return (rules_matches, rules_matches_sp, rules_measures)


def calculateTF(items, itemsets, rules, transcription_factor):
    # Transcription Factor
    rules_matches = {}
    rules_measures = {}

    for r in rules:
        rules_measures[r] = {"TF": 0} # initialize measure accumulate
        rules_matches[r] = []
        ant = rules[r][0]
        con = rules[r][1]
        for a in itemsets[ant]:
            a_item = items[a]
            if a_item["gen"] in transcription_factor:
                rules_measures[r]["TF"] += 1
                dianas = transcription_factor[a_item["gen"]]["diana_gen"]
                for c in itemsets[con]:
                    c_item = items[c]
                    if c_item["gen"] in dianas:
                        rules_measures[r]["TF"] += 1
                        if dianas[c_item["gen"]] == "Activation" and a_item["expression"]==c_item["expression"]:
                            rules_measures[r]["TF"] += 1
                            rules_matches[r].append([a, c, dianas[c_item["gen"]]])
                        elif dianas[c_item["gen"]] == "Repression" and a_item["expression"]!=c_item["expression"]:
                            rules_measures[r]["TF"] += 1
                            rules_matches[r].append([a, c, dianas[c_item["gen"]]])


    return (rules_matches, rules_measures)

def saveRulesToPMMLFile(originalpath, rules, rules_measures):
    root = xml.etree.ElementTree.parse(originalpath).getroot()
    ns = {
        "d": root.tag[root.tag.find("{") + 1:root.tag.find("}")]
    }

    for item in root.findall("d:AssociationModel/d:AssociationRule", ns):
        rid = "{}:{}".format(item.get("antecedent"), item.get("consequent"))
        if rid in rules_measures:
            for measure in rules_measures[rid]:
                item.set(measure, str(rules_measures[rid][measure]))

    output = originalpath.replace(".pmml", "_biologicmeasures.pmml")
    f = open(output, "w")
    f.write(xml.etree.ElementTree.tostring(root).replace("ns0:", ""))
    f.close()

def saveOncologyMatchesToXML(originalpath, items, itemsets, rules, rules_matches):
    output = originalpath.replace(".pmml", "_GO_matches.xml")
    f = open(output, "w")
    f.write('<data>\n')

    for r in rules_matches:

        f.write('<rule id="{}" antecedent="{}" consequent="{}">\n'.format(r, rules[r][0], rules[r][1]))

        antecedent = itemsets[rules[r][0]]
        consequent = itemsets[rules[r][1]]
        for onc in rules_matches[r]:
            f.write('\t<oncology id="{}">\n'.format(onc))
            for go in rules_matches[r][onc]:
                if len(rules_matches[r][onc][go]["ant"]) == 0 and len(rules_matches[r][onc][go]["con"]) == 0:
                    continue
                f.write('\t\t<go id="{}">\n'.format(go))
                for i in rules_matches[r][onc][go]["ant"]:
                    f.write('\t\t\t<antecedent id="{}" value="{}/{}={}" />\n'.format(str(i), items[i]["probe"], items[i]["gen"], items[i]["expression"]))
                for i in rules_matches[r][onc][go]["con"]:
                    f.write('\t\t\t<consequent id="{}" value="{}/{}={}" />\n'.format(str(i), items[i]["probe"], items[i]["gen"], items[i]["expression"]))

                f.write('\t\t</go>\n')
            f.write('\t</oncology>\n')

        f.write('</rule>\n')
    f.write('</data>')

    f.close()

def saveSPMatchesToXML(originalpath, items, itemsets, rules, rules_matches):
    output = originalpath.replace(".pmml", "_SP_matches.xml")
    f = open(output, "w")
    f.write('<data>\n')

    for r in rules_matches:
        f.write('<rule id="{}" antecedent="{}" consequent="{}">\n'.format(r, rules[r][0], rules[r][1]))
        antecedent = itemsets[rules[r][0]]
        consequent = itemsets[rules[r][1]]
        for p in rules_matches[r]:
            f.write('\t<path id="{}">\n'.format(p))
            for i in rules_matches[r][p]["ant"]:
                f.write('\t\t<antecedent id="{}" value="{}/{}={}" />\n'.format(str(i), items[i]["probe"], items[i]["gen"], items[i]["expression"]))
            for i in rules_matches[r][p]["con"]:
                f.write('\t\t<consequent id="{}" value="{}/{}={}" />\n'.format(str(i), items[i]["probe"], items[i]["gen"], items[i]["expression"]))

            f.write('\t</path>\n')
        f.write('</rule>\n')

    f.write('</data>')
    f.close()

def saveTFMatchesToXML(originalpath, items, itemsets, rules, rules_matches):
    output = originalpath.replace(".pmml", "_TF_matches.xml")
    f = open(output, "w")
    f.write('<data>\n')

    for r in rules_matches:
        f.write('<rule id="{}" antecedent="{}" consequent="{}">\n'.format(r, rules[r][0], rules[r][1]))
        antecedent = itemsets[rules[r][0]]
        consequent = itemsets[rules[r][1]]

        for m in rules_matches[r]:
            f.write('\t<match type="{}">\n'.format(rules_matches[r][m][2]))

            i = rules_matches[r][m][0]
            f.write('\t\t<antecedent id="{}" value="{}/{}={}" />\n'.format(str(i), items[i]["probe"], items[i]["gen"], items[i]["expression"]))

            i = rules_matches[r][m][1]
            f.write('\t\t<consequent id="{}" value="{}/{}={}" />\n'.format(str(i), items[i]["probe"], items[i]["gen"], items[i]["expression"]))

            f.write('\t</match>\n')
        f.write('</rule>\n')
    f.write('</data>')
    f.close()


def mergeDicts(x , y):
    for k in y:
        for kk in y[k]:
            x[k][kk] = y[k][kk]

    return x

def navigateDir(dir):
    pmmlFiles = []
    for f in os.listdir(dir):
        filePath = os.path.join(dir, f)
        if os.path.isdir(filePath):
            #print("Directory: {}".format(filePath))
            pmmlFiles = pmmlFiles + navigateDir(filePath)
        else:
            if f.endswith(".pmml"):
                pmmlFiles.append(filePath)

    return pmmlFiles


def main():

    if len(sys.argv) < 4:
        print("Usage: python {} <Directory_OR_PMMLFile_of_rules> <annotations_csv> <TRRUST_csv>".format(sys.argv[0]))
        return

    filePath = sys.argv[1]
    pmmlFiles = [filePath]
    if os.path.isdir(filePath):
        pmmlFiles = navigateDir(filePath)

    # Read BDs with biologic info
    print("Reading {}..".format(sys.argv[2]))
    probes_go, gos = readGOCSV(sys.argv[2])
    print("Reading {}..".format(sys.argv[3]))
    transcription_factor = readTRRUSTCSV(sys.argv[3])

    for path in pmmlFiles:
        # Read rules
        print("Reading {}..".format(path))
        items, itemsets, rules = readPMMLFile(path)

        # Calculate measures
        print("Calculating BP,MF,CC,SP measures..")
        rules_matches, rules_matches_sp, rules_measures = calculateProbeMeasures(items, itemsets, rules, probes_go, gos)
        print("Calculating TF measure..")
        rules_matches_tf, rules_measures_tf = calculateTF(items, itemsets, rules, transcription_factor)

        # Merge all measures on one dict
        rules_measures = mergeDicts(rules_measures, rules_measures_tf)

        print("Saving results files..")
        # Save results
        saveRulesToPMMLFile(path, rules, rules_measures)
        # Save matches
        saveOncologyMatchesToXML(path, items, itemsets, rules, rules_matches)
        saveTFMatchesToXML(path, items, itemsets, rules, rules_matches_tf)

main()
