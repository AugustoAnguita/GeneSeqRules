#GNU GENERAL PUBLIC LICENSE - Version 3
#
# pmml.py Copyright (C) 2019
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


import sys

if len(sys.argv) < 3:
	print "usage: python", sys.argv[0], "<items_names> <rules_file>"
	exit(-1)

# read correspondences file
f = open(sys.argv[1], "r")
lines = f.readlines()
correspondences = {}
for l in lines:
	elements = l.replace("\n", "").replace('\r', "").replace('"', "").split("\t")
	number = elements[2]
	name = elements[1]
	correspondences[number] = name

# read rules file
f = open(sys.argv[2], "r")
lines = f.readlines()

itemsets = []
rules = []

for l in lines:
	elements = l.replace("\r\n", "").split("#")
	# 1016,1029,1104,1247 ==> 1368,2135 #SUP: 11 #CONF: 1.0 #LIFT: 2.0
	# the first element is the rule, the rest of elements are the measures.
	rule = elements[0].replace(" ", "").split("==>")
	antecedent = rule[0].split(",")
	consequent = rule[1].split(",")
	measures = {}

	for i in range(1, len(elements)):
		m = elements[i].replace(" ", "").split(":")
		name = m[0]
		value = m[1]
		measures[name] = value

	itemsets.append(antecedent)
	itemsets.append(consequent)

	pmmlRule = {
		"antecedent": len(itemsets) - 2,
		"consequent": len(itemsets) - 1,
		"measures": measures
	}
	rules.append(pmmlRule)

# print itemsets
# print rules
# print correspondences

print '''
<PMML xmlns="http://www.dmg.org/PMML-4_3" version="4.3">
<Header copyright="www.dmg.org" description="example model for association rules"/>
<AssociationModel functionName="associationRules" numberOfTransactions="0" numberOfItems="{}" minimumSupport="0" minimumConfidence="0" numberOfItemsets="{}" numberOfRules="{}">
<MiningSchema>
</MiningSchema>
'''.format(len(correspondences), len(itemsets), len(rules))

print "<!-- Items -->"
for i in correspondences:
	print '<Item id="{}" value="{}=1"/>'.format("1"+i, correspondences[i])
	print '<Item id="{}" value="{}=2"/>'.format("2"+i, correspondences[i])


print "\n<!-- Itemsets -->"
for i in range(0, len(itemsets)):
	print '<Itemset id="{}" numberOfItems="{}">'.format(i, len(itemsets[i]))
	for ii in itemsets[i]:
		print '<ItemRef itemRef="{}"/>'.format(ii)
	print '</Itemset>'

print "\n<!-- Rules -->"
for r in rules:
	measures = ''.join(['{}="{}" '.format(name, r['measures'][name]) for name in r['measures']])
	print '<AssociationRule {}antecedent="{}" consequent="{}"/>'.format(measures, r['antecedent'], r['consequent'])


print '''
</AssociationModel>
</PMML>
'''

