import json
import os, sys
import gzip
from math import log

from goatools import obo_parser, GOEnrichmentStudy
import Bio.UniProt.GOA as GOA
import requests
import matplotlib.pyplot as plt

import utils
import counter

DATA_PATH = 'data'
FILE_NAME = 'go-basic.obo'

current_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(current_path, DATA_PATH, FILE_NAME)
if not utils.check_file_exists(data_path):
    print("Error: file does not exist.")
    sys.exit(1)

go = obo_parser.GODag(data_path)

print("2.1a: ", go.query_term('GO:0043534').name)

print("2.1b: ", " | ".join([i.name for i in go.query_term('GO:0043534').parents]))

print("2.1c: ", " | ".join([i.name for i in go.query_term('GO:0043534').children]))

all_parents = []
parent_queue = go.query_term('GO:0043534').parents
while len(parent_queue) > 0:
    parent = parent_queue[0]
    parent_queue = parent_queue[1:]
    parent_queue += parent.parents
    all_parents.append(parent.name)
print("2.1d (all parents): ", " | ".join(all_parents))

all_children = []
children_queue = go.query_term('GO:0043534').children
while len(children_queue) > 0:
    child = children_queue[0]
    children_queue = children_queue[1:]
    children_queue += child.children
    all_children.append(child.name)
print("2.1d (all children): ", " | ".join(all_children))

growth = sum(["growth" in i[1].name for i in list(go.items())]) # Number of True's
print("2.1e: ", growth)

a = set([i.id for i in go.query_term('GO:0048527').parents])
b = set([i.id for i in go.query_term('GO:0097178').parents])
last_expanded = 'b'
while len(a.intersection(b)) < 1:
    if last_expanded == 'b':
        a_parents = []
        for p in list(a):
            a_parents += [i.id for i in go.query_term(p).parents]
        last_expanded = 'a'
        a = set(a_parents)
    else:
        b_parents = []
        for p in list(b):
            b_parents += [i.id for i in go.query_term(p).parents]
        last_expanded = 'b'
        b = set(b_parents)
print('2.1f: ', go.query_term(list(a.intersection(b))[0]).name)

#file_path = os.path.join(current_path, DATA_PATH, '2_2_a_lineage.png')
#go.draw_lineage('GO:0097190', lineage_img=file_path)
# 2.2, 2.3 do not work due to pygraphviz not working.

# XML nolonger available, using JSON.
json_result = requests.get("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0048527/complete").json()["results"][0]
print('2.3a: ', json_result["name"], " ", json_result["definition"]["text"])
json_result = requests.get("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0097178/complete").json()["results"][0]
print('2.3b: ', json_result["synonyms"])

with gzip.open('./goa_arabidopsis.gaf.gz', 'rt') as arab_gaf_fp:
    arab_funcs = {}
    for entry in GOA.gafiterator(arab_gaf_fp):
        uniprot_id = entry.pop('DB_Object_ID')
        arab_funcs[uniprot_id] = entry

with_not = 0
have_annotation = 0
with_growth = []
growth_dict = {}
code_count = {}
a = 0
for i in arab_funcs:
    if 'NOT' in arab_funcs[i]['Qualifier']:
        with_not += 1
    if arab_funcs[i]['GO_ID'] == 'GO:0048527':
        have_annotation += 1
    if 'growth' in arab_funcs[i]['DB_Object_Name']:
        with_growth.append(arab_funcs[i]['GO_ID'])
        growth_dict[i] = arab_funcs[i]
    code_count[arab_funcs[i]['Evidence']] = code_count.get(arab_funcs[i]['Evidence'], 0) + 1
print('3.1a: {} with NOT ({}%).'.format(with_not, float(with_not)/len(arab_funcs)*100))
print('3.1b: {} have annotation GO:0048527.'.format(have_annotation))
print('3.1c: the following have "growth" in name:', with_growth)
print('3.1d: codes: ')
for i in code_count.items():
    print(i[0], ': ', i[1], '; ')
print('3.1f: plotting with matplotlib, close the window to continue.')
labels = [i[0] for i in code_count.items()]
sizes = [i[1] for i in code_count.items()]
patches, _ = plt.pie(sizes, startangle=90)
plt.legend(patches, labels, loc="best")
plt.axis('equal')
plt.tight_layout()
plt.show()

pop = arab_funcs.keys()
assoc = {}

for x in arab_funcs:
    if x not in assoc:
        assoc[x] = set()
    assoc[x].add(str(arab_funcs[x]['GO_ID']))

study = growth_dict.keys()

enrichment = GOEnrichmentStudy(pop, assoc, go)
results = enrichment.run_study(study)
enrichment.print_summary(results)
print("4.1a: GO:0030154 has the least depth and therefore most enriched.")

enrichment = GOEnrichmentStudy(pop, assoc, go, alpha=0.01, methods=['bonferroni'])
results = enrichment.run_study(study)
enriched = len([i for i in results if i.p_bonferroni <= 0.01])
print("4.1b: {} terms enriched under Bonferroni p=0.01.".format(enriched))

enrichment = GOEnrichmentStudy(pop, assoc, go, alpha=0.01, methods=['fdr'])
results = enrichment.run_study(study)
enriched = len([i for i in results if i.p_fdr <= 0.01])
print("4.1c: {} terms enriched under false discovery rate p=0.01.".format(enriched))

term_counts = counter.TermCounts(go, arab_funcs)

print("5.1a: Figure 1 has changed, but based on https://i.imgur.com/usDZhxl.jpg, with GO:0044707 merged into GO:0032501, at 4 branches the semantic similarity is 1/4.")
print("5.1b: Information content of GO:0048364 = {}".format(-1*log(term_counts.get_term_freq('GO:0048364'), 2)))

a = set([i.id for i in go.query_term('GO:0048364').parents])
b = set([i.id for i in go.query_term('GO:0044707').parents])
last_expanded = 'b'
while len(a.intersection(b)) < 1:
    if last_expanded == 'b':
        for p in list(a):
            a_parents += [i.id for i in go.query_term(p).parents]
        last_expanded = 'a'
        a = set(a_parents)
    else:
        for p in list(b):
            b_parents += [i.id for i in go.query_term(p).parents]
        last_expanded = 'b'
        b = set(b_parents)
common_ancestor = go.query_term(list(a.intersection(b))[0]).id
print("5.1c: Deepest common ancestor of GO:0044707 and GO:0048364: ", common_ancestor)
print("Therefore the Resnik similarity is {}".format(-1*log(term_counts.get_term_freq(common_ancestor), 2)))
