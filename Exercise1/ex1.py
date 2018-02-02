import os, sys
from goatools import obo_parser

import utils

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

file_path = os.path.join(current_path, DATA_PATH, '2_2_a_lineage.png')
go.draw_lineage('GO:0097190', lineage_img=file_path)
