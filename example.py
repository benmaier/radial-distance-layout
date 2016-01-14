from radial_distance_layout import radial_distance_layout
import matplotlib.pyplot as pl
import networkx as nx


paths  = [ [ 'a','b','c'] ]
paths += [ [ 'a','b','d'] ]
paths += [ [ 'a','e','f','g'] ]
paths += [ [ 'a','e','f','h'] ]
paths += [ [ 'a','e','i'] ]
paths += [ [ 'a','j','k'] ]
paths += [ [ 'a','j','l'] ]

dists = {'a': 0, 
         'b':1.1, 'e': 1.2, 'j': 1.4,
         'c':2.1, 'd': 2.2, 'f': 2.1, 'i': 2.34, 'k':3.8, 'l':2.5,
         'g': 3.9, 'h': 3.8}
T = nx.DiGraph()

for p in paths:
    T.add_path(p)

keystr = 'dist'

nx.set_node_attributes(T,keystr,dists)

fig,ax = pl.subplots(1,2,figsize=(15,8))

pos = radial_distance_layout(T,keystr,mode='soph')
nx.draw_networkx(T,pos,ax=ax[0])
pos = radial_distance_layout(T,keystr,mode='normal')
nx.draw_networkx(T,pos,ax=ax[1])
pl.show()




