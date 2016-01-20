# Radial Distance Layout

Generates a radial layout for trees whose nodes are associated with a distance to the root, similar to how it has been done in [1]. You can choose the basic method or a more sophisticated version which makes a more efficient use of space.

[1] *The Hidden Geometry of Complex, Network-Driven Contagion Phenomena*, D. Brockmann, D. Helbing, Science Vol. 342, Issue 6164, pp. 1337-1342 (2013)

## Install 

    $ sudo python setup.py install

## Example

    $ python example.py

or look here:

```
#!python
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

#The Tree has to be a DiGraph! The root is always the one with distance 0.
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

```