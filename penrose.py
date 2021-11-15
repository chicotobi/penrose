import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from math import sin,cos,atan2,pi
from matplotlib.patches import Arc

@dataclass
class Point:
  x: float
  y: float

@dataclass
class Poly:
  t: str
  p: [int]

@dataclass
class Vertex(Point):
  ts: str
  l: int
  r: int
  ps: [int]
  n: int
  cont: str

free = 0.6

def get_coords(poly):
  x = [vertices[i].x for i in poly.p]
  y = [vertices[i].y for i in poly.p]
  return x, y

def plot_poly(ax,poly):
  x, y = get_coords(poly)
  ax.plot(x+[x[0]],y+[y[0]],'k')

def plot_arcs(ax,poly):
    x, y = get_coords(poly)
    t1 = atan2(y[1] - y[0], x[1] - x[0]) / pi * 180
    t2 = atan2(y[3] - y[0], x[3] - x[0]) / pi * 180
    t3 = atan2(y[3] - y[2], x[3] - x[2]) / pi * 180
    t4 = atan2(y[1] - y[2], x[1] - x[2]) / pi * 180
    if poly.t == 'k':
      r12 = free
      c12 = 'green'
      c34 = 'red'
    elif poly.t == 'd':
      r12 = .5*(5**.5-3) + free
      c12 = 'red'
      c34 = 'green'
    r34 = 1 - free
    ax.add_patch(Arc((x[0], y[0]), 2*r12, 2*r12, theta1=t1, theta2=t2, color=c12))
    ax.add_patch(Arc((x[2], y[2]), 2*r34, 2*r34, theta1=t3, theta2=t4, color=c34))

def plot(polys,vertices):
  fig, ax = plt.subplots(1)
  for poly in polys:
    plot_poly(ax,poly)
    #plot_arcs(ax,poly)
  ax.axis('equal')
  ax.axis('off')
  plt.show()
  plt.pause(0.1)

def add_vertex_if_new(p_new):
  dists = [(p.x-p_new.x)**2+(p.y-p_new.y)**2 for p in vertices]
  idx = np.argmin(dists)
  if dists[idx] > 1e-4:
    # New vertex
    vertex_new = Vertex(x=p_new.x,y=p_new.y,ts='',l=-1,r=-1,ps=[],n=-1,cont='')
    vertices.append(vertex_new)
    idx = len(vertices)-1
    unfinished_vertices.append(idx)
  return idx

def update_vertex_info(l,idx,r,typ,idx_poly):
  if vertices[idx].ts == '':
    # This is a new vertex
    vertices[idx].ts = typ
    vertices[idx].l = l
    vertices[idx].r = r
    vertices[idx].ps = [idx_poly]
  elif vertices[idx].l == r:
    # This vertex was already there before
    vertices[idx].ts = vertices[idx].ts + typ
    vertices[idx].ps = vertices[idx].ps + [idx_poly]
    vertices[idx].l = l
  elif vertices[idx].r == l:
    # This vertex was already there before
    vertices[idx].ts = typ + vertices[idx].ts
    vertices[idx].ps = [idx_poly] + vertices[idx].ps
    vertices[idx].r = r

def add_poly_at_vertex(idx):
  ptype = vertices[idx].cont[1]
  orig1 = int(vertices[idx].cont[0])

  if ptype == 'k':
    l = pts_kite
  elif ptype=='d':
    l = pts_dart
  a = l[orig1]
  b = l[(orig1+1)%4]
  c = vertices[idx]
  d = vertices[vertices[idx].l]

  phi = atan2(d.y - c.y, d.x - c.x) - atan2(b.y - a.y, b.x - a.x)

  def my_trafo(p):
    x1 = cos(phi) * (p.x-a.x) - sin(phi) * (p.y-a.y) + c.x
    y1 = sin(phi) * (p.x-a.x) + cos(phi) * (p.y-a.y) + c.y
    return Point(x=x1,y=y1)

  tmp = [add_vertex_if_new(my_trafo(i)) for i in l]

  polys.append(Poly(ptype,tmp))
  idx_poly = len(polys)-1

  for i in range(4):
    update_vertex_info(tmp[(i-1)%4],tmp[i],tmp[(i+1)%4],str(i)+ptype,idx_poly)

  return tmp

def update_possible_continuations(idx):
  s = vertices[idx].ts
  vertices[idx].n = 0
  if s in possible_vertices:
    vertices[idx].cont = ''
    unfinished_vertices.remove(idx)
    return
  k = len(s)
  for pv in possible_vertices:
    if s == pv[:k]:
      vertices[idx].n += 1
      vertices[idx].cont = pv[k:k+2]

# Data for kites and darts
x1 = .25*(1+5**.5)
x2 = .5*(1+5**.5)
y1 = .5*(.5*(5-5**.5))**.5
pts_kite = [Point(x=0,y=0), Point(x=x1,y=-y1), Point(x=1 ,y=0), Point(x=x1,y=y1)]
pts_dart = [Point(x=1,y=0), Point(x=x1,y=-y1), Point(x=x2,y=0), Point(x=x1,y=y1)]

# Possible connections encoded
possible_vertices = [
  '2d2d2d3k1k', '2d2d3k1k2d', '2d3k1k2d2d', '3k1k2d2d2d', '1k2d2d2d3k',
  '3d0k0k1d2k', '0k0k1d2k3d', '0k1d2k3d0k', '1d2k3d0k0k', '2k3d0k0k1d',
  '3k1k3k1k2d', '1k3k1k2d3k', '3k1k2d3k1k', '1k2d3k1k3k', '2d3k1k3k1k',
  '2k2k3d1d', '2k3d1d2k', '3d1d2k2k', '1d2k2k3d',
  '0d1k3k', '1k3k0d', '3k0d1k',
  '2d2d2d2d2d',
  '0k0k0k0k0k']

# Initial
vertices = [Vertex(x=p.x,y=p.y,ts='',l=-1,r=-1,ps=[],n=-1,cont='') for p in pts_kite]
polys = [Poly('k',[0,1,2,3])]
for i in range(4):
  update_vertex_info((i-1)%4,i,(i+1)%4,str(i)+'k',0)
unfinished_vertices = [0,1,2,3]
changed_vertices = [0,1,2,3]

import time
t_calc = 0
t_plt = 0

# Check if there is a node which is determined
for j in range(5000):
  t = time.time()
  for idx in changed_vertices:
    update_possible_continuations(idx)
  for idx in unfinished_vertices:
    if vertices[idx].n == 1:
      break
  changed_vertices = add_poly_at_vertex(idx)
  t_calc += time.time() - t
  if j%1000==0:
    t = time.time()
    plot(polys,vertices)
    t_plt += time.time() - t
    print(j,t_calc,t_plt)
    t_calc = 0
    t_plt = 0