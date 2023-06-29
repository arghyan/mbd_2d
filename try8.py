from sympy.physics.mechanics import ReferenceFrame
from kinmod4 import *
from sympy import symbols
import numpy as np

val=np.loadtxt('pos.txt')
aval=np.loadtxt('lin_acc.txt')
angavalt=np.loadtxt('ang_acc.txt')

N=ReferenceFrame('N')

Ox,Oy,P1x,P1y,P2x,P2y=symbols('Ox,Oy,P1x,P1y,P2x,P2y')
m1=36
m2=36
L=1
Ig1=(1/12)*m1*L**2
Ig2=(1/12)*m2*L**2
g=10

# Time step
i=0

[m,n]=np.shape(val)
n1=np.floor_divide(n,2)
O=ourpoint('O',N,[val[i,0],val[i,0+n1]])
G1=ourpoint('G1',N,[val[i,1],val[i,1+n1]])
P1a=ourpoint('P1a',N,[val[i,2],val[i,2+n1]])
P1b=ourpoint('P1b',N,[val[i,2],val[i,2+n1]])
G2=ourpoint('G2',N,[val[i,3],val[i,3+n1]])
P2=ourpoint('P2',N,[val[i,4],val[i,4+n1]])

O.setfor_xy(N,[Ox,Oy])
P1a.setfor_xy(N,[P1x,P1y])
P1b.setfor_xy(N,[-P1x,-P1y])
P2.setfor_xy(N,[P2x,P2y])

G1.setfor_xy(N,[0,-m1*g])
G2.setfor_xy(N,[0,-m2*g])

Link1=ourbody('OP1',[O,G1,P1a],G1)
Link2=ourbody('P1P2',[P1b,G2,P2],G2)

Link1.set_for_dynamics(G1,m1,G1,Ig1,N,[aval[i,0],aval[i,2]],angavalt[i,0])
Link2.set_for_dynamics(G2,m2,G2,Ig2,N,[aval[i,1],aval[i,3]],angavalt[i,1])

Link1.generate_eqn(N)
Link2.generate_eqn(N)

group1 = group_of_bodies('group1',[Link1,Link2])
group1.gr_solve([Ox,Oy,P1x,P1y,P2x,P2y])