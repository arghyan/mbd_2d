from sympy.physics.mechanics import *
from sympy import trigsimp,symbols,cos,sin,lambdify
from math import radians,degrees
import numpy as np
from scipy import linalg, integrate

q1,q2 = dynamicsymbols('q1,q2')
#q1d,q2d = dynamicsymbols('q1 q2',1)
u1,u2 = dynamicsymbols('q1,q2',1)
u1d,u2d = dynamicsymbols('q1,q2',2)
t = dynamicsymbols._t

m1,m2,L1,L2,g=symbols('m1,m2,L1,L2,g')
R,S=symbols('R,S')
a1=symbols('a1')

N = ReferenceFrame('N')
O=Point('O')
B1=N.orientnew('B1','Axis',[q1, N.z])
P1 = O.locatenew('P1',L1*B1.x)
G1 = O.locatenew('G1',(L1/2)*B1.x)
B2=B1.orientnew('B2','Axis',[q2, B1.z])
P2 = P1.locatenew('P2',L2*B2.x)
G2 = P1.locatenew('G2',(L2/2)*B2.x)

B1.set_ang_vel(N,u1*N.z)
B2.set_ang_vel(N,(u1+u2)*N.z)
O.set_vel(N,0)
P1.v2pt_theory(O,N,B1)
G1.v2pt_theory(O,N,B1)
P2.v2pt_theory(P1,N,B2)
G2.v2pt_theory(P1,N,B2)

Tt=(1/2)*m1*dot(G1.vel(N),G1.vel(N)) +(1/2)*m2*dot(G2.vel(N),G2.vel(N))
Tr=(1/2)*(m1*L1**2/12)*dot(B1.ang_vel_in(N),B1.ang_vel_in(N)) + (1/2)*(m2*L2**2/12)*dot(B2.ang_vel_in(N),B2.ang_vel_in(N))
T=Tt+Tr
#V1=m*g*dot(G1.pos_from(O),N.y)
#V2=m*g*dot(G2.pos_from(O),N.y)
#V=-R*q1+S*q2
V=0

cons1=L1*cos(q1)*u1+L2*cos(q1+q2)*(u1+u2)-2*a1*t


Forces=[(G1,-m1*g*N.y),(G2,-m2*g*N.y)]
#Forces=[(G1,-m*g*N.y),(G2,-m*g*N.y)]
L=T-V
LM = LagrangesMethod(L, [q1,q2], forcelist = Forces, nonhol_coneqs=[cons1], frame=N)
#LM = LagrangesMethod(L, [q1,q2], forcelist = Forces,frame=N)
LM.form_lagranges_equations()


trigsimp(LM.mass_matrix)
trigsimp(LM.forcing)


V

mm=LM.mass_matrix_full
ff=LM.forcing_full

subs_dict={L1:1,L2:1,m1:36,m2:36,g:10,a1:1}
mmq=mm.subs(subs_dict)
forcingq=ff.subs(subs_dict)

mmq
forcingq

mmf=lambdify([q1,q2],mmq)
forcingf=lambdify([q1,q2,u1,u2],forcingq)

mmf(1,1)

forcingf(1,1,0.5,0.5)

def dydt(stv,t):
    amat=mmf(stv[0],stv[1])
    bvec=forcingf(stv[0],stv[1],stv[2],stv[3])
    xvec=linalg.solve(amat,bvec)
    return xvec.flatten()

stv_in=np.array([radians(70),radians(-140),0,0,0],dtype=float)
nval=9; ts = np.linspace(0, 0.8, nval)

out=integrate.odeint(dydt,stv_in,ts)

out;

qlist=[q1,q2]
plist=[O,G1,P1,G2,P2]
slist=out

xvalt=np.array([],float)
yvalt=np.array([],float)
for ii in range(0,np.shape(slist)[0]):
    xval=np.array([],float)
    yval=np.array([],float)
    for pdum in plist:
        #print(pdum)
        posv=pdum.pos_from(O).subs(subs_dict)
        #print(posv)
        for i in range(0,len(qlist)):
            posv=posv.express(N).subs({qlist[i]:slist[ii,i]})
        posvx=dot(posv,N.x).evalf()
        posvy=dot(posv,N.y).evalf()
        #print(posvx)
        xval=np.append(xval,posvx)
        yval=np.append(yval,posvy)
    xvalt=np.append(xvalt,xval)
    yvalt=np.append(yvalt,yval)
xvalt=np.reshape(xvalt,(len(slist),len(plist)))
yvalt=np.reshape(yvalt,(len(slist),len(plist)))
val=np.block([xvalt,yvalt])

ts
np.reshape(ts,(9,1))
ts
val[:,5];

np.savetxt('pos.txt',val)


stv_deriv=np.array([],float)
for ii in range(0,np.shape(slist)[0]):
    amat=mmf(slist[ii,0],slist[ii,1])
    bvec=forcingf(slist[ii,0],slist[ii,1],slist[ii,2],slist[ii,3])
    xvec=linalg.solve(amat,bvec)
    stv_deriv=np.append(stv_deriv,xvec)
stv_deriv=np.reshape(stv_deriv,(len(slist),np.shape(slist)[1]))

stv_deriv;

np.savetxt('stv_deriv.txt',stv_deriv)

P1.a2pt_theory(O,N,B1)
G1.a2pt_theory(O,N,B1)
G2.a2pt_theory(P1,N,B2);

cmlist=[G1,G2]
qlist=[q1,q2]
qdlist=[u1,u2,u1d,u2d]

axvalt=np.array([],float)
ayvalt=np.array([],float)
for ii in range(0,np.shape(slist)[0]):
    axval=np.array([],float)
    ayval=np.array([],float)
    qdict=subs_dict
    for i in range(0,len(qlist)):
        qdict.update({qlist[i]:slist[ii,i]})
    for i in range(0,len(qdlist)):
        qdict.update({qdlist[i]:stv_deriv[ii,i]})        
    for cmdum in cmlist:
        #print(cmdum)
        accv=cmdum.acc(N).express(N).subs(qdict)
        #print(accv)
        accvx=dot(accv,N.x).evalf()
        accvy=dot(accv,N.y).evalf()
        axval=np.append(axval,accvx)
        ayval=np.append(ayval,accvy)
    axvalt=np.append(axvalt,axval)
    ayvalt=np.append(ayvalt,ayval)
axvalt=np.reshape(axvalt,(len(slist),len(cmlist)))
ayvalt=np.reshape(ayvalt,(len(slist),len(cmlist)))
aval=np.block([axvalt,ayvalt])

aval;

np.savetxt('lin_acc.txt',aval)

B2.ang_acc_in(N)

bdlist=[B1,B2]
qdlist=[u1,u2,u1d,u2d]
angavalt=np.array([],float)
for ii in range(0,np.shape(slist)[0]):
    angaval=np.array([],float)
    qdict={}
    for i in range(0,len(qdlist)):
        qdict.update({qdlist[i]:stv_deriv[ii,i]})        
    for bddum in bdlist:
        angav=bddum.ang_acc_in(N).subs(qdict)
        angavz=dot(angav,N.z).evalf()
        angaval=np.append(angaval,angavz)
    angavalt=np.append(angavalt,angaval)
angavalt=np.reshape(angavalt,(len(slist),len(bdlist)))


angavalt;


np.savetxt('ang_acc.txt',angavalt)










