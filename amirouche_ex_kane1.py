from sympy.physics.mechanics import *
from sympy import trigsimp,symbols, sin, cos, shape, lambdify
from math import pi
init_vprinting(pretty_print=False)

q0,q1,q2 = dynamicsymbols('q0 q1,q2')
q0d,q1d,q2d = dynamicsymbols('q0 q1 q2',1)
u0,u1,u2 = dynamicsymbols('u0,u1,u2')


m0,m,I1,I2,L,g=symbols('m0,m,I1,I2,L,g')
R,S=symbols('R,S')

N = ReferenceFrame('N')
A = N.orientnew('A','Axis',[0,N.z])
O=Point('O')
B1=N.orientnew('B1','Axis',[-q1, N.z])
P1 = O.locatenew('P1',L*B1.y)
G1 = O.locatenew('G1',(L/2)*B1.y)
B2=B1.orientnew('B2','Axis',[-q2, B1.z])
P2 = P1.locatenew('P2',L*B2.y)
G2 = P1.locatenew('G2',(L/2)*B2.y)
P = O.locatenew('P',q0*N.x)

P.set_vel(N,u0*N.x)
A.set_ang_vel(N,0)
B1.set_ang_vel(N,-u1*N.z)
B2.set_ang_vel(N,-(u1+u2)*N.z)
O.set_vel(N,0)
P1.v2pt_theory(O,N,B1)
G1.v2pt_theory(O,N,B1)
P2.v2pt_theory(P1,N,B2)
G2.v2pt_theory(P1,N,B2)

kde=[u0-q0d,u1-q1d,u2-q2d]

cons1=L*sin(q1)+L*sin(q1+q2)-q0
cons2=L*cos(q1)+L*cos(q1+q2)
cons=[cons1,cons2]
dcons=[time_derivative(cons1,N), time_derivative(cons2,N)]
coordinates=[q0,q1,q2]
speeds=[u0,u1,u2]


kane= KanesMethod(N,q_ind=[q0],u_ind=[u0],q_dependent=[q1,q2],u_dependent=[u1,u2],
                 configuration_constraints=cons, velocity_constraints=dcons,kd_eqs=kde)



it1=inertia(B1,0,0,m*L**2/12) 
it2=inertia(B2,0,0,m*L**2/12)
i1=(it1,G1)
i2=(it2,G2)

par=Particle('particle',P,m0)
link1 = RigidBody('link1', G1, B1, m, i1)
link2 = RigidBody('link2', G2, B2, m, i2)

Bodies=[par,link1,link2]
Forces=[(G1,-m*g*N.y),(G2,-m*g*N.y),(B1,-(R+S)*N.z),(B2,S*N.z)]

fr,frstar = kane.kanes_equations(Bodies,Forces)
mm_cons= kane.mass_matrix
ff_cons= kane.forcing

subs_dict={L:1,m0:1,m:1,g:10,R:1,S:0}

mm_halfq=mm_cons.subs(subs_dict)
mm_half=lambdify([q1,q2],mm_halfq)

ff_halfq=ff_cons.subs(subs_dict)
ff_half=lambdify([q0,q1,q2,u0,u1,u2],ff_halfq)


