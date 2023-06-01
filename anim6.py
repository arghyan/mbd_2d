import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def our_rotate(theta,xvalt,yvalt):
    m,n1=np.shape(xvalt)
    for i in range(m):
        for j in range(n1):
            xdum=xvalt[i,j]
            ydum=yvalt[i,j]
            xvalt[i,j]=np.cos(theta)*xdum+np.sin(theta)*ydum
            yvalt[i,j]=-np.sin(theta)*xdum+np.cos(theta)*ydum
    return xvalt,yvalt

def our_reflect(refl_abt,xvalt,yvalt):
    m,n1=np.shape(xvalt)
    for i in range(m):
        for j in range(n1):
            xdum=xvalt[i,j]
            ydum=yvalt[i,j]
            if refl_abt=='x':
                yvalt[i,j]=-ydum
            if refl_abt=='y':
                xvalt[i,j]=-xdum
    return xvalt,yvalt

val=np.loadtxt('db_check_anim2.txt')

m,n=np.shape(val)
n1=n//2
xvalt=val[:,0:n1]
yvalt=val[:,n1:n]

#xvalt,yvalt= our_rotate(np.pi/2,xvalt,yvalt)
xvalt,yvalt= our_reflect('y',xvalt,yvalt)

xmax=np.max(xvalt)
xmin=np.min(xvalt)
ymax=np.max(yvalt)
ymin=np.min(yvalt)
ext=max(xmax-xmin,ymax-ymin)
ext=ext*1.1
xc=(xmin+xmax)/2
yc=(ymin+ymax)/2

llist=[[0,1,2]]

def init():
    for dumline in lines:
        dumlines.set_data([], [])
    return lines

def our_animate1(ii):
    for i in range(0,len(lines)):
        lines[i].set_data(xvalt[ii,llist[i]],yvalt[ii,llist[i]])
    return lines

fig = plt.figure()
axis = plt.axes(xlim =(xc-ext/2, xc+ext/2),
                ylim =(yc-ext/2, yc+ext/2))
 
line1, = plt.plot([], [], lw = 5)
lines=[line1]

#line=our_animate1(2)
#plt.show()

anim = animation.FuncAnimation(fig, our_animate1,                         
                            frames = m,
                            interval = m,
                            blit = True)

plt.show()