import numpy as np
import matplotlib.pyplot as plt
#Define analytic function
def ode_analytic(x,a,fa,D,omega):
    so=np.sqrt(omega)
    return (D/so)*np.sinh(so*(x-a))+fa*np.cosh(so*(x-a))
#********************************************************
#Set domain
#********************************************************
a=1.0
b=2.0
#********************************************************
#Set boundary conditions
#********************************************************
fa=1.0
D=10.0
omega=10.0
n=range(11,101,10)
nN=np.size(n)
dx=np.zeros(nN)
err=np.zeros(nN)
#********************************************************
#Cycle through discretisations
#********************************************************
for j in range(0,nN):
    #initialise solution and discretisation
    #dx[j]=(b-a)/(n[j]-1)
    #xi=np.arange(a,b,dx[j])
    xi,dx[j]=np.linspace(a,b,num=n[j],retstep=True)
    #Calculate numerical solution
    fnum=np.zeros(n[j]) #initialise numerical vector
    fnum[0]=fa
    for i in range(0,n[j]-1):
        if i==0:
            fnum[i+1]=0.5*fa*(2.0+omega*dx[j]**2)+D*dx[j]
        elif i>0:
            fnum[i+1]=fnum[i]*(2.0+omega*dx[j]**2)-fnum[i-1]
    #Calculate analytic solution
    fan=ode_analytic(xi,a,fa,D,omega)
    err[j]=max(np.abs(fnum-fan))
    if j==0:
        #Display solution
        plt.figure(1)
        plt.plot(xi,fan,'b')
        plt.plot(xi,fnum,'ro')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.grid()
#********************************************************
#calculate order
#********************************************************
em=err[nN-1]
dxm=dx[nN-1]
p=np.zeros(nN-1)
for i in range(0,nN-1):
    p[i]=np.log(em/err[i])/np.log(dxm/dx[i])
#********************************************************
#Display order with number of nodes
#********************************************************
plt.figure(2)
plt.plot(n[0:nN-1],p)
plt.xlabel('Number of nodes, n')
plt.ylabel('Calculated order, p')
plt.grid()
plt.show()
