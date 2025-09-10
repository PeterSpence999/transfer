import numpy as np
import matplotlib.pyplot as plt
#Define analytic function
def ode_analytic(x,a,fa,D,omega):
    so=np.sqrt(omega)
    return (D/so)*np.sinh(so*(x-a))+fa*np.cosh(so*(x-a))
#********************************************************
#Set domain
#********************************************************
a=1
b=2
#********************************************************
#Set boundary conditions
#********************************************************
fa=1
D=10
omega=10
n=range(100,1000,100)
dx=np.zeros(len(n))
err=np.zeros(len(n))
#********************************************************
#Cycle through discretisations
#********************************************************
for j in range(0,len(n)):
    #initialise solution and discretisation
    dx[j]=(b-a)/n[j]
    xi=np.arange(a,b,dx[j])
    #Calculate numerical solution
    f_num=np.zeros(n[j]) #initialise numerical vector
    f_num[0]=fa
    for i in range(0,n[j]-1):
        if i==0:
            f_num[i+1]=0.5*fa*(2.0+omega*dx[j]**2)+D*dx[j]
        elif i>0:
            f_num[i+1]=f_num[i]*(2.0+omega*dx[j]**2)-f_num[i-1]
    #Calculate analytic solution
    f_an=ode_analytic(xi,a,fa,D,omega)
    err[j]=max(np.abs(f_num-f_an))
    if j==len(n)-1:
        #Display solution
        plt.figure(1)
        plt.plot(xi,f_an,'b')
        plt.plot(xi,f_num,'ro')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.grid()
#********************************************************
#calculate order
#********************************************************
em=err[len(n)-1]
dxm=dx[len(n)-1]
p=np.zeros(len(n)-2)
for i in range(0,len(n)-2):
    p[i]=np.log(err[i]/em)/np.log(dx[i]/dxm)
#********************************************************
#Display order with number of nodes
#********************************************************
plt.figure(2)
plt.plot(n[0:len(n)-2],p)
plt.xlabel('Number of nodes, n')
plt.ylabel('Calculated order, p')
plt.grid()
plt.show()



