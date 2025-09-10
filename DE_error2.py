#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 14:29:25 2025

@author: edenye
"""
import matplotlib.pyplot as plt
import numpy as np
import math

plt.plot(())

a = 1
b = 2
N = 1000
h = (b-a)/(N-1) #dx

D = 10
omega = 10
fa = 1

def f(x):
    """ analytic solution
    returns the y value at a given x value"""
    return fa*np.cosh(np.sqrt(abs(omega))*(x-a))+(D/np.sqrt(abs(omega)))*np.sinh(np.sqrt(abs(omega))*(x-a))

x = np.arange(a, b+h, h)
y = f(x)
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
#show the analytic solution
"""plt.plot(x,y)
plt.show()"""

h_values=[] # list of all the dx values
n_val = [] # list of all the N values
j=0
p=[] # list to store arrays of points for the numeric solution, each array varying in length
analytic_points = [] # list to store arrays of points for the analytic solution, each array varying in length

while j<50:
    h = (b-a)/(N-1) #dx
    x = np.arange(a, b+h, h) #points from a to b (inclusive) in steps of h
    points = np.zeros_like(x) 
    
    points[0] = fa #initial conditions
    points[1] = D*h + points[0]
    
    for i in range(1, len(x)-1):
        points[i+1] = (2+omega*h**2)*points[i] - points[i-1] #finite differences
    #print('h', h, 'points', points)
    
    # appending to each list
    analytic_points.append(f(x))    
    p.append(points)   
    h_values.append(h)
    n_val.append(N)
    # change N
    N = N-0.05*N
    j+=1

# find the errors for each h value   
max_error = [] # list of the maximum absolute error at each N
for i in range(len(p)):
    array = p[i] #extract an array from the list
    analytic = analytic_points[i] #extract an array from the list
    error=[]   
    for k in range(len(array)): #iterate through the current array and find absolute difference
        dif = abs(array[k]-analytic[k])
        error.append(dif)
    error_max = np.max(error) # extract the max error
    max_error.append(error_max) # store the max error for this N

#finding the power
power=[] # list of each power
# error and dx of the starting N
e1= max_error[0]
h1= h_values[0]

for i in range(1, len(max_error)):
    e2=max_error[i]
    h2=h_values[i]
    powers = (math.log(e1/e2))/(math.log(h1/h2)) # calculate power
    power.append(powers) #store power

plt.xlabel('N')
plt.ylabel('error')
plt.plot(n_val, max_error) #plot error
plt.show()

plt.plot(power) #plot power
   
print(power) #output power



# print graphs at each N
"""m=0
for i in range(len(p)):
    plt.title(m)
    plt.plot(p[i], color='black', label='numeric')
    plt.plot(analytic_points[i], ':',color='red', markersize=3, label='analytic')
    plt.legend()
    plt.show()
    m+=1"""

# alternative way to see whats wrong
"""print(max_error)
e1= max_error[0]
h1= h_values[0]
e2=max_error[-1]
h2=h_values[-1]
log_error=math.log(e2/e1)

log_h = math.log(h2/h1)
powers=log_error/log_h

print(powers)"""
    
"""power=[]

for i in range(len(max_error)-1):
    e1 = max_error[i] #first node
    e2 = max_error[i+1] #next node
    h1 = h_values[i] # h value first node
    h2 = h_values[i+1] # h value next node
    
    e = abs(e1/e2)
    log_error = math.log(e)
    log_h = math.log(h1/h2)
    powers = log_error/log_h
        
    power.append(powers)

print('powers 2:', power)"""

"""
#MSE
error_sq = []
for i in range(len(points2)):
    error = abs(points2[i]-y[i])
    error_sq.append(error**2)
    
MSE = sum(error_sq)/len(points2)
print('MSE', MSE)"""






