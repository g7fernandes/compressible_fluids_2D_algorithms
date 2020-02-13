#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 11:56:01 2018

@author: gabriel
"""
import numpy as np
import matplotlib.pyplot as plt
from math import pi

def hs(a):
    a = a>=0
    return a

caso = input('Caso: ')

#%% define os pontos da malha

if caso == '1':
    xi = np.linspace(0,1,72)
    eta = np.linspace(0,1,40)
    
    x = lambda xi,eta: xi*3.4-.6
    y = lambda xi,eta: eta*(1.8-np.tan(10*pi/180)*hs(x(xi,eta))*x(xi,eta))+ \
    np.tan(10*pi/180)*hs(x(xi,eta))*x(xi,eta)
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha1.vtk"
    
elif caso == '2':
    theta = 39.31*pi/180
    
    xi = np.linspace(0,1,60)
    eta = np.linspace(0,1,40)
    
    x = lambda xi,eta: xi*3.4-0.6 + eta*1.8/np.tan(theta)
    y = lambda xi,eta: eta*(1.8-np.tan(10*pi/180)*hs((x(xi,eta)- \
    eta*1.8/np.tan(theta)))*(x(xi,eta)- eta*1.8/np.tan(theta)))+ \
    np.tan(10*pi/180)*hs((x(xi,eta)- eta*1.8/np.tan(theta)))*(x(xi,eta)- \
    eta*1.8/np.tan(theta))        
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
#    for k in range(len(xi)):
#        plt.plot(X[:,k],Y[:,k],'b')
#    for k in range(len(eta)):
#        plt.plot(X[k,:],Y[k,:],'b')
#    plt.show()       
    file_name = "malha2p1.vtk"
    
    
elif caso == '3':    
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,10)

    x = lambda xi,eta: xi*100
    y = lambda xi,eta: eta*10    
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)

    file_name = "malha3.vtk"
    
elif caso == '4':    
    xi = np.linspace(0,1,20)
    eta = np.linspace(0,1,50)

    x = lambda xi,eta: xi*10
    y = lambda xi,eta: eta*100    
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)

    file_name = "malha5.vtk"

elif caso == '5':          
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,10)  

    x = lambda xi,eta: xi*100
    y = lambda xi,eta: eta*20+xi*100*np.tan(10*pi/180) 
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha6.vtk"    
    
elif caso == '6':          
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,10)  

    x = lambda xi,eta: xi*100
    y = lambda xi,eta: eta*(10+xi*10)
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha7.vtk"      
    
elif caso == '7':          
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,10)  

    x = lambda xi,eta: xi*100
    y = lambda xi,eta: 10-eta*(20-xi*10)
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha8.vtk"     
    
elif caso == '8':
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,10)

    x = lambda xi,eta: xi*100
    y = lambda xi,eta: eta*20    
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    for a in range(40):
        for b in range(2,8):
            Y[b,a] = (Y[b,a])-(Y[b,a]*3)/((50-a))

    file_name = "malha9.vtk"        
    
elif caso == '9':
    theta = 39.31*pi/180
        
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,10)

    x = lambda xi,eta: xi*10-0.6 + eta*1.8/np.tan(theta)
    y = lambda xi,eta: eta*20    
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)

    file_name = "malha10.vtk"      
    
elif caso == '10':
    xi = np.linspace(0,1,60)
    eta = np.linspace(0,1,30)
    
    x = lambda xi,eta: xi*4-.6
    y = lambda xi,eta: (eta*(1.8-np.tan(10*pi/180)*hs(x(xi,eta))*x(xi,eta))+ \
    np.tan(10*pi/180)*hs(x(xi,eta))*x(xi,eta))*(hs(2.8-x(xi,eta))) + \
            hs(x(xi,eta)-2.8)*(eta*(1.8-np.tan(10*pi/180)*2.8)+np.tan(10*pi/180)*2.8)            
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha11.vtk"      
    
    
if caso == '11':
    xi = np.linspace(0,1,36)
    eta = np.linspace(0,1,20)
    
    x = lambda xi,eta: xi*3.4
    y = lambda xi,eta: eta*(1.8-np.tan(10*pi/180)*x(xi,eta))+ \
    np.tan(10*pi/180)*x(xi,eta)
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha12.vtk"    
    
elif caso == '13':    

      
    xi = np.linspace(0,1,40)
    eta = np.linspace(0,1,20)  

    x = lambda xi,eta: xi*3.4
    y = lambda xi,eta: eta*1.8+xi*3.4*np.tan(10*pi/180) 
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha13.vtk"         
    
if caso == '14':
    xi = np.linspace(0,1,36)
    eta = np.linspace(0,1,20)
    
    x = lambda xi,eta: xi*3.4*20-.6*20
    y = lambda xi,eta: eta*(1.8*20-np.tan(10*pi/180)*hs(x(xi,eta))*x(xi,eta))+ \
    np.tan(10*pi/180)*hs(x(xi,eta))*x(xi,eta)
    
    Xi,Eta = np.meshgrid(xi,eta)
    
    X = x(Xi,Eta)
    Y = y(Xi,Eta)
    
    file_name = "malha14.vtk"    
#%% faz malha     

celulas = np.array([[4,0,0,0,0]]*((len(xi)-1)*(len(eta)-1)))
pontos = np.zeros(((len(xi)*len(eta)),3))
npontos = 0;
c = 0;

jlen = len(xi)
ilen = len(eta)

for i in range(ilen):
    for j in range(jlen):
        pontos[npontos,:] = np.array([X[i,j],Y[i,j],0])
        npontos += 1

for i in range(ilen-1):
    for j in range(jlen-1):
         celulas[c,1:5] = np.array([j+i*jlen,j+1+i*jlen,j+jlen+1+i*jlen,j+jlen+i*jlen])
         c +=1

#for i in range(len(eta)-1):
#    for j in range(len(xi)-1):
#        if j == 0:
#            if i == 0:
#                pontos[npontos,:] = np.array([X[i,j],Y[i,j],0])
#                npontos+=1                       
#            pontos[npontos+3,:] = np.array([X[i+1,j],Y[i+1,j],0])
#            npontos += 1
#        if i == 0:        
#            pontos[npontos+1,:] = np.array([X[i,j+1],Y[i,j+1],0])
#            npontos += 1            
#        pontos[npontos+2,:] = np.array([X[i+1,j+1],Y[i+1,j+1],0])
#        npontos += 1;
#        
#        celulas[c,1:5] = np.array([npontos,(npontos+1),(npontos+2),(npontos+3)]) 
#        c += 1        

for i in pontos:
    plt.plot(i[0],i[1],'k*')
plt.show()        

for k in celulas:
    L = np.array([pontos[k[1]],pontos[k[2]],pontos[k[3]],pontos[k[4]],pontos[k[1]]])
    plt.plot(L[:,0],L[:,1],'b')
plt.show()
#%% escreve arquivo vtk
    
file = open(file_name,'w')   

file.write("# vtk DataFile Version 3.0\n")
file.write("malha\n")           
file.write("ASCII\n")
file.write("DATASET UNSTRUCTURED_GRID\n")
file.write("POINTS {} float\n".format(npontos)) 
for pto in pontos:
    file.write("{} {} {}\n".format(pto[0],pto[1],pto[2]))
file.write("CELLS {} {}\n".format(c,5*c)) # 
for cel in celulas:
    file.write("{} {} {} {} {}\n".format(cel[0],cel[4],cel[3],cel[2],cel[1]))
file.write("CELL_TYPES {}\n".format(c))
for i in range(c):
    file.write("{}\n".format(9))  # INDICA CELULA QUADRADA vtk_quad

file.close()
print("Malha {} salva!".format(file_name))
    








#vtk output
#ASCII
#DATASET UNSTRUCTURED_GRID
#POINTS 4 float
#0 0 0
#1 0 0
#0 1 0
#1.1 1.1 0
#CELLS 1 5
#4 0 1 3 2
#CELL_TYPES 1
#9
#CELL_DATA 1
#POINT_DATA 4
#FIELD FieldData 1
#nodal 1 4 float
#0 1 1.1 2
#















