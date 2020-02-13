# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 19:48:00 2018

@author: Gabriel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 10:04:47 2018



Este programa lê o arquivo VTK. E armazena os pontos nas linhas da variável 
'pontos' e as células dos volumes na variável 'celulas'. Cada linha em 'celula' 
corresponde a um volume formado pelos pontos indicados nas colunas das linhas 
ligados no sentido horário. Os pontos nas linhas da variável 'célula' cores-
pondem aos índices na variável 'pontos'.

Os contornos são indicados. 

Lê-se um arquivo BC.txt com as condições de contorno nos contornos. 

Invoca-se o solver. 


@author: gabriel

"""

 

#Pre-processamento
#serve pra malha nao estruturada

# identificar contorno e vizinhos

import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from evtk.hl import unstructuredGridToVTK
from evtk.vtk import VtkQuad

def exportVTK(pontos,celulas,info,arquivo):
# info deve ser algum dado em formato array(ex. pressão, densidade, etc) e 
# nome_info string com seu rótulo (ex. "pressure")
 
    conn = np.empty(len(celulas)*4)
    ctype = np.ones(len(celulas))*VtkQuad.tid
    conn = np.ravel(celulas)
    offset = np.array([i*4 for i in range(1,len(celulas)+1)])
    unstructuredGridToVTK(arquivo, np.ravel(pontos[:,0]), np.ravel(pontos[:,1]),np.ravel(pontos[:,2]), \
                          connectivity = conn, offsets = offset, cell_types = ctype, \
                          cellData = info, pointData = None)
   
def PolygonArea(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def estacheia(A):
    if A[0] == -1:
        return 0
    elif A[1] == -1:
        return 0
    elif A[2] == -1:
        return 0
    elif A[3] == -1:
        return 0
    else:
        return 1

#%% le vtk
file_name = input('Entre a malha\n')

count = 0;
with open(file_name) as file:
    l = file.readline().split()
    while l:
        if l[0] == 'POINTS': # coleta os pontos
            npoints = int(l[1])
            pontos = np.zeros((npoints,3))
            while count < npoints: #
                l = file.readline().split()
           
                pontos[count,:] = np.array([float(i) for i in l])
                count += 1
            count = 0    

        if l[0] == 'CELLS':
            ncells = int(l[1])
            celulas = np.zeros((ncells,4),dtype=np.int)
            while count < ncells:
                l = file.readline().split()               
                celulas[count,:] = np.array([float(i) for i in l][1:5])
                count += 1
            count = 0      
        l = file.readline().split()            
                    
             
#%% encontra vizinhos

viz = np.ones((ncells,4),dtype=int)*-1
Vol = np.empty((ncells,1))
SR = np.empty((ncells,2))
SL= np.empty((ncells,2))
SU = np.empty((ncells,2))
SD = np.empty((ncells,2))  #superficies em volta do elemento
fantasma = [] #nao necessario mas interessante ver quais sao as celulas no contorno
elementos = []
Br = []
Bl = []
Bu = []
Bd = []

for cact in range(ncells):
    cel = celulas[cact,:]
    # Volumes
    Vol[cact] = PolygonArea([(pontos[cel[0],0],pontos[cel[0],1]), \
                     (pontos[cel[1],0],pontos[cel[1],1]), \
                     (pontos[cel[2],0],pontos[cel[2],1]), \
                     (pontos[cel[3],0],pontos[cel[3],1])])
    # superfícies
    SU[cact,:] = np.array([(pontos[cel[1],0]-pontos[cel[0],0]) , (pontos[cel[1],1]-pontos[cel[0],1])])
    SR[cact,:] = np.array([(pontos[cel[1],0]-pontos[cel[2],0]) , (pontos[cel[1],1]-pontos[cel[2],1])])
    SD[cact,:] = np.array([(pontos[cel[2],0]-pontos[cel[3],0]) , (pontos[cel[2],1]-pontos[cel[3],1])])
    SL[cact,:] = np.array([(pontos[cel[0],0]-pontos[cel[3],0]) , (pontos[cel[0],1]-pontos[cel[3],1])])
#    S[cact,:] = np.array([np.sqrt((pontos[cel[0],0]-pontos[cel[1],0])**2+(pontos[cel[0],1]-pontos[cel[1],1])**2), \
#                         np.sqrt((pontos[cel[1],0]-pontos[cel[2],0])**2+(pontos[cel[1],1]-pontos[cel[2],1])**2), \
#                         np.sqrt((pontos[cel[2],0]-pontos[cel[3],0])**2+(pontos[cel[2],1]-pontos[cel[3],1])**2), \
#                         np.sqrt((pontos[cel[3],0]-pontos[cel[0],0])**2+(pontos[cel[3],1]-pontos[cel[0],1])**2)])    # cima, direita, baixo esquerda



    # sao vizinhos se compartilharem uma aresta
    edgeY2 = np.array([cel[1],cel[0]])
    edgeX2 = np.array([cel[2],cel[1]])
    edgeY1 = np.array([cel[3],cel[2]])
    edgeX1 = np.array([cel[0],cel[3]])
    

    for cviz in range(cact+1,ncells):
#        print('{} {}'.format(cact,cviz))
#        print('{}'.format(viz[cact,:]))
#        time.sleep(.25)
        if estacheia(viz[cact,:]):
#            print('ch')
            break #ja encontrou todos os vizinhos
        
        if (edgeY2 == celulas[cviz,2:4]).all():
            viz[cact,0] = cviz
            viz[cviz,2] = cact
#            print('a')
        elif (edgeY1 == celulas[cviz,0:2]).all():
            viz[cact,2] = cviz
            viz[cviz,0] = cact
#            #print('b')            
        elif (edgeX1 == celulas[cviz,1:3]).all():
            viz[cact,3] = cviz
            viz[cviz,1] = cact
#            print('c')            
        elif (edgeX2 == np.array([celulas[cviz,3],celulas[cviz,0]])).all():
            viz[cact,1] = cviz
            viz[cviz,3] = cact
#            print('{} vz {}'.format(viz[cact],viz[cviz]))            
#            print('d')                          
    
    if estacheia(viz[cact]):
        elementos.append(cact)  
    else:
        fantasma.append(cact)           
#        print('fan {}'.format(viz[cacy]))  

dx_min = np.sqrt(min(Vol)) 


for i in elementos:
    k = celulas[i,:]
    L = np.array([pontos[k[0]],pontos[k[1]],pontos[k[2]],pontos[k[3]],pontos[k[0]]])
    plt.plot(L[:,0],L[:,1],'k')    
    
for i in fantasma:
    k = celulas[i,:]
    L = np.array([pontos[k[0]],pontos[k[1]],pontos[k[2]],pontos[k[3]],pontos[k[0]]])
    if (viz[i,0]) == -1 and (viz[i,3]) != -1 and (viz[i,1]) != -1:
        Bu.append(i) #cima
        cor = 'b'
    elif (viz[i,1]) == -1 and (viz[i,0] >= 0 and viz[i,2] >= 0):
        Br.append(i) # direita
        cor = 'r'
    elif (viz[i,2]) == -1 and (viz[i,1] >= 0 and viz[i,3] >= 0):
        Bd.append(i)  # baixo
        cor = 'g'
    elif (viz[i,3]) == -1 and (viz[i,0] >= 0 and viz[i,2] >= 0):
        Bl.append(i)  # esquerda
        cor = 'c'
    else:
        continue
        
        
        
    plt.plot(L[:,0],L[:,1],cor)   
green_patch = mpatches.Patch(color='green', label='B1')
blue_patch = mpatches.Patch(color='blue', label='B2')
cyan_patch = mpatches.Patch(color='cyan', label='B3')
red_patch = mpatches.Patch(color='red', label='B4')
plt.legend(handles=[green_patch,blue_patch,cyan_patch,red_patch])    
plt.title('Lista de contornos')                     
plt.show()                
print('Entrada: B3\n Saida:B4,Paredes: B1, B2\n')
print('Lendo condiçoes iniciais e de contorno do arquivo BC.txt ...')



#%% Condiçao de Contorno




icss = 0 # ic supersônica
inss = 0 # entrada supersônica

with open('BC.txt') as bc:
    aux = True 
    line = bc.readline()
    while(aux):     
        line = bc.readline()
        l = line.split()
        if l[0] == 'PRESSURE':
            if l[1] == 'IC':
                line = bc.readline()
                pIC = float(line)
            elif l[1] == 'ENTRADA':
                line = bc.readline()
                pIN = float(line)
            elif l[1] == 'SAIDA':
                line = bc.readline()
                pOUT = float(line)
        elif l[0] == 'VELOCITY_U':
            if l[1] == 'IC':
                line = bc.readline().split()
                if line[0] == 'c':
                    icss = 1
                    uIC = float(line[1])
                else:
                    uIC = float(line[0])
            elif l[1] == 'ENTRADA':
                line = bc.readline().split()
                if line[0] == 'c':
                    inss = 1
                    uIN = float(line[1])
                else:
                    uIN = float(line[0])            
        elif l[0] == 'VELOCITY_V':
            
            if l[1] == 'IC':
                line = bc.readline().split()
                if line[0] == 'c':
                    icss = 1
                    vIC = float(line[1])
                else:
                    vIC = float(line[0])
            elif l[1] == 'ENTRADA':
                line = bc.readline().split()
                if line[0] == 'c':
                    inss = 1
                    vIN = float(line[1])
                else:
                    vIN = float(line[0])   
        elif l[0] == 'DENSITY':
            if l[1] == 'IC':
                line = bc.readline()
                rhoIC = float(line)
            elif l[1] == 'ENTRADA':
                line = bc.readline()
                rhoIN = float(line)                            
        elif l[0] == 'TEMPERATURE':
            if l[1] == 'IC':
                line = bc.readline()
                T_IC = float(line)
            elif l[1] == 'ENTRADA':
                line = bc.readline()
                T_IN = float(line)                      
        elif l[0] == 'k':        
            line = bc.readline()
            k = float(line)
        elif l[0] == 'CFL':
            line = bc.readline()
            cfl = float(line)
        elif l[0] == '%':
            aux = False
             
if not ('rhoIC' in locals()):
    print('Densidade inicial não especificada,\n será utilizado do ar como gás perfeito')
    rhoIC = pIC*.028947/(8.314*T_IC)
if not ('rhoIN' in locals()):
    print('Densidade da entrada não especificada,\n será utilizado do ar como gás perfeito')
    rhoIN = pIN*.028947/(8.314*T_IN)    

cIN = np.sqrt(k*pIN/rhoIN) #velocidade do som na entrada
cIC = np.sqrt(k*pIC/rhoIC) #velocidade do som de IC

if icss:
    vIC = cIN*vIC
    uIC = cIN*uIC
if inss:
    vIN = cIN*vIN
    uIN = cIN*uIN        
if not inss or max([vIN,uIN]) < cIN:
    print('condição de bocal subsônico')
    if not ('pOUT' in locals()):
        print('ERRO: Especifique a pressão na saída!\n')
else:
    print('condição de bocal supersônico')

eIN = pIN/(k-1) + (1/(2*rhoIN))*(vIN**2+uIN**2)
eIC = pIC/(k-1) + (1/(2*rhoIC))*(vIC**2+uIC**2)


#%% Solver
    
dt = .05
# prealoca vetores de fluxo
F,G = np.zeros((len(celulas),4)),np.ones((len(celulas),4))
# prealoca vetor de estados
Q = np.zeros((len(celulas),4))
Qp = np.zeros((len(celulas),4))
# pressões (pode ajudar...)
P = np.zeros((len(celulas)))
U,V,Rho = np.ones((len(celulas))),np.ones((len(celulas))),np.ones((len(celulas)))

# impõe condição inicial
U = U*0
V = V*0

for i in range(len(celulas)):
    if pontos[celulas[i,1],0] > 50:
       Q[i,0] = 1
       Q[i,3] = 1/(k-1)
       P[i] = 1
       Rho[i] = 1
    else:
       Q[i,0] = 2
       Q[i,3] = 2/(k-1)
       P[i] = 2
       Rho[i] = 2

       


#Q[:,0],Q[:,1],Q[:,2],Q[:,3] = rhoIC,rhoIC*uIC,rhoIC*vIC,eIC
#P[:] = pIC
F[:,0],F[:,1],F[:,2],F[:,3] = Q[:,1], Q[:,1]*0+P[:], Q[:,1]*0, (eIC+pIC)*0
G[:,0],G[:,1],G[:,2],G[:,3] = Q[:,1], Q[:,2]*0, Q[:,2]*0+P[:], (eIC+pIC)*0
#Rho = Rho*rhoIC

# Aqui começa o solver
it = int(input('Entre o numero máximo de iterações:\n'))
nres = int(input('Entre o número de resultados para imprimir:\n'))
print('resolvendo...')
aux = 1;

exportar = [int(it*i/nres) for i in range(nres-1,-1,-1)]

info = {"densidade" : Rho,"pressao" : P,"U":U,"V":V}

# exporta o estado inicial

arquivo = 'resTeste0' 
exportVTK(pontos,celulas,info,arquivo)  
           
# solver supersônico 

if True:
    while it:
        #
        it-=1
        # preditor de estrado com os elementos reais # ARRUMAR
        for el in elementos:
            Qp[el] = Q[el] - (dt/Vol[el])*( \
              F[el,:]*SR[el,1] - F[viz[el,3],:]*SL[el,1] \
              + G[el,:]*SU[el,0] - G[viz[el,2],:]*SD[el,0] \
              - F[el,:]*SU[el,1] + F[viz[el,2],:]*SD[el,1] \
              - G[el,:]*SR[el,0] + G[viz[el,3],:]*SL[el,0] \
              )
        
              
          
        # preditor de fluxos   
   
        for el in elementos:
            Rho[el] = Qp[el,0]
            P[el] = (k-1)*Rho[el]*( Qp[el,3]/Rho[el] - (1/2)*((Qp[el,1]/Qp[el,0])**2+ (Qp[el,2]/Qp[el,0])**2))
            
            F[el,0] = Qp[el,1]
            F[el,1] = (Qp[el,1]**2/Rho[el]) + P[el]
            F[el,3] = (Qp[el,3]+P[el])*Qp[el,1]/Qp[el,0]

            G[el,0] = Qp[el,2]
            G[el,2] = (Qp[el,2]**2/Rho[el]) + P[el]
            G[el,3] = (Qp[el,3]+P[el])*Qp[el,2]/Qp[el,0]

            F[el,2] = F[el,0]*G[el,0]/Rho[el]
            G[el,1] = F[el,2]

        #print('b')

        #preditor de fantasmas
        for fa in Bu:

            #angulo com a parede
            a = np.arctan(SD[fa,1]/SD[fa,0])
            P[fa] = P[viz[fa,2]]                       
            
            F[fa,0] = (a**2*(-F[viz[fa,2],0])+2*a*G[viz[fa,2],0]+F[viz[fa,2],0])/(a**2+1)
            F[fa,1] = (a**2*(-F[viz[fa,2],1])+2*a*G[viz[fa,2],1]+F[viz[fa,2],1])/(a**2+1)+ 2*P[fa]*a**2/(a**2+1)
            F[fa,2] = (a**2*(-F[viz[fa,2],2])+2*a*G[viz[fa,2],2]+F[viz[fa,2],2])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            F[fa,3] = (a**2*(-F[viz[fa,2],3])+2*a*G[viz[fa,2],3]+F[viz[fa,2],3])/(a**2+1)

            G[fa,0] = (a**2*(G[viz[fa,2],0])+2*a*F[viz[fa,2],0]-G[viz[fa,2],0])/(a**2+1)
            G[fa,1] = (a**2*(G[viz[fa,2],1])+2*a*F[viz[fa,2],1]-G[viz[fa,2],1])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            G[fa,2] = (a**2*(G[viz[fa,2],2])+2*a*F[viz[fa,2],2]-G[viz[fa,2],2])/(a**2+1)+ 2*P[fa]/(a**2+1)
            G[fa,3] = (a**2*(G[viz[fa,2],3])+2*a*F[viz[fa,2],3]-G[viz[fa,2],3])/(a**2+1)             
            
#            F[fa,0] = F[viz[fa,2],0]
#            F[fa,1] = F[viz[fa,2],1] 
#            F[fa,2] = F[viz[fa,2],2] 
#            F[fa,3] = F[viz[fa,2],3] 
#
#            G[fa,0] = -G[viz[fa,2],0]
#            G[fa,1] = -G[viz[fa,2],1] 
#            G[fa,2] = -G[viz[fa,2],2] + 2*P[viz[fa,2]]
#            G[fa,3] = -G[viz[fa,2],3] 

        for fa in Bd:
            
            #angulo com a parede
            a = np.arctan(SU[fa,1]/SU[fa,0])
            P[fa] = P[viz[fa,0]]                       
            
            F[fa,0] = (a**2*(-F[viz[fa,0],0])+2*a*G[viz[fa,0],0]+F[viz[fa,0],0])/(a**2+1)
            F[fa,1] = (a**2*(-F[viz[fa,0],1])+2*a*G[viz[fa,0],1]+F[viz[fa,0],1])/(a**2+1)+ 2*P[fa]*a**2/(a**2+1)
            F[fa,2] = (a**2*(-F[viz[fa,0],2])+2*a*G[viz[fa,0],2]+F[viz[fa,0],2])/(a**2+1) - 2*P[fa]*a/(a**2+1)
            F[fa,3] = (a**2*(-F[viz[fa,0],3])+2*a*G[viz[fa,0],3]+F[viz[fa,0],3])/(a**2+1)

            G[fa,0] = (a**2*(G[viz[fa,0],0])+2*a*F[viz[fa,0],0]-G[viz[fa,0],0])/(a**2+1)
            G[fa,1] = (a**2*(G[viz[fa,0],1])+2*a*F[viz[fa,0],1]-G[viz[fa,0],1])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            G[fa,2] = (a**2*(G[viz[fa,0],2])+2*a*F[viz[fa,0],2]-G[viz[fa,0],2])/(a**2+1)+ 2*P[fa]/(a**2+1)
            G[fa,3] = (a**2*(G[viz[fa,0],3])+2*a*F[viz[fa,0],3]-G[viz[fa,0],3])/(a**2+1)            

#            F[fa,0] = F[viz[fa,0],0]
#            F[fa,1] = F[viz[fa,0],1] 
#            F[fa,2] = F[viz[fa,0],2] 
#            F[fa,3] = F[viz[fa,0],3] 
#
#            G[fa,0] = -G[viz[fa,0],0]
#            G[fa,1] = -G[viz[fa,0],1] 
#            G[fa,2] = -G[viz[fa,0],2] + 2*P[viz[fa,0]]
#            G[fa,3] = -G[viz[fa,0],3] 

        for fa in Br:
  
            F[fa,0] = -F[viz[fa,3],0]
            F[fa,1] = -F[viz[fa,3],1] + 2*P[viz[fa,3]]
            F[fa,2] = -F[viz[fa,3],2] 
            F[fa,3] = -F[viz[fa,3],3] 

            G[fa,0] = G[viz[fa,3],0]
            G[fa,1] = G[viz[fa,3],1] 
            G[fa,2] = G[viz[fa,3],2] 
            G[fa,3] = G[viz[fa,3],3]            
            
        for fa in Bl:
            F[fa,0] = -F[viz[fa,1],0]
            F[fa,1] = -F[viz[fa,1],1] + 2*P[viz[fa,1]]
            F[fa,2] = -F[viz[fa,1],2] 
            F[fa,3] = -F[viz[fa,1],3] 

            G[fa,0] = G[viz[fa,1],0]
            G[fa,1] = G[viz[fa,1],1] 
            G[fa,2] = G[viz[fa,1],2] 
            G[fa,3] = G[viz[fa,1],3]     

        #print('b') 
    
        # Corretor de estados 
        for el in elementos:
            Q[el] = .5*(Q[el]+Qp[el]) - (dt/(2*Vol[el]))*( \
                     F[viz[el,1],:]*SR[el,1] - F[el,:]*SL[el,1]  \
                     +G[viz[el,0],:]*SU[el,0] - G[el,:]*SD[el,0] \
                     +F[el,:]*SD[el,1] - F[viz[el,0],:]*SU[el,1] \
                     +G[el,:]*SL[el,0] - G[viz[el,1],:]*SR[el,0] \
                     )

        # corretor de fluxos    
        for el in elementos:
            Rho[el] = Q[el,0]
            P[el] =  (k-1)*Rho[el]*( Q[el,3]/Rho[el] - (1/2)*((Q[el,1]/Qp[el,0])**2+ (Q[el,2]/Q[el,0])**2))
            F[el,0] = Q[el,1]
            F[el,1] = (Q[el,1]**2/Rho[el]) + P[el]
            F[el,3] = (Q[el,3]+P[el])*Q[el,1]/Q[el,0]

            G[el,0] = Q[el,2]
            G[el,2] = (Q[el,2]**2/Rho[el]) + P[el]
            G[el,3] = (Q[el,3]+P[el])*Q[el,2]/Rho[el]

            F[el,2] = F[el,0]*G[el,0]/Rho[el] 
            G[el,1] = F[el,2]
            
            U[el] = F[el,0]/Rho[el]
            V[el] = G[el,0]/Rho[el]
#            if el == 200:
#                print('{}')

        #print('b')
        #corretor de fantasmas
        for fa in Bu:
            #angulo com a parede
            a = np.arctan(SD[fa,1]/SD[fa,0])
            P[fa] = P[viz[fa,2]]                       
            
            F[fa,0] = (a**2*(-F[viz[fa,2],0])+2*a*G[viz[fa,0],0]+F[viz[fa,2],0])/(a**2+1)
            F[fa,1] = (a**2*(-F[viz[fa,2],1])+2*a*G[viz[fa,0],1]+F[viz[fa,2],1])/(a**2+1)+ 2*P[fa]*a**2/(a**2+1)
            F[fa,2] = (a**2*(-F[viz[fa,2],2])+2*a*G[viz[fa,0],2]+F[viz[fa,2],2])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            F[fa,3] = (a**2*(-F[viz[fa,2],3])+2*a*G[viz[fa,0],3]+F[viz[fa,2],3])/(a**2+1)

            G[fa,0] = (a**2*(G[viz[fa,2],0])+2*a*F[viz[fa,0],0]-G[viz[fa,2],0])/(a**2+1)
            G[fa,1] = (a**2*(G[viz[fa,2],1])+2*a*F[viz[fa,0],1]-G[viz[fa,2],1])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            G[fa,2] = (a**2*(G[viz[fa,2],2])+2*a*F[viz[fa,0],2]-G[viz[fa,2],2])/(a**2+1)+ 2*P[fa]/(a**2+1)
            G[fa,3] = (a**2*(G[viz[fa,2],3])+2*a*F[viz[fa,0],3]-G[viz[fa,2],3])/(a**2+1)             

            U[fa] = F[fa,0]/Rho[viz[fa,2]]
            V[fa] = G[fa,0]/Rho[viz[fa,2]]
            
        for fa in Bd:
            
            #angulo com a parede
            a = np.arctan(SU[fa,1]/SU[fa,0])
            P[fa] = P[viz[fa,0]]                       
            
            F[fa,0] = (a**2*(-F[viz[fa,0],0])+2*a*G[viz[fa,0],0]+F[viz[fa,0],0])/(a**2+1)
            F[fa,1] = (a**2*(-F[viz[fa,0],1])+2*a*G[viz[fa,0],1]+F[viz[fa,0],1])/(a**2+1)+ 2*P[fa]*a**2/(a**2+1)
            F[fa,2] = (a**2*(-F[viz[fa,0],2])+2*a*G[viz[fa,0],2]+F[viz[fa,0],2])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            F[fa,3] = (a**2*(-F[viz[fa,0],3])+2*a*G[viz[fa,0],3]+F[viz[fa,0],3])/(a**2+1)

            G[fa,0] = (a**2*(G[viz[fa,0],0])+2*a*F[viz[fa,0],0]-G[viz[fa,0],0])/(a**2+1)
            G[fa,1] = (a**2*(G[viz[fa,0],1])+2*a*F[viz[fa,0],1]-G[viz[fa,0],1])/(a**2+1)- 2*P[fa]*a/(a**2+1)
            G[fa,2] = (a**2*(G[viz[fa,0],2])+2*a*F[viz[fa,0],2]-G[viz[fa,0],2])/(a**2+1)+ 2*P[fa]/(a**2+1)
            G[fa,3] = (a**2*(G[viz[fa,0],3])+2*a*F[viz[fa,0],3]-G[viz[fa,0],3])/(a**2+1)            

            U[fa] = F[fa,0]/Rho[viz[fa,0]]
            V[fa] = G[fa,0]/Rho[viz[fa,0]]

        for fa in Br:
            P[fa] = P[viz[fa,3]]
            
            F[fa,0] = -F[viz[fa,3],0]
            F[fa,1] = -F[viz[fa,3],1] + 2*P[fa]
            F[fa,2] = -F[viz[fa,3],2] 
            F[fa,3] = -F[viz[fa,3],3] 

            G[fa,0] = G[viz[fa,3],0]
            G[fa,1] = G[viz[fa,3],1] 
            G[fa,2] = G[viz[fa,3],2] 
            G[fa,3] = G[viz[fa,3],3]   
            P[fa] = P[viz[fa,3]]

            U[fa] = F[fa,0]/Rho[viz[fa,3]]
            V[fa] = G[fa,0]/Rho[viz[fa,3]]
            
            
            
        for fa in Bl:
            P[fa] = P[viz[fa,1]]
            
            F[fa,0] = -F[viz[fa,1],0]
            F[fa,1] = -F[viz[fa,1],1] + 2*P[fa]
            F[fa,2] = -F[viz[fa,1],2] 
            F[fa,3] = -F[viz[fa,1],3] 

            G[fa,0] = G[viz[fa,1],0]
            G[fa,1] = G[viz[fa,1],1] 
            G[fa,2] = G[viz[fa,1],2] 
            G[fa,3] = G[viz[fa,1],3]    

            U[fa] = F[fa,0]/Rho[viz[fa,1]]
            V[fa] = G[fa,0]/Rho[viz[fa,1]]
        #print('b')
     #   print('-')
        # como exportar
        if exportar[aux] == it:
            arquivo = 'resTeste' + str(aux)  
            aux += 1
            exportVTK(pontos,celulas,info,arquivo)            
           
#%%