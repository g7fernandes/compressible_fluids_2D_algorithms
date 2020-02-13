#!/uSW/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 10:04:47 2018



Este programa lê o arquivo VTK. E armazena os pontos nas linhas da variável 
'pontos' e as células dos volumes na variável 'celulas'. Cada linha em 'celula' 
corresponde a um volume formado pelos pontos indicados nas colunas das linhas 
ligados no sentido horário. Os pontos nas linhas da variável 'célula' corres-
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
import os
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
    
def mulmatriz(M):
    Ans = np.eye(4)
    for m in range(8):
        Ans = np.matmul(Ans,M[m])
    return Ans    
        

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
SW = np.empty((ncells,2))
SE= np.empty((ncells,2))
SN = np.empty((ncells,2))
SS = np.empty((ncells,2))  #superficies em volta do elemento
fantasma = [] #nao necessario mas interessante ver quais sao as celulas no contorno
elementos = []
Bw = []
Be = []
Bn = []
Bs = []

for cact in range(ncells):
    cel = celulas[cact,:]
    # Volumes
    Vol[cact] = PolygonArea([(pontos[cel[0],0],pontos[cel[0],1]), \
                     (pontos[cel[1],0],pontos[cel[1],1]), \
                     (pontos[cel[2],0],pontos[cel[2],1]), \
                     (pontos[cel[3],0],pontos[cel[3],1])])
    # superfícies (normais às)
    SN[cact,:] = np.array([-(pontos[cel[1],1]- pontos[cel[0],1]),(pontos[cel[1],0]- pontos[cel[0],0])])
    SW[cact,:] = np.array([-(pontos[cel[0],1]- pontos[cel[3],1]),(pontos[cel[0],0]- pontos[cel[3],0])])
    SS[cact,:] = np.array([(pontos[cel[2],1]-pontos[cel[3],1]) , -(pontos[cel[2],0]-pontos[cel[3],0])])
    SE[cact,:] = np.array([(pontos[cel[1],1]-pontos[cel[2],1]) , -(pontos[cel[1],0]-pontos[cel[2],0])])
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
        if estacheia(viz[cact,:]):
#            print('ch')
            break #ja encontrou todos os vizinhos
        
        if (edgeY2 == celulas[cviz,2:4]).all():
            viz[cact,0] = cviz
            viz[cviz,2] = cact

        elif (edgeY1 == celulas[cviz,0:2]).all():
            viz[cact,2] = cviz
            viz[cviz,0] = cact
          
        elif (edgeX1 == celulas[cviz,1:3]).all():
            viz[cact,3] = cviz
            viz[cviz,1] = cact
           
        elif (edgeX2 == np.array([celulas[cviz,3],celulas[cviz,0]])).all():
            viz[cact,1] = cviz
            viz[cviz,3] = cact
                       
    
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
        Bn.append(i) #cima
        cor = 'b'
    elif (viz[i,1]) == -1 and (viz[i,0] >= 0 and viz[i,2] >= 0):
        Be.append(i) # direita
        cor = 'r'
    elif (viz[i,2]) == -1 and (viz[i,1] >= 0 and viz[i,3] >= 0):
        Bs.append(i)  # baixo
        cor = 'g'
    elif (viz[i,3]) == -1 and (viz[i,0] >= 0 and viz[i,2] >= 0):
        Bw.append(i)  # esquerda
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

eIN = pIN/(k-1) + (rhoIN/2)*(vIN**2+uIN**2)
eIC = pIC/(k-1) + (rhoIC/2)*(vIC**2+uIC**2)


#%% Solver
    
dt = cfl*dx_min/(8*34300)

# prealoca vetor de estados
Q = np.zeros((len(celulas),4))
# pressões (pode ajudar...)
P = np.zeros((len(celulas)))
U,V,Rho = np.ones((len(celulas))),np.ones((len(celulas))),np.ones((len(celulas)))
U = U*uIC
V = V*vIC
c = np.ones((len(celulas)))*cIC
Mach = np.sqrt(U**2+V**2)/c


# impõe condição inicial
Q[:,0],Q[:,1],Q[:,2],Q[:,3] = rhoIC,rhoIC*uIC,rhoIC*vIC,eIC
P[:] = pIC

Rho = Rho*rhoIC

A,B = np.empty((4,4)), np.empty((4,4))
GammaA,GammaB = np.zeros((4,4)), np.zeros((4,4))


###  Calcula as matrizes   ###

# Matriz de rotação RA *só precisa dos cossenos e senos

beta = k-1

RA,RB = np.empty((len(celulas),2)), np.empty((len(celulas),2))

RAe,RAw,RBn,RBs = np.eye(4),np.eye(4),np.eye(4),np.eye(4)
invRAe,invRAw,invRBn,invRBs = np.eye(4),np.eye(4),np.eye(4),np.eye(4)

RA[:,1] =  SE[:,0]/np.sqrt(SE[:,0]**2 + SE[:,1]**2) #cos aalpha
RA[:,0] =  SE[:,1]/np.sqrt(SE[:,0]**2 + SE[:,1]**2)  #sen aalpha

RB[:,0] =  SN[:,0]/np.sqrt(SN[:,0]**2 + SN[:,1]**2) #sin abeta
RB[:,1] =  SN[:,1]/np.sqrt(SN[:,0]**2 + SN[:,1]**2) #cos abeta

# matrizes S 
S,invS = np.zeros((ncells,4,4)),np.zeros((ncells,4,4))
S[:,0,0],invS[:,0,0] = 1,1
S[:,3,3],invS[:,3,3] = beta, 1/beta

# matrizes C
Ca,invCa = np.zeros((ncells,4,4)),np.zeros((ncells,4,4))
Ca[:,0,0],invCa[:,0,0] = 1,1
Ca[:,2,2],invCa[:,2,2] = 1,1
Ca[:,1,3],Ca[:,3,3] = 1,1
invCa[:,3,:] = np.array([0,.5,0,.5])

Cb, invCb = np.zeros((ncells,4,4)),np.zeros((ncells,4,4))
Cb[:,0,0],Cb[:,1,1],Cb[:,2,3],Cb[:,3,3] = 1,1,1,1
invCb[:,0,0],invCb[:,1,1] = 1,1
invCb[:,3,:] = np.array([0, 0, .5, .5])

GammaA, GammaB = np.zeros((ncells,4,4)),np.zeros((ncells,4,4))
GammaPA, GammaPB = np.zeros((ncells,4,4)),np.zeros((ncells,4,4))
GammaMA, GammaMB = np.zeros((ncells,4,4)),np.zeros((ncells,4,4))

fn,fs,fe,dw = np.empty((3,1)),np.empty((3,1)),np.empty((3,1)),np.empty((3,1))

# preset das condições de contorno nas paredes

for fa in Bn:
    a = -(SS[fa,0]/SS[fa,1])
    P[fa] = P[viz[fa,2]] 

    Q[fa,0] = Q[viz[fa,2],0]
    Q[fa,1] = (a**2*(-Q[viz[fa,2],1])+2*a*Q[viz[fa,2],2]+Q[viz[fa,2],1])/(a**2+1)            
    Q[fa,2] = (a**2*(Q[viz[fa,2],2])+2*a*Q[viz[fa,2],1]-Q[viz[fa,2],2])/(a**2+1)                    
    Q[fa,3] = Q[viz[fa,2],3]

    Rho[fa] = Q[fa,0]            
    U[fa] = Q[fa,1]/Q[fa,0]
    V[fa] = Q[fa,2]/Q[fa,0]                

for fa in Bs:
    a = -(SN[fa,0]/SN[fa,1])
    P[fa] = P[viz[fa,0]] 
    
    Q[fa,0] = Q[viz[fa,0],0]
    Q[fa,1] = (a**2*(-Q[viz[fa,0],1])+2*a*Q[viz[fa,0],2]+Q[viz[fa,0],1])/(a**2+1)            
    Q[fa,2] = (a**2*(Q[viz[fa,0],2])+2*a*Q[viz[fa,0],1]-Q[viz[fa,0],2])/(a**2+1) 
    Q[fa,3] = Q[viz[fa,0],3]
    
    Rho[fa] = Q[fa,0]
    U[fa] = Q[fa,1]/Q[fa,0]
    V[fa] = Q[fa,2]/Q[fa,0]      

# preset das matrizes que variam

for el in range(ncells):
    alpha = (U[el]**2+V[el]**2)/2
    
    Up = RA[el,1]*U[el]+RA[el,0]*V[el]  #u linha
    Vp = RB[el,0]*U[el]+RB[el,1]*V[el] #v linha
    
    S[el,1,0],S[el,1,1] = -Up/Rho[el], 1/Rho[el]
    S[el,2,0],S[el,2,2] = -Vp/Rho[el], S[el,1,1]
    S[el,3,0],S[el,3,1],S[el,3,2] = alpha*beta,-Up*beta,-Vp*beta

    invS[el,1,0],invS[el,1,1] = Up,Rho[el]
    invS[el,2,0],invS[el,2,2] = Vp,Rho[el]
    invS[el,3,0],invS[el,3,1],invS[el,3,2] = alpha, Rho[el]*Up,Rho[el]*Vp                                      

    Ca[el,0,3],Ca[el,1,1] = -1/c[el]**2, Rho[el]*c[el]
    Ca[el,3,1] = -Ca[el,1,1]
    
    invCa[el,0,1], invCa[el,1,1] = 1/(2*c[el]**2), 1/(2*Rho[el]*c[el])
    invCa[el,0,3], invCa[el,1,3] = invCa[el,0,1], -invCa[el,1,1]
    
    Cb[el,0,3], Cb[el,2,2] = -1/c[el]**2, Rho[el]*c[el] 
    Cb[el,3,2] = - Cb[el,2,2]

    invCb[el,0,2], invCb[el,2,2] = 1/(2*c[el]**2), 1/(2*Rho[el]*c[el])
    invCb[el,0,3], invCb[el,2,3] = invCb[el,0,2], -invCb[el,2,2]
    
    GammaA[el,0,0],GammaA[el,1,1]  = Up,Up+c[el] 
    GammaA[el,2,2],GammaA[el,3,3] = Up,Up-c[el]

    GammaB[el,0,0],GammaB[el,1,1]  = Vp,Vp 
    GammaB[el,2,2],GammaB[el,3,3] = Vp+c[el],Vp-c[el]
    
    GammaPA = GammaA*(GammaA>0)
    GammaPB = GammaB*(GammaB>0)
    
    GammaMA = GammaA - GammaPA
    GammaMB = GammaB - GammaPB





##### Aqui começa o solver #####
it = int(input('Entre o numero máximo de iterações:\n'))
nres = int(input('Entre o número de resultados para imprimir:\n'))
arquivo_nome= input('Entre o nome do arquivo de saída:\n')


try:
    os.stat(arquivo_nome)
except:
    os.mkdir(arquivo_nome)       



print('Resolvendo...')
aux = 1;

exportar = [int(it*i/nres) for i in range(nres-1,-1,-1)]

info = {"Mach" : Mach,"pressao" : P,"U":U,"V":V}

# exporta o estado inicial
arquivo_nome = './' + arquivo_nome + '/' + arquivo_nome
arquivo = arquivo_nome + '_0' 
exportVTK(pontos,celulas,info,arquivo)  
           

# solver supersônico 
Qp = np.copy(Q)

#dt = cfl*np.min(Vol/(np.abs(U*SE[:,0]+V*SE[:,1])+np.abs(U*SN[:,0]+V*SN[:,1]) \
#             + c*np.sqrt(np.abs(SE[:,1]*SE[:,0]) + np.abs(SN[:,1]*SN[:,0]))))


if inss:
    while it:
        it-=1
        #

        for el in elementos:
            RAe[1,1],RAe[2,2],RAe[2,1],RAe[1,2] = RA[el,1],RA[el,1],-RA[el,0],RA[el,0]
            
            RAw[1,1],RAw[2,2],RAw[2,1],RAw[1,2] = RA[viz[el,3],1],RA[viz[el,3],1],-RA[viz[el,3],0],+RA[viz[el,3],0]
            
            RBn[1,1],RBn[2,2],RBn[2,1],RBn[1,2] = RB[el,1],RB[el,1],RB[el,0],-RB[el,0]
            
            RBs[1,1],RBs[2,2],RBs[2,1],RBs[1,2] = RB[viz[el,2],1],RB[viz[el,2],1],RB[viz[el,2],0],-RB[viz[el,2],0]
            #######
            invRAe[1,1],invRAe[2,2],invRAe[2,1],invRAe[1,2] = RA[el,1],RA[el,1],RA[el,0],-RA[el,0]
            
            invRAw[1,1],invRAw[2,2],invRAw[2,1],invRAw[1,2] = RA[viz[el,3],1],RA[viz[el,3],1],RA[viz[el,3],0],-RA[viz[el,3],0]
            
            invRBn[1,1],invRBn[2,2],invRBn[2,1],invRBn[1,2] = RB[el,1],RB[el,1],-RB[el,0],RB[el,0]
            
            invRBs[1,1],invRBs[2,2],invRBs[2,1],invRBs[1,2] = RB[viz[el,2],1],RB[viz[el,2],1],-RB[viz[el,2],0],RB[viz[el,2],0]
            

            
            fe = mulmatriz((invS[el],invRAe,invCa[el],GammaPA[el],Ca[el],RAe,S[el],Q[el])) \
                + mulmatriz((invS[viz[el,1]],invRAe,invCa[viz[el,1]],GammaMA[viz[el,1]],Ca[viz[el,1]],RAe,S[viz[el,1]],Q[viz[el,1]]))

            fw = +mulmatriz((invS[el],invRAw,invCa[el],GammaMA[el],Ca[el],RAw,S[el],Q[el])) \
                + mulmatriz((invS[viz[el,3]],invRAw,invCa[viz[el,3]],GammaPA[viz[el,3]],Ca[viz[el,3]],RAw,S[viz[el,3]],Q[viz[el,3]]))                
            
            fn = mulmatriz((invS[el],invRBn,invCb[el],GammaPB[el],Cb[el],RBn,S[el],Q[el])) \
                + mulmatriz((invS[viz[el,0]],invRBn,invCb[viz[el,0]],GammaMB[viz[el,0]],Cb[viz[el,0]],RBn,S[viz[el,0]],Q[viz[el,0]]))

            fs = mulmatriz((invS[el],invRBs,invCb[el],GammaMB[el],Cb[el],RBs,S[el],Q[el])) \
                + mulmatriz((invS[viz[el,2]],invRBs,invCb[viz[el,2]],GammaPB[viz[el,2]],Cb[viz[el,2]],RBs,S[viz[el,2]],Q[viz[el,2]]))        

      
            Qp[el] = Q[el] - (dt/Vol[el])*(fe-fw+fn-fs)
            U[el] = Qp[el,1]/Qp[el,0]
            V[el] = Qp[el,2]/Qp[el,0]
            Rho[el] = Qp[el,0]
            P[el] = (k-1)*Rho[el]*(Qp[el,3]/Rho[el] - (1/2)*(U[el]**2+ V[el]**2))
            
        Q[:] = Qp[:]    
            
        # atualiza as variáveis características dos fantasmas parede
        for fa in Bn:
            a = -(SS[fa,0]/SS[fa,1])
            P[fa] = P[viz[fa,2]] 

            Q[fa,0] = Q[viz[fa,2],0]
            Q[fa,1] = (a**2*(-Q[viz[fa,2],1])+2*a*Q[viz[fa,2],2]+Q[viz[fa,2],1])/(a**2+1)            
            Q[fa,2] = (a**2*(Q[viz[fa,2],2])+2*a*Q[viz[fa,2],1]-Q[viz[fa,2],2])/(a**2+1)                    
            Q[fa,3] = Q[viz[fa,2],3]

            Rho[fa] = Q[fa,0]            
            U[fa] = Q[fa,1]/Q[fa,0]
            V[fa] = Q[fa,2]/Q[fa,0]            
        
        for fa in Bs:
            a = -(SN[fa,0]/SN[fa,1])
            P[fa] = P[viz[fa,0]] 
            
            Q[fa,0] = Q[viz[fa,0],0]
            Q[fa,1] = (a**2*(-Q[viz[fa,0],1])+2*a*Q[viz[fa,0],2]+Q[viz[fa,0],1])/(a**2+1)            
            Q[fa,2] = (a**2*(Q[viz[fa,0],2])+2*a*Q[viz[fa,0],1]-Q[viz[fa,0],2])/(a**2+1) 
            Q[fa,3] = Q[viz[fa,0],3]
            
            Rho[fa] = Q[fa,0]
            U[fa] = Q[fa,1]/Q[fa,0]
            V[fa] = Q[fa,2]/Q[fa,0]   
        
# atualiza valores de saída
        for fa in Be:
            
             P[fa] = P[viz[fa,3]] 
             Q[fa] = Q[viz[fa,3]]
             U[fa] = Q[fa,1]/Q[fa,0]
             V[fa] = Q[fa,2]/Q[fa,0]    
             Rho[fa] = Q[fa,0]      

# Atualização dos valores        
        for el in range(ncells):
            alpha = (U[el]**2+V[el]**2)/2
            
            Up = RA[el,1]*U[el]+RA[el,0]*V[el]  #u linha
            Vp = RB[el,0]*U[el]+RB[el,1]*V[el] #v linha
            
            S[el,1,0],S[el,1,1] = -Up/Rho[el], 1/Rho[el]
            S[el,2,0],S[el,2,2] = -Vp/Rho[el], S[el,1,1]
            S[el,3,0],S[el,3,1],S[el,3,2] = alpha*beta,-Up*beta,-Vp*beta

            invS[el,1,0],invS[el,1,1] = Up,Rho[el]
            invS[el,2,0],invS[el,2,2] = Vp,Rho[el]
            invS[el,3,0],invS[el,3,1],invS[el,3,2] = alpha, Rho[el]*Up,Rho[el]*Vp                                      

            Ca[el,0,3],Ca[el,1,1] = -1/c[el]**2, Rho[el]*c[el]
            Ca[el,3,1] = -Ca[el,1,1]
            
            invCa[el,0,1], invCa[el,1,1] = 1/(2*c[el]**2), 1/(2*Rho[el]*c[el])
            invCa[el,0,3], invCa[el,1,3] = invCa[el,0,1], -invCa[el,1,1]
            
            Cb[el,0,3], Cb[el,2,2] = -1/c[el]**2, Rho[el]*c[el] 
            Cb[el,3,2] = - Cb[el,2,2]

            invCb[el,0,2], invCb[el,2,2] = 1/(2*c[el]**2), 1/(2*Rho[el]*c[el])
            invCb[el,0,3], invCb[el,2,3] = invCb[el,0,2], -invCb[el,2,2]
            
            GammaA[el,0,0],GammaA[el,1,1]  = Up,Up+c[el] 
            GammaA[el,2,2],GammaA[el,3,3] = Up,Up-c[el]

            GammaB[el,0,0],GammaB[el,1,1]  = Vp,Vp 
            GammaB[el,2,2],GammaB[el,3,3] = Vp+c[el],Vp-c[el]
            
            GammaPA = GammaA*(GammaA>0)
            GammaPB = GammaB*(GammaB>0)
            
            GammaMA = GammaA - GammaPA
            GammaMB = GammaB - GammaPB

        c = np.sqrt(k*P/Rho)
#        dt = cfl*np.min(Vol/(np.abs(U*SW[:,0]+V*SW[:,1])+np.abs(U*SN[:,0]+V*SN[:,1]) \
#                     + c*np.sqrt(np.abs(SW[:,1]*SW[:,0]) + np.abs(SN[:,1]*SN[:,0]))))

        Mach[:] = np.sqrt(U**2+V**2)/c           
        # como exportar
        if exportar[aux] == it:
            arquivo = arquivo_nome + '_' + str(aux)  
            aux += 1
            exportVTK(pontos,celulas,info,arquivo)            

#%%