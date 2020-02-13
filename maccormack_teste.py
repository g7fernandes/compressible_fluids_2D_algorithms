import numpy as np
import matplotlib.pyplot as plt
import time

## Setup
c =  343
lx = 10
CFL = 0.9
T = 300 # temperatura
cv = 0.718
k_1 = .4 # k-1 onde k e cp/cv
n = 1001 # No de elementos
steps = 5000

# saida 
P = np.zeros((steps,n+2))
Rho = np.zeros((steps,n+2))
U = np.zeros((steps,n+2))
## Volumes
# (n+2) adicionara 2 volumes fantasma
elem = np.array([[0]*2]*(n+2) ,dtype=np.float64)
Q = np.array([[0]*3]*(n+2) ,dtype=np.float64)
F = np.array([[0]*3]*(n+2) ,dtype=np.float64)

# indice linha = Nº do volume


# se vizinho = -1, entao e' contorno (fantasma)
elem[0,0] = -1 
elem[n-1,1] = -1

for i in range(1,n-1):
    elem[i,0] = i-1; elem[i,1] = i+1;

dx = lx/n;
vol = dx*1; # volumes sempre iguais

## C.I. 

# p = rho = 1 para x < 1, p = rho = 2 para x >= 1, u = 0

fu = lambda x: 300/(1+25*x**2)

rho = 1;
p = 1;

for i in range(n+2):
    
    u = fu(dx*i-5) 
#    rho = 1+u/100
#    p = 1+ u/400
    
    e = p/(k_1)
        
    Q[i,0] = rho
    Q[i,1] = rho*u
    Q[i,2] = e
    F[i,0] = Q[i,1]
    F[i,1] = rho*u**2+p
    F[i,2] = (e+p)*u
    P[0,i] = p;
    Rho[0,i] = rho;

    
#fantasma
F[0,0] = -F[1,0]
F[0,1] = -F[1,1]+2*1 #2x a pressao 1
F[0,2] = -F[1,2]    
F[n+1,0] = -F[n,0]
F[n+1,1] = -F[n,1]+2*2
F[n+1,2] = -F[n,2]

    
Qp = Q #passo preditor
    

# MacCormack
    
dt = .9*dx/(2*c) # alterar?
a1 = dt/vol
a2 = dt/(2*vol)
S = 1 # superficie unitaria

for step in range(1,steps):
    # i denota o elemento
    
    #preditor
    Qp[0,:] = Q[0,:] - a1*(F[0,:]*S-F[n+1,:]*S) 
    #corretor 
    Q[0,:] = 0.5*(Q[0,:]+Qp[0,:])-a2*(F[1,:]*S-F[0,:]*S) 
    
    #preditor
    Qp[n+1,:] = Q[n+1,:] - a1*(F[n+1,:]*S-F[n,:]*S) 
    #corretor 
    Q[n+1,:] = 0.5*(Q[n+1,:]+Qp[n+1,:])-a2*(F[0,:]*S-F[n+1,:]*S) 

    for i in range(1,n+1):
        #preditor
        Qp[i,:] = Q[i,:] - a1*(F[i,:]*S-F[(i-1),:]*S) 
        #corretor 
        Q[i,:] = 0.5*(Q[i,:]+Qp[i,:])-a2*(F[i+1,:]*S-F[i,:]*S) 
#        Q[i,:] = Q[i,:]-a1*(F[i+1,:]*S-F[i-1,:]*S)     
    # atualiza valores de fluxo
    for i in range(0,n+2):
 #       print('{} {} {} \n'.format(Q[i,0],Q[i,1],Q[i,2]))
        Rho[step,i] = Q[i,0]
#        print('{} {} {}'.format( k_1*Rho[step,i] , Q[i,2]/Rho[step,i] , (1/2)*(Q[i,1]/Q[i,0] )**2))
        P[step,i] = k_1*Rho[step,i]*( Q[i,2]/Rho[step,i] - (1/2)*(Q[i,1]/Q[i,0])**2 )
        F[i,0] = Q[i,1]
        F[i,1] = (Q[i,1]**2/Q[i,0]) + P[step,i]
        F[i,2] = (Q[i,2]+P[step,i])*Q[i,1]/Q[i,0]
        U[step,i] = Q[i,1]/Q[i,0]
        
#    dt = CFL*dx/(c+np.amax(np.abs(U[step,:])))
#    a1 = dt/vol
#    a2 = dt/(2*vol)
        

  #  print('{} {} {}'.format(np.amax(np.abs(Q)),np.amax(U),dt))
    if np.amax(Q) > 10**6:
        print('{}'.format(step))
        break
    
#    if not(step % 100):
#        plt.plot(Q[:,0])
#        plt.plot(Q[:,1])        
#        plt.plot(Q[:,2])
#        plt.show()
         
plt.plot(U[9,1:n+1],'k')
plt.plot(U[90,1:n+1],'b')
plt.plot(U[500,1:n+1],'g')
plt.plot(U[900,1:n+1],'r')
plt.show()    
        
plt.plot(P[9,1:n+1],'k')
plt.plot(P[90,1:n+1],'b')
plt.plot(P[500,1:n+1],'g')
plt.plot(P[900,1:n+1],'r')
plt.show()    

    
