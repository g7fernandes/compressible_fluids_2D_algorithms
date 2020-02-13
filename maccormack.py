import numpy as np
import matplotlib.pyplot as plt

## Setup
lx = 2
CFL = 0.9
T = 300 # temperatura
cv = 0.718
k_1 = .4 # k-1 onde k e cp/cv
n = 100 # No de elementos
steps = 40

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

u = 0;

for i in range(n+2):
    if i*dx > 1:
        rho = 1; p = 1;
    else:
        rho = 2; p = 2;
    
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

    
Qp = np.copy(Q) #passo preditor
    

# MacCormack
print('{} {}'.format(p,rho))
c = np.sqrt(1.4*p/rho)
    
dt = CFL*dx/(c) # alterar?
a1 = dt/vol
a2 = dt/(2*vol)
S = 1 # superficie unitaria

dtmin = 10
dtmax = 0

for step in range(1,steps):
#    print(dt)
    # i denota o elemento
    for i in range(1,n+1):         #preditor
        Qp[i,:] = Q[i,:] - a1*(F[i,:]*S-F[(i-1),:]*S) 
        
    # atualiza valores de fluxo preditor
    for i in range(1,n+1):
        Rho[step,i] = Qp[i,0]
        P[step,i] = k_1*Rho[step,i]*( Qp[i,2]/Rho[step,i] - (1/2)*(Qp[i,1]/Qp[i,0])**2 )
   
        F[i,0] = Qp[i,1]
        F[i,1] = (Qp[i,1]**2/Q[i,0]) + P[step,i]
        F[i,2] = (Qp[i,2]+P[step,i])*Qp[i,1]/Q[i,0]

    # atualiza fantasmas     
    F[0,0] = -F[1,0]
    F[0,1] = -F[1,1]+2*P[step,1] #2x a pressao 1
    F[0,2] = -F[1,2]    
    F[n+1,0] = -F[n,0]
    F[n+1,1] = -F[n,1]+2*P[step,n]
    F[n+1,2] = -F[n,2]

    for i in range(1,n+1):      #corretor   
        Q[i,:] = 0.5*(Q[i,:]+Qp[i,:])-a2*(F[i+1,:]*S-F[i,:]*S) 
        
    # atualiza valores de fluxo coretor
    for i in range(1,n+1):
 #       print('{} {} {} \n'.format(Q[i,0],Q[i,1],Q[i,2]))
        Rho[step,i] = Q[i,0]
#        print('{} {} {}'.format( k_1*Rho[step,i] , Q[i,2]/Rho[step,i] , (1/2)*(Q[i,1]/Q[i,0] )**2))
        U[step,i] = Q[i,1]/Q[i,0]
        P[step,i] = k_1*Rho[step,i]*( Q[i,2]/Rho[step,i] - (1/2)*(U[step,i])**2 )
        F[i,0] = Q[i,1]
        F[i,1] = (U[step,i]**2*Q[i,0]) + P[step,i]
        F[i,2] = (Q[i,2]+P[step,i])*U[step,i]
        
        
    dt = CFL*dx/(c+np.amax(np.abs(U[step,:])))

    if dt < dtmin:
        dtmin = dt
    if dt > dtmax:
        dtmax = dt    
    a1 = dt/vol
    a2 = dt/(2*vol)
        
    # atualiza fantasmas     
    F[0,0] = -F[1,0]
    F[0,1] = -F[1,1]+2*P[step,1] #2x a pressao 1
    F[0,2] = -F[1,2]    
    F[n+1,0] = -F[n,0]
    F[n+1,1] = -F[n,1]+2*P[step,n]
    F[n+1,2] = -F[n,2]
  #  print('{} {} {}'.format(np.amax(np.abs(Q)),np.amax(U),dt))
    if np.amax(Q) > 10**6:
        print('{}'.format(step))
        break
    
#    if not(step % 100):
#        plt.plot(Q[:,0])
#        plt.plot(Q[:,1])        
#        plt.plot(Q[:,2])
#        plt.show()
#        time.sleep(1)
        
    
LX = np.linspace(0,lx,n)    
        
plt.plot(LX,U[0,1:n+1],'k')
plt.plot(LX,U[int(steps/3),1:n+1],'b')
plt.plot(LX,U[int(steps*2/3),1:n+1],'g')
plt.plot(LX,U[steps-1,1:n+1],'r')
plt.show()    
        
plt.plot(LX,P[0,1:n+1],'k')
plt.plot(LX,P[int(steps/3),1:n+1],'b')
plt.plot(LX,P[int(2*steps/3),1:n+1],'g')
plt.plot(LX,P[steps-1,1:n+1],'r')
plt.show()    

    
