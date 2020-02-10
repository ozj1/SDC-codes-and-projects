import os
import matplotlib.pyplot as plt
import numpy as np
import math
from pylab import *
import pandas as pd

barfont = {'fontname':'Times New Roman','weight' : 'bold','size'   : 20}
csfont = {'fontname':'Times New Roman','weight' : 'bold','size'   : 32}
myfont = matplotlib.font_manager.FontProperties(family='times new roman', size=28)

fig=plt.figure(figsize=(24, 18))
#Log=np.genfromtxt('log.lammps')
#Time=Log[:,0].tolist()
allPE=[]
allT=[]
for x in range(0,31): 
  print(x)
  #Replica=Log[:,x+1].tolist()
  Thermo=np.genfromtxt('thermo.{0}.lammps'.format(x))[:26000]
  #PE_ave=np.mean(Thermo[200000:206000,3])
  #T_ave=np.mean(Thermo[200000:206000,1])
  #print (PE_ave);allPE.append(PE_ave);allT.append(T_ave)
  plt.plot(Thermo[:,0],Thermo[:,3],alpha=0.3,label='real')
  #Temp=np.arange(300,910,25)
  
  #[300.0,310.0,320.0,330.0,340.0,350.0,360.0,370.0,380.0,390.0,400.0,425.0,450.0,475.0,500.0,525.0,550.0,575.0,600.0,625.0,650.0,675.0,700.0,725.0,750.0,775.0,800.0]
  #Alltemp=[Temp[int(y)] for y in Replica] 
  #print(Replica);print(len(Temp));print(Alltemp)
  
#plt.scatter(allT,allPE,label='set')

plt.yticks(**barfont)
plt.xticks(**barfont)
plt.ylabel('PE (ev)',**csfont)
plt.xlabel('Timestep',**csfont)  #Temperature
#plt.legend(fontsize=20)
#plt.title('Temperature (K) '.format(x),**csfont)
plt.savefig('PEvstime.jpg', dpi=300) 
plt.show()
#plt.savefig('graph/Temp{0}.jpg'.format(x)) 
  #plt.clf()
#Replica_total=[]
#E_all=[]
#Log=np.genfromtxt('log.lammps')



#  for x in range(0,27): 
#    Replica=Log[4:,x+1].tolist()
#    Replica_total=Replica_total+Replica
#    Thermo=np.genfromtxt('thermo.{0}.lammps'.format(x))[:10000]
#    Thermo=Thermo[15::4][:,3].tolist()
#    #print(Thermo) 
#    E_all=E_all+Thermo
#    #print(len(E_all[:20]),len(Replica_total))

Combine=np.append([Replica_total],[E_all], axis=0)

Combine=np.transpose(Combine)
print(Combine[:20,:])
#print(Combine[np.where(Combine[:20,0]==0)])

#np.savetxt('test.out', surface, fmt= '%i ' '%1.3f')
# 
fig=plt.figure(figsize=(24, 18))
for x in range(0,27):
  print(x)
  Einbin=Combine[np.where(Combine[:,0]==x)][:,1]
  print(Einbin)
  
  N, bins, patches =plt.hist(Einbin,bins=100,alpha=0)
  plt.plot(bins[:-1],N)
#plt.savefig('Energy.jpg', dpi=300) 
#plt.show()  

 
 
 
 
  