import os
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from pylab import *


barfont = {'fontname':'Times New Roman','weight' : 'bold','size'   : 20}
csfont = {'fontname':'Times New Roman','weight' : 'bold','size'   : 32}
myfont = matplotlib.font_manager.FontProperties(family='times new roman', size=28)
# j : Numberof sharednearest-neighbours
# k: Numberof bonds between sharedneighbours
# l : Numberofbonds in longestbond-chainformed by sharedneighbours

def calculateDistance(a,b):
     dist = math.sqrt((float(a[1])-float(b[1]))**2 +(float(a[2])-float(b[2]))**2+(float(a[3])-float(b[3]))**2)
     return dist	 	

def CommonNeighbor(a,b):
    if (a[0]!=b[0] and calculateDistance(a,b)<3.5):      #3.3 for cubo
       return int(b[0]-1)
    else:
       return 1000000 	
	   
def Generatelist(a,List):
    NN=[]
    for b in List: 
       if CommonNeighbor(a,b)!=1000000:
         NN.append(CommonNeighbor(a,b))
    return(NN)   
	
def Listneighbor(List):
    LN=[]
    for a in range(0,len(List)):
       for b in range(0,len(List)):	
         if (a!=b and CommonNeighbor(List[a],List[b])!=1000000): # and CommonNeighbor(a,b)!=1000000):
           LN.append([a,b])
    return(LN)	
	
	
if not os.path.exists('steps'):
    os.makedirs('steps')
  
f = open('file.lammpstrj', 'r') #Cu_N.lammps
array = []
for line in f.read().split('\n'):
  if line !='':
    array.append(line.split())
f.close()

L=len(array)
print (L)
j=-1
for x in range (0,L):
 if len(array[x])>1:
  if array[x][1]=='TIMESTEP':
    j=j+1 
    print (x+1)
    g = open('steps/coords{0}'.format(j), 'w')
    #g.write('561 \n')
    i=x+9
    print (i)
    while i<x+141:
      g.write('{0} '.format(array[i][0])+'{0} '.format(array[i][1])+'{0} '.format(array[i][2])+'{0} \n'.format(array[i][3]))
      i=i+1
    g.close()
	
	
INDEX=np.genfromtxt('CNAindex').tolist()
E=open('Rod8_class.txt','w')	
for z in range(0,31):
  #print (z)
  data=np.genfromtxt('steps/coords{0}'.format(z)) 
  hcp=[]
  #g = open('CNA{0}'.format(z),'w')
  f = open('type{0}.lammpstrj'.format(z),'w')
  f.write('ITEM: TIMESTEP  \n')
  f.write('{0}\n'.format(z))
  f.write('ITEM: NUMBER OF ATOMS \n')
  f.write('132 \n')
  f.write('ITEM: BOX BOUNDS pp pp pp \n')
  f.write('0.0000000000000000e+00 8.0000000000000000e+01 \n')
  f.write('0.0000000000000000e+00 8.0000000000000000e+01 \n')
  f.write('0.0000000000000000e+00 8.0000000000000000e+01 \n')
  f.write('ITEM: ATOMS id type  x y z  \n')
  TYPE=[]
  for x in range(0, len(data)):
     Allneighbor = Generatelist(data[x],data)
     #print('NN')  #print(x,len(Allneighbor),Allneighbor)
     Allneighbor_list=[data[i] for i in Allneighbor]  
     #print(Allneighbor_list)   # print('CommonNN')
     L=[0,0,1,2,1,5,5]
     LEN=[]
     LEN1=[]
     LEN2=[]
     for y in range(0,len(Allneighbor_list)):
        Nextneighbor = Generatelist(Allneighbor_list[y],Allneighbor_list);LEN.append(len(Nextneighbor));Nextneighbor_list=[data[i] for i in Nextneighbor]  
        bonds=Listneighbor(Nextneighbor_list); preL=len(np.unique(bonds)) #print(Nextneighbor_list)
        LEN1.append(len(bonds)//2);LEN2.append(L[preL]); #print(bonds,np.unique(bonds));
     #g.write(str(x)+'  '+str(len(Allneighbor))+'  '+str(max(LEN2))+'  '+str(min(LEN2))+'\n')
     G=[len(Allneighbor),max(LEN2),min(LEN2)]
     if G in INDEX:
        TYPE.append(INDEX.index(G));m=INDEX.index(G)+1
     else:
        TYPE.append(15); m=16 #print (x,G);
     f.write(str(x)+'  '+str(int(m))+'  '+str(data[x][1])+'  '+str(data[x][2])+'  '+str(data[x][3])+'\n') 
  f.close()	
  print(z,TYPE.count(0),TYPE.count(1),TYPE.count(2),TYPE.count(3),TYPE.count(15))
  E.write(str(z)+'  '+str(TYPE.count(0)) +'  '+str(TYPE.count(1)) +'  '+str(TYPE.count(2))+'  '+str(TYPE.count(3)) +'  '+str(TYPE.count(15))+'\n') 
E.close()  
  
data=np.genfromtxt('Rod8_class.txt')	
df = pd.DataFrame(data)            #,index=data[:,0]) #,columns=['a','b','c','d','e','f','g'])
df['label']=0 
df.loc[df[5]>50,'label'] = 3
df.loc[(df[5]<50)& (df[4]==1 ), 'label'] = 1
df.loc[(df[5]<50)& (df[3]< 2), 'label'] = 2
matrix = df.values
np.savetxt('test.out', matrix, fmt='%i')  
cla=matrix[:,-1].tolist()
print(cla.count(0),cla.count(1), cla.count(2),cla.count(3))



print(matrix)

os.system('mkdir ../Rod_data')
os.system('mkdir ../Dis_data')
os.system('mkdir ../Ico_data')
os.system('mkdir ../TO_data')



for x in range(0,len(matrix)):
  if matrix[x,-2]>95:
    print(x)
    os.system('cp ../Rod2_132.{0}.data ../Dis_data'.format(x))
  if matrix[x,-1]==0:
    print(x)
    os.system('cp ../Rod2_132.{0}.data ../Rod_data'.format(x))
  if matrix[x,-1]==1:
    print(x)
    os.system('cp ../Rod2_132.{0}.data ../Ico_data'.format(x))
  if matrix[x,-1]==2:
    print(x)
    os.system('cp ../Rod2_132.{0}.data ../TO_data'.format(x))
 


#print(data)
#df = pd.DataFrame(data,index=data[:,0],columns=['a','b','c','d','e','f'])
#organize=df[df['f']<400]
#print(organize[(organize['d']>2) & (organize['e']==1)])
#print(organize[(organize['d']>2) & (organize['e']==0)])

   
  #g.close()  
  #print(TYPE)
    #   performance =[TYPE.count(i) for i in range(0,16)]
    #   print(performance)
    #   objects = ('fcc bulk','icosahedral internal twinning plane or HCP','icosahedral spine','icosahedral center','decahedral notch edge','fcc (111) surface','fcc (100) surface','icosahedral surface edge','fcc (111)-(100) edge','fcc (111)-(111) edge','decahedral notch vertex','icosahedral surface vertex or decahedral axial vertex','tetrahedral edge','truncated octahedron vertex','Cuboctahedron corner','Other')
    #   y_pos = np.arange(len(objects))
    #   fig=plt.figure(figsize=(24, 18))
    #   rc('axes', linewidth=3)
    #   plt.axes().set_aspect(15)
    #   plt.barh(y_pos, performance, align='center', alpha=0.5)
    #   plt.yticks(y_pos, objects,**barfont,rotation=10)
    #   plt.xticks(**barfont)
    #   plt.xlabel('count',**csfont)
    #   plt.title('Types of atoms in a imperfect Cuboctahedron559 cluster',**csfont)
    #   plt.savefig('TypeAtom_559Cubo.jpg'.format(z)) 
    #   
    #   #plt.title('Types of atoms in a Ag568 cluster_replica{0}'.format(z),**csfont)
    #   #plt.savefig('TypeAtom{0}.jpg'.format(z)) 
    #   plt.show()
    #   
  
   



