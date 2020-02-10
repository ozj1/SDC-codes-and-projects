import os
from shutil import copyfile
import numpy as np
# print os.getcwd()
   
   
#for x in range (31):
#   os.system('mv TO{0}_127.{0}.lammpstrj TO0_127.{0}.lammpstrj'.format(x))
#   os.system('mv TO{0}_127.{0}.data TO0_127.{0}.data'.format(x))
    #os.system('tail  {0}_window/rst/out.colvars.traj'.format(x))
for x in range (1,11):
    os.system('mkdir Pd_127_{0}'.format(x))
	os.system('cp in.Au script.lammps Pd.eam  Pd_127_{0}'.format(x))
	
   #for y in range(31):
#print(organize[(organize['d']>2) & (organize['e']==1)])
 #   print(x)
    #os.system('tail -n 1  Ico{0}/log.lammps  >> alltail '.format(x))  
      #os.system('cp -r Ico{0}/mini/steps/' .format(x)+'coords{0}'.format(y)+ ' 117Ico/117Ico{0}'.format(31*x-31+y))
#for x in range (1,11):

 #   y=np.genfromtxt('Ico{0}/mini/test.out'.format(x))
 #   print(y[:,-1].tolist())
  
  #os.system('mkdir  Ico{0}/mini'.format(x))
    
    #os.system('sed -i Ico{0}/mini/ \'s/dis1/dis{0}/\' in.min_Ag  CNA.py '.format(x) )


# print ('0.5psCubo.{0}.lammpstrj'.format(x))
#    os.system('tail  {0}_win_step3A/Interactions.lammps'.format(x*0.1))
    #os.system('sed -i \'s/{0}/#{0}/\' 2x2density3050000.out'.format(x*500+3050000) )
   # os.system('sed -i \'s/{0}/#{0}/\' 2x2velocity3050000.out'.format(x*500+3050000) )
    #os.system('sed -i \'s/{0} 999108 1i08000/#{0} 999108 108000/\' 3d2x2x6massdensity3050000.out'.format(x*500+3050000) )
   # os.system('sed -i \'s/{0}/#{0}/\' 4x4velocity3050000.out'.format(x*500+3050000) )
    ##os.system('cat head coords{0}'.format(x*10000) + '> data{0}'.format(x))
   # os#.system('tail -n 570 2psCubo561.{0}.lammpstrj > last/last{0}'.format(x))
    #os.system('cat title.{0} CNA.{0}.lammpstrj >  head.{0}.lammpstrj'.format(x))
    
   # os.system('tail -n 561 min.{0}.lammpstrj > nohead.{0}'.format(x))
#if os.path.isfile('{0}_window/rst/out.colvars.traj'.format(x)):
    #os.system('tail  {0}_window/rst/out.colvars.traj'.format(x))  
