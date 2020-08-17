import numpy as np
import math as m
import os
#from vpython import *
import itertools
import random
from scipy.interpolate import griddata
from operator import itemgetter
from FS_functions_corrected import *
Hartree_to_eV=27.21138602

try: PROJ=raw_input('Do you want a coloured FS (possible if you have atom-proj.xml file)? [y/n]')
except: PROJ=input('Do you want a coloured FS (possible if you have atom-proj.xml file)? [y/n]')
mode='BZ' #BZ or WIGNER
name=raw_input('Name of system (name of .save directory): ')
allk=[]
noneqk=[]
whichk=[]
e=[]

if PROJ=='y':
 try:
  h=open('sum','r')
 except:
  print ('Look at proj.out and see which states you want to sum up. Make file "sum" with ranges of states. Eg. "1 5\n 6 10\n 11 15"\n means: make 3 figures: FS with contributions of states 1-5, next with states 6-10 and next with states 11-15. \n Run program again.')
  exit()
 tmp=h.readlines()
 STATES=[ [ int(m)-1 for m in i.split()[0:2]] for i in tmp if i.split()!=[]]
 print STATES


########READ PARAMETERS
print('Reading parameters from '+name+'.save/data-file-schema.xml file...')
[EF,e,lattice_type,NK_noneq,nk1,nk2,nk3,CELL]=read_parameters(name)

###############BOUNDARIES OF BZ
if lattice_type ==1:
 make_cubic_NN(e,lattice_type)
elif lattice_type ==3:
 make_bcc_NN(e)
elif lattice_type in [2,5]:
 make_trig_NN(e)
elif lattice_type ==6:
 make_tetrag_NN(e)
elif lattice_type==9:
 make_ortorhombic_NN(e)

#########MAKING K-GRID
allk2= make_wigner_cell(e,name)
allk=[ np.array(i[0:3]) for i in allk2]
whichk=[int(i[3]) for i in allk2]
NK=len(allk)
print ('Reciprocal lattice vectors: '),
print e

###############MAKING WIGNER OR BZ
if mode=='WIGNER':
 for i in range(len(allk)):
  allk[i].append(whichk[i])
 allk2=allk


else: 
 print 'Making k-points: brillouin zone...'
 allk2=make_BZ(e,allk,whichk)
 print(' BZ- Done!')

 print(' No of points in BZ:'),
 print len(allk2)

nkp=len(allk2)


########READ BANDS
print('Reading eigenvalues...')
if PROJ=='y':
 print('Reading eigenvalues and colours from '+name+'.save/atomic_proj.xml file...')
 [BAND1, COLOUR1]=read_bands(name,STATES,EF)
else: 
 print('Reading eigenvalues from '+name+'.save/data_file_schema.xml file...')
 BAND1=read_notproj_bands(name,EF)
nb=len(BAND1)

print('Print Fermi surface to file...')
h=open('FS.dat','w')
h.write(str(nkp)+' '+str(nb)+' '+str(EF)+' ')
for b in range(nb):
 h.write('0 ')
h.write('\n')
for i in range(len(allk2)):
 h.write(str(allk2[i][0])+' '+str(allk2[i][1])+' '+str(allk2[i][2])+' ')
 for b in range(nb):
  h.write(str(BAND1[b][allk2[i][3]])+' ')
 h.write('\n')
h.close()

if PROJ=='y':
 for i in range(len(COLOUR1)):
  h=open('FS_colour'+str(i)+'.dat','w')
  for j in range(len(allk2)):  
 # for i in range(ns):
   for b in range(nb):
    h.write(str((COLOUR1[i][b][allk2[j][3]]))+' ')
  #  h.write(str(int(COLOUR1[i][b][allk2[j][3]]*256))+' ')
   h.write('\n')
  h.close()
 print len(COLOUR1[0])

n=name.replace('/','-')
print('Copying results to dir '+name+'.save/../'+n+'-results')
os.system('mkdir '+name+'.save/../'+n+'-results; cp sum noneqs* FS.dat BZ.dat FS_colo* FS_NN.dat WRZ.dat  '+name+'.save/../'+n+'-results')



