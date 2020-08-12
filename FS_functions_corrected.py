import numpy as np
import os
from operator import itemgetter
import  xml.etree.ElementTree as ET
PRECIS=4
Hartree_to_eV=27.21138602
##########################################################################################
def sorting(allk2):
 Xall=[]
 allk2=sorted(allk2, key=itemgetter(0))
 i=0
 while i<len(allk2): 
  X=[]
  x=allk2[i]
  while i<len(allk2) and x[0]==allk2[i][0]:
   X.append(allk2[i])
   i=i+1
  if len(X)>1: X=sorted(X, key=itemgetter(1))
  Xall.append(X)

 Yall=[]
 for X in Xall:
  x=X[0]
  i=0
  while i<len(X): 
   Y=[]
   x=X[i]
   while i<len(X) and x[1]==X[i][1]:
    Y.append(X[i])
    i=i+1
   Y=sorted(Y, key=itemgetter(2))
   Yall.append(Y)

 allk=[]
 for i in Yall:
  for j in i:
   allk.append(j)

 print ' Sorting - Done!'

 return allk

##########################################################################################
def make_wigner_cell(e,name):
 pm=[-1,0.,1]
 A_vectors=[(h*e[0]+k*e[1]+l*e[2]) for h in pm for k in pm for l in pm if (not (h==0 and k==0 and  l==0))] 
 e=np.array([[round(m,PRECIS) for m in k] for k in e])
 einv=np.linalg.inv(e)


 #symmetry operations
 symm=[]
 tree = ET.parse(name+'.save/data-file-schema.xml')
 root = tree.getroot()
 for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     symm.append(np.dot((einv),(np.dot(tmp2,e))))

 symm.append(np.array([[1,0,0],[0,1,0],[0,0,1]]))


 #kpoints
 nb=(root.find('output/band_structure/nks').text)
 noneqs=[]
 for neighbor in root.iter('k_point'):
     tmp=neighbor.text.split()
     tmp2=np.array( [round(float(m),PRECIS) for m in tmp[0:3]])
     noneqs.append(tmp2)



 print(' Number of nonequivalent k-points:'+str(len(noneqs))+' should equal '+str(nb))
 h=open('noneqs0.dat','w')
 h.write('#reciprocal lattice vectors:\n# ')
 for n in range(3): 
  for m in range(3): h.write(str(e[n][m])+' ')
  h.write('\n# ')
 h.write(str(len(noneqs))+' kpoints:\n')
 for x in noneqs:
   for m in range(3): h.write(str(x[m])+' ')
   h.write('\n')
 h.close()


 #equivalent kpoints
 maxe=[ round(max([ e[i][m] for i in range(3)]),PRECIS) for m in range(3)]

 allk=[]
 mmm=0
 for i in noneqs:
  k_tmp=[i]
  for j in symm:
   x=np.dot(j,i) 
   #if x[0]<maxe[0] and x[1]<maxe[1] and x[2]<maxe[2] and x[0]>-maxe[0] and x[1]>-maxe[1] and x[2]>-maxe[2]:
   allk.append(np.array([ round(x[0],PRECIS),round(x[1],PRECIS) ,round(x[2],PRECIS),int(mmm)]))
  mmm=mmm+1


 allk2=sorting(allk)

 allk=[] 
 for i in range(len(allk2)-1):
  x=allk2[i]
  y=allk2[i+1]
  if not(x[0]==y[0] and x[1]==y[1] and x[2]==y[2]):
   allk.append(x)

 h=open('noneqs.dat','w')
 h.write('#reciprocal lattice vectors:\n# ')
 for n in range(3): 
  for m in range(3): h.write(str(e[n][m])+' ')
  h.write('\n# ')
 h.write(str(len(allk))+' kpoints:\n')
 for x in allk:
   for m in range(4): h.write(str(x[m])+' ')
   h.write('\n')
 h.close()
 
 return allk


##########################################################################################
def make_BZ(e,allk,whichk):
 e=np.array(e)
 pm=[-1,0.,1]
 A_vectors=[0.5*(h*e[0]+k*e[1]+l*e[2]) for h in pm for k in pm for l in pm if (not (h==0 and k==0 and  l==0))] 
 print ' A_vectors:',
 for i in A_vectors:
  print i,
 k_list=[]

 mmm=0
 for k_point in allk:  
  for h1 in pm:
   for k1 in pm:
    for l1 in pm:
     k_point2=k_point-(h1*e[0]+k1*e[1]+l1*e[2]) 
     sign0=0
     for a in A_vectors:
#      if k_point2[0]*a[0]>=0 and k_point2[1]*a[1]>=0 and k_point2[2]*a[2]>=0:
        ra=k_point2-a
        scalar_product=sum([a[i]*ra[i] for i in range(3)])
        if scalar_product>0: 
         sign0=1
         break
     if sign0==0:
      k_list.append([round(kk,PRECIS) for kk in k_point2])
      k_list[-1].append(whichk[mmm])
  mmm=mmm+1


 k_list=sorting(k_list)
 nkp=len(k_list)


###WRITE TO FILE
 h=open('BZ.dat','w') 
 for i in k_list:
  for m in range(3): h.write(str(i[m])+' ')
  h.write('\n')
 h.close()


 return k_list

#####################################################################################################
#MAKE_NN:
def round_table(a):
 return [round(m,PRECIS) for m in a]

def check_el(a,b):
 sign=0
 for i in range(len(a)):
  for j in range(len(a[i])):
   if a[i][j]!=b[i][j]:
    sign=1
    break
 if sign==1 and  a[0][0]==b[1][0] and a[0][1]==b[1][1] and a[0][2]==b[1][2] and a[1][0]==b[0][0] and a[1][1]==b[0][1] and a[1][2]==b[0][2]:
    sign=0
 if sign==0: return 1 #if are equal
 else: return 0

def check(NN):
 sign=0
 MM=NN[:-1]
 a=NN[-1]
 if a[0][0]==a[1][0] and a[0][1]==a[1][1] and a[0][2]==a[1][2]: return 1 #if point-thesamepoint connection
 for i in MM:
  if check_el(i,a):
   sign=1
   break
 if sign==1: return 1 #if is in list
 else: return 0


def make_ortorhombic_NN(e):
 WRZ=[]
 NN=[]
 a=0.5*(e[0]+e[1])
 b=0.5*e[1]
 y=(a[1]*a[1]+a[0]*a[0]-a[0]*b[0]-a[0]*b[1]*b[1]/b[0])/(a[1]-a[0]*b[1]/b[0])
 x=b[0]-b[1]*(y-b[1])/b[0]
 a0=np.array([x,y,-e[2][2]/2.])
 a1=np.array([-e[1][1]*2,0,-e[2][2]/2.])
 pm0=np.array([[1,1,1],[-1,1,1],[-1,-1,1],[1,-1,1],[1,1,-1],[-1,1,-1],[-1,-1,-1],[1,-1,-1],[-1,-1,-1]])
 pm=np.array([[1,1,1],[-1,1,1],[1,1,-1]])
 WRZ=[]
 for j in pm0:
  beg0=len(WRZ)
  for i in pm:
   a=a0*j
   WRZ.append(i*a)
   if len(WRZ)>1:   
    NN.append([round_table(WRZ[beg0]),round_table(WRZ[-1])])
   if len(NN)>2 and check(NN): NN=NN[:-1]
  a=a1*j
  beg=len(WRZ)
  for i in pm:
   WRZ.append(i*a)
   NN.append([round_table(WRZ[beg]),round_table(WRZ[-1])])
   if len(NN)>2 and check(NN): NN=NN[:-1]
  NN.append([round_table(WRZ[beg0]),round_table(WRZ[beg])])
  if check(NN): NN=NN[:-1]

 h=open('WRZ.dat','w')
 for i in WRZ:
   for m in range(3):
    h.write(str(i[m])+' ')
   h.write('\n')
 h.close()
 h=open('FS_NN.dat','w')
 for i in NN:
  for m in range(3):
   h.write( str(i[0][m])+' '+str(i[1][m])+' ' )
  h.write('\n')
 h.close()

def make_cubic_NN(e,lattice_type):

 WRZ=[]
 NN=[]
 G=np.linalg.norm(e[0]+e[1]+e[2])+1e-2
 pm=[0,1]
 A_vectors1=[(h*e[0]+k*e[1]+l*e[2]) for h in pm for k in pm for l in pm if ( not (h==0 and k==0 and  l==0))]  
 A_vectors=[]
 for i in A_vectors1:
  A_vectors.append(i)
  A_vectors.append(-i)


 for ii in range(len(A_vectors)):
  WRZi=[]
  for jj in range(len(A_vectors)):
   for kk in range(jj+1,len(A_vectors)):
    n1=A_vectors[ii]
    n2=A_vectors[jj]
    n3=A_vectors[kk]
    if abs(np.dot(n1,n2))<1e-3 or abs(np.dot(n1,n3))<1e-3 or abs(np.dot(n3,n2))<1e-3: continue
   #SOLVE THE EQUATION [n1; n2; n3]*[x,y,z]=[0.5n1**2; 0.5n2**2,0.5n3**2]
    W=np.array([n1,n2,n3])
    A=np.array([0.5*np.linalg.norm(n1)**2, 0.5*np.linalg.norm(n2)**2, 0.5*np.linalg.norm(n3)**2])
    dW=np.linalg.det(W)
    if dW==0: continue
    detW123=[]
    for i in range(3):
     W1=np.transpose(np.array([n1,n2,n3]))
     W1[i]=A
     W1=np.transpose(W1)
     detW123.append(np.linalg.det(W1)/dW)
    a=np.array(detW123)
    sign0=0
    #check if is in BZ 
    for av in A_vectors:
      ra=a-0.5*av
      scalar_product=sum([0.5*av[i]*ra[i] for i in range(3)])
      if scalar_product>1e-2: 
       sign0=1
       break
    if sign0==0: 
     WRZ.append(a)
     WRZi.append(a)
  if lattice_type==5:
   if ii in [4,5,8,9,10,11]: n_NN=2
   else: n_NN=1
  else: n_NN=2
  for i in WRZi:
   odl=[[np.linalg.norm(i-j),j] for j in WRZi if  np.linalg.norm(i-j)>1e-2]
   odl=sorted(odl, key=itemgetter(0))

   for j in odl[:n_NN]:
    NN.append([i,j[1]])

 h=open('WRZ.dat','w')
 for i in WRZ:
   for m in range(3):
    h.write(str(i[m])+' ')
   h.write('\n')
 h.close()
 h=open('FS_NN.dat','w')
 for i in NN:
  for m in range(3):
   h.write( str(i[0][m])+' '+str(i[1][m])+' ' )
  h.write('\n')
 h.close()
########


def make_tetrag_NN(e):
 e2=[ i/2. for i in e]
 WRZ=[ np.array([0,0,0]),e[0],e[1],e[2],e[0]+e[1],e[0]+e[2],e[1]+e[2],e[0]+e[1]+e[2]]
 WRZ=[ i-e2[0]-e2[1]-e2[2] for i in WRZ]
 h=open('FS_NN.dat','w')
 for i in [0,1,4,2,0,3,6,2,4,7,6,3,5,7,4,1,5]:
  for m in range(3):
   h.write( str(WRZ[i][m])+' ' )
  h.write('\n')
 h.close()

#################
def calc_point_of_three_planes(a,b,c):
 A=[a,b,c]
 D=[np.dot(a,a),np.dot(b,b),np.dot(c,c)]
 return np.dot(np.linalg.inv(A),D)

###############
def make_trig_NN(e):
 einv=np.linalg.inv(e)
 a=0.5*e[0]+0.5*e[2]
 b=0.5*e[0]+0.5*e[1]+0.5*e[2]
 c=0.5*e[0]
 [a,b,c]=calc_point_of_three_planes(a,b,c)
 [a,b,c]=a*einv[0]+b*einv[1]+c*einv[2]

 NN=[[a,b,c],[c,b,a],[b,-b,c],[c,-b,b],\
[a,b,c],[a,c,b],[c,b,-b],[b,-b,-c],[b,-c,-b],[c,-b,b],\
[c,-b,b],[b,-b,c],[-b,-c,b],[-c,-a,-b,],[-b,-a,-c],[b,-c,-b],\
[b,-b,-c],[-b,-c,-a],[-b,-a,-c],\
[-b,-c,-a],[-c,-b,-a],[-b,b,-c],[-b,b,-c],[b,c,-b],[c,b,-b],\
[a,c,b],[c,a,b],[b,c,-b],\
[c,a,b],[b,a,c],[-b,c,b],[-c,b,-b],[-b,b,-c],\
[-c,-b,-a],[-a,-b,-c],[-c,b,-b],\
[-a,-b,-c],[-a,-c,-b],[-c,-b,b],[-b,b,c],[-b,c,b],[b,a,c],[b,c,a],\
[c,b,a],[b,c,a],[-b,b,c],[-c,-b,b],[-b,-c,b],[-c,-a,-b],[-a,-c,-b]]
 h=open('FS_NN.dat','w') 
 for j in NN:
   jp=np.dot(j,e)
   for k in jp:
    h.write(str(k)+' ')
   h.write('\n')
 h.close()

 '''
 WRZ=[[a,b,c],[c,b,a],[b,-b,c],[c,-b,b],\
 [a,c,b],[c,b,-b],[b,-b,-c],[b,-c,-b],\
 [-b,-c,b],[-c,-a,-b],[-b,-a,-c],\
 [-b,-c,-a],[-c,-b,-a],[-b,b,-c],[b,c,-b],\
 [c,a,b],[-c,b,-b],[-a,-b,-c],[-a,-c,-b],\
 [-c,-b,b],[-b,b,c],[-b,c,b],[b,a,c],[b,c,a]]
 h=open('FS_WRZ.dat','w') 
 for j in WRZ:
   jp=np.dot(j,e)
   for k in jp:
    h.write(str(k)+' ')
   h.write('\n')
 h.close()
 '''
##########################################################################################
def read_parameters(name):
 e=[]
 os.system('rm aaa')
 os.system('grep  "<nelec>" '+name+'.save/data-file-schema.xml >aaa')
 os.system('grep "<fermi_energy>" '+name+'.save/data-file-schema.xml  >>aaa')
 os.system('grep -A 1 " <starting_k_points>" '+name+'.save/data-file-schema.xml |tail -n 1 >>aaa')
 os.system('grep  "<b1>" '+name+'.save/data-file-schema.xml >>aaa')
 os.system('grep  "<b2>" '+name+'.save/data-file-schema.xml >>aaa')
 os.system('grep  "<b3>" '+name+'.save/data-file-schema.xml >>aaa')
 os.system('grep -m 1  "<nks>" '+name+'.save/data-file-schema.xml >>aaa')
 os.system('grep -m 1  "bravais_index=" '+name+'.save/data-file-schema.xml >>aaa')
 os.system('grep -m 1  "<atomic_structure" '+name+'.save/data-file-schema.xml >>aaa')


 h=open('aaa','r')
 NEL=float(h.readline().split()[0][7:-8])
 EF=float(h.readline().split()[0][14:-15])*Hartree_to_eV
 tmp=h.readline().split()
 [nk1,nk2,nk3]=[int(tmp[1][5:-1]),int(tmp[2][5:-1]),int(tmp[3][5:-1])]

 for i in range(3):
  tmp=h.readline().replace('>',' ').replace('<',' ').split()
  e.append(np.array([float(i) for i in tmp[1:4]]))

 NK_noneqk=int(h.readline().replace(">",' ').replace("<",' ').split()[1])

 lattice_type=int(h.readline().replace('"', ' ').split()[-2]) 

 NAT=int(h.readline().replace('"', ' ').split()[2])
 h.readline()

 AT_POS=[]
 os.system('grep  -m '+str(NAT)+'  "<atom " '+name+'.save/data-file-schema.xml >>aaa')
 for i in range(NAT):
  AT_POS.append(np.array(h.readline().replace('>',' ').replace('<',' ').split()[3:6]))

 os.system('grep -m 1 -A 4 "<cell>" '+name+'.save/data-file-schema.xml >>aaa')
 h.readline()
 CELL=[]
 for i in range(4):
  CELL.append([float(m) for m in h.readline().replace('>',' ').replace('<',' ').split()[1:4] ])

 print (' Number of el.='+str(NEL))
 print (' Fermi energy='+str(EF))
 print (' No of kpoints= '+str(nk1)+' '+str(nk2)+' '+str(nk3))
 print (' Reciprocal lattice vectors'),
 for i in e:
  print(i),
 print ('\n Lattice type: '+str(lattice_type))
 print (' No of nonequivalent atoms: '+str(NAT))
 print (' Positions of atoms: '),
 for i in AT_POS:
  print(i),
 print('\n Cell vectors: '),
 for i in CELL:
  print(i),
 print('\n')

 return [EF,e,lattice_type,NK_noneqk,nk1,nk2,nk3,CELL]

##########################################################################################
def read_notproj_bands(name,EF):
 BAND=[]
 BAND1=[]
 tree = ET.parse(name+'.save/data-file-schema.xml')
 root = tree.getroot()
 for neighbor in root.iter('eigenvalues'):
     BAND.append([ float(m)*Hartree_to_eV for m in neighbor.text.split()])
 BAND=np.transpose(np.array(BAND))
 for b in range(len(BAND)):
  if min(BAND[b])<EF and max(BAND[b])>EF:
    BAND1.append(BAND[b])
    continue
 print(' I found '+str(len(BAND1))+' bands of Fermi surface')   
 return BAND1  
##########################################################################################
##########################################################################################
def read_bands(name,STATES,EF):
 
 BAND1=[]
 COLOUR1=[]

 tree = ET.parse(name+'.save/atomic_proj.xml')
 root = tree.getroot()
 nk=int(root.find('HEADER/NUMBER_OF_K-POINTS').text)
 nb=int(root.find('HEADER/NUMBER_OF_BANDS').text)
 nwfc=int(root.find('HEADER/NUMBER_OF_ATOMIC_WFC').text)
 print(nk,nb,nwfc)
 COLOUR=[[ [] for j in range(nwfc) ] for k in range(nk)]
 COLOUR2=[[ [] for j in range(nb) ] for k in range(nk)]
 for i in range(nk):
  for j in range(nwfc):
    tmp=root.find('PROJECTIONS/K-POINT.'+str(i+1)+'/ATMWFC.'+str(j+1)).text
    tmp=tmp.replace(',',' ').replace('\n',' ').split()
    tmp=[ [tmp[2*b],tmp[2*b+1]] for b in range(nb)]
#    tmp=[ float(m[0]) for m in tmp]
    tmp=[ (float(m[0])**2+float(m[1])**2) for m in tmp]
#    COLOUR[i][j]=[ sum(tmp[st[0]:st[1]+1]) for st in STATES]
    COLOUR[i][j]=tmp
  COLOUR[i]=np.transpose(np.array(COLOUR[i])) # [nk][nb][nwfc]
  for k in range(nb): 
   COLOUR2[i][k]=[ sum(COLOUR[i][k][st[0]:st[1]+1]) for st in STATES]
 
 BAND=[[ ] for k in range(nk)]
 for i in range(nk):
  tmp=root.find('EIGENVALUES/K-POINT.'+str(i+1)+'/EIG').text
  tmp=tmp.replace(',',' ').replace('\n',' ').split()
  BAND[i]=[ float(m)*Hartree_to_eV/2. for m in tmp ]

 [BAND,COLOUR]=[np.transpose(np.array(BAND)),np.transpose(np.array(COLOUR2))]
 ns=len(COLOUR)
 print(ns)
 COLOUR1=[ [] for i in range(ns)]
#search  bands  crossing the fermi energy
 for b in range(len(BAND)):
  if min(BAND[b])<EF and max(BAND[b])>EF:
    BAND1.append(BAND[b])
    for i in range(ns):
      COLOUR1[i].append(COLOUR[i][b])
    continue
 print(' I found '+str(len(BAND1))+' bands of Fermi surface')
 print(' I have '+str(len(COLOUR1))+' colours')
 return [BAND1,COLOUR1]
##########################################################################################
'''
def read_bands(name,STATES,EF):

 BAND=[]
 BAND1=[]
 BAND2=[]
 COLOUR=[]
 COLOUR1=[]

 h=open(name+'.save/atomic_proj.xml','r')
 s=0
 for line in h:
  if '<EIGENVALUES>' in line: s=2
  elif '<PROJECTIONS>' in line: s=1

  elif s==2:
   if '<K-PO' in line:
    BAND.append([])
   elif 'EIG' not in line and '/K-PO' not in line:
    BAND[-1].append(float(line.split()[0])*Hartree_to_eV/2.) 


  elif s==1:
   if '<K-PO' in line:
    COLOUR.append([])
   elif '<ATM' in line:
    COLOUR[-1].append([])
   elif '</ATM' in line or '</K-P' in line:
    continue
   elif '</PROJ' in line:
    break
   else: 
#    COLOUR[-1][-1].append( abs(float(line.replace(',',' ').split()[0])) )
    COLOUR[-1][-1].append(sum([float(wfc)**2 for wfc in line.replace(',',' ').split()])**0.5)
 h.close()


 BAND=np.transpose((BAND))
 COLOUR2=[]
 for i in COLOUR:
  tmp=np.transpose(i)
  COLOUR2.append([])
  for k in tmp:
   COLOUR2[-1].append([])
   for j in STATES:
    COLOUR2[-1][-1].append(sum( k[j[0]:j[1]+1] ))
 COLOUR=np.transpose(COLOUR2)


 print(str(len(COLOUR))+' '+str(len(COLOUR[0]))+' '+str(len(COLOUR[0][0]))+' '+str(len(BAND))+' '+str(len(BAND[0])) )

 ns=len(COLOUR)
 COLOUR1=[ [] for i in range(ns) ]
#search  bands  crossing the fermi energy
 for b in range(len(BAND)):
  if min(BAND[b])<EF and max(BAND[b])>EF:
    BAND1.append(BAND[b])
    for i in range(ns):
      COLOUR1[i].append(COLOUR[i][b])
    continue
 print(' I found '+str(len(BAND1))+' bands of Fermi surface')

 BAND1=np.array(BAND1)

 print(str(len(COLOUR1))+' '+str(len(COLOUR1[0]))+' '+str(len(COLOUR1[0][0]))+' '+str(len(BAND1))+' '+str(len(BAND1[0])) )

 return [BAND1,COLOUR1]
'''




