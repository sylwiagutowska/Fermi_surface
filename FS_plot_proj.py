#!/usr/bin/env mayavi2
import numpy as np
import math as m
import os
#from vpython import *
import random
from scipy.interpolate import griddata
import mayavi
from mayavi import mlab

Hartree_to_eV=27.21138602

try: ifcolors=raw_input('Coloured FS? [y]/n ')
except: ifcolors=input('Coloured FS? [y]/n ')
if ifcolors!='n': ktorycolor=input('Which states you want to see as colours (which line from sum)?') 
print('Got it!')

filename = 'FS.dat'
f_plot = open(filename, 'r')
plotData = np.loadtxt(f_plot)
[nkp,nb]=[int(i) for i in plotData[0,0:2]]
print([nkp,nb])
EF=plotData[0,2]
'''
maxx=max([ m for m in plotData[:,0]])
maxz=max([ m for m in plotData[:,2]])
for i in range(2,len(plotData)):
 if  plotData[i,1]!=plotData[i-1,1]: 
  for j in range(nb): plotData[i,3+j]=EF*0.9
  for j in range(nb): plotData[i-1,3+j]=EF*0.9
 elif  plotData[i,0]==maxx:
  for j in range(nb): plotData[i,3+j]=EF*0.9
 elif  plotData[i,2]==maxz:
  for j in range(nb): plotData[i,3+j]=EF*0.9
'''
plotX = plotData[1:, 0]
plotY = plotData[1:, 1]
plotZ = plotData[1:, 2]
Energy=[plotData[1:, 3+b] for b in range(nb)]
R = plotData[1:, 0:3]
f_plot.close()
if ifcolors!='n': 
 filename = 'FS_colour'+str(ktorycolor)+'.dat'
 h= open(filename, 'r')
 plotData = np.loadtxt(h)
 Colours=[plotData[0:, b] for b in range(nb)]
 colormax=max([ max(i) for i in Colours])
 h.close()
 print (colormax)
print('Reading - done')

pts =48*1j
ptss=48
#dx=[(max(plotX)-min(plotX))/pts,(max(plotY)-min(plotY))/pts,(max(plotZ)-min(plotZ))/pts]
X, Y, Z = np.mgrid[min(plotX):max(plotX):pts, min(plotY):max(plotY):pts, min(plotZ):max(plotZ):pts]

F = [griddata(R, Energy[b],  (X, Y, Z), method='nearest') for b in range(nb)]
if ifcolors!='n': C = [griddata(R, Colours[b], (X, Y, Z), method='nearest') for b in range(nb)]


##############################


aa=mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))


##############edges of BZ
filename = 'FS_NN.dat'
f_plot = open(filename, 'r') 
NN=np.loadtxt(f_plot)
f_plot.close()
pts=mlab.plot3d(NN[:, 0],NN[:, 1],NN[:, 2],\
color=(0,0,0),tube_radius=0.01)




###############Fermi surface
if ifcolors!='n':
 sf=[]
 for iband in range(nb):
 
  src =  mlab.pipeline.scalar_field(X,Y,Z,F[iband])
  src.mlab_source.dataset.point_data.add_array(C[iband].T.ravel())
  src.mlab_source.dataset.point_data.get_array(1).name = 'band'+str(iband)
  src.update()

  kontur = mlab.pipeline.contour(src)
  kontur.filter.contours= [EF,]
  mlab.pipeline.set_active_attribute(kontur,
                                    point_scalars='band'+str(iband))
  sf.append(mlab.pipeline.surface(kontur, colormap='blue-red',opacity=1 ))# ,vmin=0))#,vmax=colormax))
else:
 FS=[mlab.contour3d(X, Y, Z, F[iband], contours=[EF] ) for iband in range(nb)]
#mlab.axes()

print('Interpolating - done')


#### END OF BRILLOIN ZONE


#aa.scene.render_window.aa_frames = 8

if ifcolors!='n': 
 colorbar=mlab.colorbar (nb_labels=3,label_fmt='%.2f',orientation='horizontal')
#colorbar.scalar_bar_representation.proportional_resize=True
 colorbar.scalar_bar_representation.position = [0.31, 0.85]
 colorbar.scalar_bar_representation.position2 = [0.38, 0.15]
#colorbar.scalar_bar.unconstrained_font_size = True
 colorbar.label_text_property.bold = 1
# colorbar.label_text_property.font_size=172


# mlab.savefig('FS_'+str(iband)+'.png')
mlab.show()


