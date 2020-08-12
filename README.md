# Fermi_surface
The program FS_read_proj.py:
- reads data from PREFIX.save directory (produced by Quantum Espresso)
- translates the data from IBZ to full Brillouin Zone using symmetry operations
- writes output ready to use by write_fermi_surface.py


The program FS_plot_proj.py
- reads file produced by read_fermi_surface.py
- vizualizes the Fermi surface using Mayavi library

How to use:
- calculate the band structure on the dense grid with help of pw.x (part of Quantum Espresso) using "nscf" type of calculations
- if you would like to colour the Fermi surface by the contribution of the atomic states, then you need run projwfc.x (part of Fermi Surface)
- run: python2 FS_read_proj.py
   *if you would like to use the data from projwfc.x, then answer 'y'
   *enter the path+PREFIX to *save directory
- if you would like to use the data from projwfc.x, then look on output of projwfc.x, look which states you would like to consider, create file "sum" and write down the number of states to consider. Each line fill with the states, which you would like to sum, eg 1-6 means, that you would like to sum up all states from 1 to 6.
- run: python3 FS_plot_proj.py
   *if you would like to use the data from projwfc.x, then answer 'y' and type which atomic states you would like to consider
- the Mayavi window appears with the Fermi surface. You can use wide menu of Mayavi: you can change colors, hide some parts of Fermi Surface, change Fermi energy etc.

