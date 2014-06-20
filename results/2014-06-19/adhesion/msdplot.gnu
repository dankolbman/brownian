
unset log
plot 'prop/msdave1.dat' u 3:4 title 'With Propulsion' w lines, \
    'rep/msdave1.dat' u 3:4 title 'With Repulsion' w lines,\
    'proprep/msdave1.dat' u 3:4 title 'With Propulsion and Repulsion' w lines, \
    'adhonly/msdave1.dat' u 3:4 title 'Adhesion Only' w lines

set title 'Adhesive, Single Species'
set key left top

set term png
set output 'msd.png'
replot

set output 'msdlog.png'
set log
replot

