
unset log
plot '2014-06-17/noproprep/msdave1.dat' u 3:4 title 'No Interaction' w lines, \
    '2014-06-17/prop/msdave1.dat' u 3:4 title 'Propulsion' w lines, \
    '2014-06-17/rep/msdave1.dat' u 3:4 title 'Repulsion' w lines, \
    '2014-06-19/adhesion/adhonly/msdave1.dat' u 3:4 title 'Adhesion Only' w lines, \
    '2014-06-17/proprep/msdave1.dat' u 3:4 title 'Propulsion Repulsion' w lines, \
    '2014-06-19/adhesion/prop/msdave1.dat' u 3:4 title 'Adhesion Propulsion' w lines, \
    '2014-06-19/adhesion/rep/msdave1.dat' u 3:4 title 'Adhesion Repulsion' w lines,\
    '2014-06-19/adhesion/proprep/msdave1.dat' u 3:4 title 'Adhesion Propulsion Repulsion' w lines

set title 'Single Species, Varying Properties'
set key left top

set term png
set output 'msd.png'
replot

set output 'msdlog.png'
set key bottom right
set log
replot

