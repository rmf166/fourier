rm *.dat *.p *.pdf
../src/fa.exe
gnuplot plot*.p
pdfunite *DD.pdf *SC.pdf *LD.pdf *LC.pdf final.pdf
rm plot*.pdf
