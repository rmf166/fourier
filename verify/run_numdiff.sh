rm diff.log
for file in result*s1*.dat
do
    numdiff ${file} ../02012021/${file} >> diff.log 
done
