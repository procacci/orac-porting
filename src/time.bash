#!/bin/bash
# do timing on PAR dirs 
if [ $# -ne 3 ]; then 
    echo "Syntax: time.bash Orac.inputfile i1 i2"
    echo ""        
    echo "        where i1, i2 are the first digits of lower bound (in Angs)"
    echo "        of the l and h intermolecular shells, respectively"
    echo "        "
    echo "        Example: "
    echo "        $ .\time.bash orac.in 5 8 "
    exit 
fi
rm -fr ?PAR
for i in 1 2 4 8 ; do
    mpiexec  -n $i ../../../src/Intel_PARALLEL/orac_Linux < $1 ; mv PAR0001 $i.PAR 
done
rm timing.$1 time.*PAR
grep step 1.PAR/fort.1000 | grep -v Tst | grep -v time | sed "s/=/#/g" > timing.$1
grep FLexible 1.PAR/fort.1000  | tail -1 >> timing.$1
for j in ?.PAR; do 
    for i in N0 N1 TIME\-0 TIME\-$2 TIME\-$3 MPI  ;do 
	grep "$i" $j/fort.1000 | awk -v shell=$i '{{a+=$6}}END {print  shell,a}' >> time.$j ; 
    done; 
done
echo "   Shell    N=1           N=2          N=4         N=8        s(1/2)         s(1/4)      s(1/8)  
" >> timing.$1
paste time.?.PAR | awk '{printf "%6s%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f\n",\
 $1, $2,$4,$6,$8 ,$2/$4,$2/$6,$2/$8 }' >> timing.$1

