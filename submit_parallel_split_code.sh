#!/bin/bash

geoms="protein complex"
outer_basis="def2-SVP def2-TZVP"
be_levels="be1"
nproc="24"
runpath=$PWD
geompath=$runpath/../geom
code=$runpath/parallel_split_code.py
charge='1'

for g in $geoms
do

for o in $outer_basis
do

mkdir -p ${g}_${o}
cd ${g}_${o}

for be in $be_levels
do
mkdir -p $be
cd $be

geom=$geompath/$g.xyz

cat > submit_${g}_${be}.sh << eof
#!/bin/bash

#SBATCH -t 1500:00:00
#SBATCH -o nsbatch.out
#SBATCH -e sbatch.err
#SBATCH -c $nproc
#SBATCH -N 1
#SBATCH --mem=16000
#SBATCH --partition=veryhigh

export PYTHONPATH=/home/q4bio/source/mol-be:$PYTHONPATH

python3 $code $geom $be $o $nproc $charge > out_file_${g}_be1_in_${be}_${nproc}.out
eof

sbatch submit_${g}_${be}.sh
cd ../
done
cd ..
done
cd ..
done
