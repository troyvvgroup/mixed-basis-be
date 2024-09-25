#!/bin/bash

file="trial"
geoms="4_poly"
outer_basis="cc-pVDZ"
inner_be_levels="be1"
outer_be_levels="be2"
nproc="1"
runpath=$PWD
geompath=$runpath/../
code=$runpath/parallel_split_code.py
charge='0'
chempotopt=False

for g in $geoms
do
for o in $outer_basis
do
mkdir -p ${g}_${o}
cd ${g}_${o}

for ibe in $inner_be_levels
do
for obe in $outer_be_levels
do
be=${ibe}_in_${obe}
mkdir -p $be
cd $be

geom=$geompath/$g.xyz
cat > submit_${g}_${be}.sh << eof
#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -o nsbatch_${file}.out
#SBATCH -e sbatch_${file}.err
#SBATCH -c $nproc
#SBATCH -N 1
#SBATCH --mem=8000
#SBATCH --partition=high

export PYTHONPATH=/home/lweisbur/source/quemb:$PYTHONPATH

python3 $code $geom $obe $ibe $o $nproc $charge $chempotopt > out_file_${g}_${be}_${nproc}_${file}.out
eof

sbatch submit_${g}_${be}.sh
cd ../
done
done
cd ..
done
done
