#!/bin/bash

geoms="2 4"
outer_basis="cc-pVDZ"
inner_be_levels="be1"
outer_be_levels="be1 be2"

nproc="4"
partition="short"
mem0="4000"
mem=`echo $mem0 $nproc | awk '{print $1*$2}'`
runpath=$PWD
geompath=$runpath/../xyz/
code=$runpath/parallel_mb_code.py
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
#SBATCH -o nsbatch_${g}_${be}.out
#SBATCH -e sbatch_${g}_${be}.err
#SBATCH -c $nproc
#SBATCH -N 1
#SBATCH --mem=$mem
#SBATCH --partition=$partition

export PYTHONPATH=/home/lweisbur/source/quemb:$PYTHONPATH

python3 $code $geom $obe $ibe $o $nproc $charge $chempotopt > out_file_${g}_${be}_${nproc}.out
eof

sbatch submit_${g}_${be}.sh
cd ../
done
done
cd ..
done
done
