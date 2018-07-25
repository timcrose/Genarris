
for i in `seq 1 56`
do
    sbatch sub_single_rcd_folder.sh
    echo $i
done
