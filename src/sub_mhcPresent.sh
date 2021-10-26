#BSUB -W 5:00
#BSUB -M 80000
#BSUB -n 1
#BSUB -J sub_mhcPresent
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


# you need to make sure necessary packages has been pre-installed in the following virtual enviroment
cd $(pwd)
module load anaconda3
conda init --all bash
source activate /data/salomonis2/LabFiles/Frank-Li/python3/mhcEnv

# get arguments
int='/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/demo.ori.txt'
task='demo'
out='/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast'
#HLA='HLA-A*2301,HLA-A*0301,HLA-B*4402,HLA-B*3501,HLA-C*0704,HLA-C*0401'
HLA='HLA-A23:01,HLA-A03:01,HLA-B44:02,HLA-B35:01,HLA-C07:04,HLA-C04:01'


#  /users/ligk2e/.local/share/mhcflurry/4/2.0.0/models_class1_presentation/

for k in {8..11..1}
do
python3 /data/salomonis2/LabFiles/Frank-Li/python3/mhcPresent.py -i ${int} \
    -t ${task} \
    -o ${out} \
    -d /data/salomonis2/LabFiles/Frank-Li/python3/data \
    -k $k -H ${HLA} \
    -s /data/salomonis2/LabFiles/Frank-Li/python3/netMHCpan-4.1/netMHCpan -M MHCI \
    -m Neoantigen -C 5 -c True \
    --cutoffPSI 0.1 --cutoffSample 0.1 --cutoffTissue 0.1
done


python3 /data/salomonis2/LabFiles/Frank-Li/python3/get_combined.py --task ${task} \
    --outdir ${out} \
    --HLA ${HLA}

python3 /data/salomonis2/LabFiles/Frank-Li/python3/getAugment.py --task ${task} \
    --outdir ${out} \
    --HLA ${HLA}




