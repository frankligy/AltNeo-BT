#BSUB -W 3:00
#BSUB -M 80000
#BSUB -n 2
#BSUB -J sub_mhcPresent
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err



# get arguments
int='/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/TCGA-A1-A0SB.ori.txt'
task='test_run'
out='/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/test_AltNeo_T'
HLA='HLA-A11:01,HLA-A25:01,HLA-B35:01,HLA-B44:02,HLA-C04:01,HLA-C05:01'



/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/mhcPresent.py -i ${int} \
    -t ${task} \
    -o ${out} \
    -d /data/salomonis2/LabFiles/Frank-Li/python3/data \
    -k 9 -H ${HLA} \
    -s /data/salomonis2/LabFiles/Frank-Li/python3/netMHCpan-4.1/netMHCpan -M MHCI \
    -m GTEx_check -C 2 -c True \
    --cutoffPSI 0.1 --cutoffSample 0.1 --cutoffTissue 0.1






