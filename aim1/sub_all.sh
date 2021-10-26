#BSUB -W 5:00
#BSUB -M 196000
#BSUB -n 2
#BSUB -R "span[ptile=2]"
#BSUB -J run_python
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


./process_gtex_count.py







