## How to run GTEx_check on CCHMC cluster?

Let's first talk about single query.

In order to check the GTEx expression level (PSI) for a single splicing event. For example, `KYAT3:ENSG00000137944:E4.1-I4.1_88965413`, all you need to do is in your terminal, **remember, no need to do any module loading, as long as you have salomonis2 drive access, you should be fine**:

```
/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/ -e KYAT3:ENSG00000137944:E4.1-I4.1_88965413 -c 0.1 -m savage -t all -p True -d /data/salomonis2/LabFiles/Frank-Li/python3/data -o /data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/GTEx_viewer
```

A full prompt:

```
Options:
-e --event : Splicing event you want to interrogate
-c --cutoff : cutoff value for PSI
-m --mode : load off-the-shelf GTEx data or generating it from scratch
-t --tissue : tissue-specific abundance you want to inspect
-p --plot : Do you want to plot all tissue-specific expression bar charts?
-d : data folder path
-o : output folder path
-h --help : check help information 
Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati
```

The most important parameter you make sure to modify is the output folder, it determine where the output result, including a html page would go into.


If you instead want to know the raw read count information for a single spicing event, run following:

```
/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/queryGTEx_counts.py -e KYAT3:ENSG00000137944:E4.1-I4.1_88965413 -c 5 -o /data/salmonis2/LabFiles/Frank-Li/mhc/test_AltNeo/GTEx_viewer_count -d /data/salomonis2/LabFiles/Frank-Li/python3/data -m savage -t all -p True
```

A full prompt:

```
Options:
-e --event : Splicing event you want to interrogate
-c --cutoff : cutoff value for read counts
-o : output folder
-d : data folder
-m --mode : load off-the-shelf GTEx data or generating it from scratch
-t --tissue : tissue-specific abundance you want to inspect
-p --plot : Do you want to plot all tissue-specific expression bar charts?
-h --help : check help information 
Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati
```

Again, make sure you change the output path.



Now let's talk about multiple query:

I've already made a shell script for job submission, the path is:

```
/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/sub.sh
```

Please copy it to your folder instead of modifying it directly. Several things to change in your specific case:

1. `int` variable to your input file. Input file needs to have a unique column named "UID". It should be a tab delimited file.
2. `task` give your task a name, convenient to track further
3. `out` output folder, make sure do not include forward slash in the end
4. `HLA` doesn't matter, you can leave it untouched
5. `-m` this `mode` parameter needs to be set to `GTEx_check`
6. `-C` determine how many cores to use, it should be consistent with your shell script `#BSUB -n` parameter, but it doesn't matter for these task, you can set it to `1`.

No matter how big your input file is, the program needs to take **40** minutes or so to load GTEx dataset I processed previously. As long as the dataset being loaded into the memory, following queries would be fast. If you have >100K events to query, this is the advantage of using this program. The single query use HDF5 to allow quick access, the multiple query (pbz2 file) is designed for quering tens of thousands events.

To understand the output files:

Terminilogy:

`True`: events that do not express in GTEx at an appreciable level
`Unknown`: events that are not detected in GTEx
`False`: events that do not pass the GTEx check
`train`: `True` + `False`
`final`: `True` + `Unknown`
`whole`: `True` + `False` + `Unknown`

1. `singleSample_check.txt`: `final` file that pass the GTEx test, will be used for further step in the program for in-silico translation and neoantigen prioritization. 

2. `df_train.txt`: files will have IW tumor specicifity score.



