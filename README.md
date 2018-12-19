Parallel version of TEtranscript using multi-process and global variable based shared memory
=======
Technically, TEtranscipt executes analysis for samples which is independent between each other, there is no need to develop 
a parallel version of TEtranscript since the concurrency execution can be simply accomplished through running multiple script at same time. 
However, the naive parallelsim can't be done properly because of its huge memory cost will break the limitation for single server. This is the main reason we develop this version of script.

Python didn't support "true" multi-thread on multi-core architectures. To increase the concurrency of this script, we should using multi process method instead of multi-thread. 
To overcome the memory limitation, large size variables for analysis should be shared between multiple process instead of copied to different processes. 

In this script, large memory cost variables such as teIdx, geneIdx, and intronIdx was declared as a global memory, 
since these variables is read only during the life cycle of the program, it is possible to make it 
accessible to all sub-processes without extra memory cost. Because Linux OS will automatically create a copy of these 
variables in the new process with a copy-on-write model,

In this script, we make a fuction named subprocessWorker, this function is used for analysis for different treatments with same references. 
The function subprocessWorker() accept tfile as its only argument. All other data which is required for function count_reads will be stored at global variables as mentiond in above paragraph. 

To run this code please do following:

```
python ./bin/TEParallel.py --t tfilelist -np 3 --GTF gtf/GKN1/GKN1.refGene.gtf --TE gtf/GKN1/GKN1.rmsk_TE.gtf --intron gtf/GKN1/GKN1.intron_TE_interval.gtf
 --exonTE gtf/GKN1/GKN1.rmsk_TE_exonic.gtf --intergenicTE gtf/GKN1/GKN1.rmsk_TE_intergenic.gtf --project GKN1test  
```
where -np is the number of child process you want to spawn (-np equal to 3 means there is 3 child process used for processing, the total number of python process
will be 4)

tfilelist is a file contains all treatment files location, I uploaded a sample file in this porject. 

Please be aware that I made a ./result folder to store all the results, and I did some modificiation on the result file name to make it related to different input files. 
Please correct it according to your process pipeline standard. 



