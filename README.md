### Introduction

deBGA-VARA is a novel variation-aware high-throughput sequencing read aligner. VARA is a variation-aware read global alignment algorithm. We developed deBGA-VARA by integrating VARA into deBGA. It is suitable for aligning various kinds of high-throughput sequencing reads to multiple similar reference genomes.

deBGA-VARA indexes the genome and the variation knowledge through de Bruijn graph framework. 

deBGA has outstanding throughput on reads mapping to genome than other conventional read aligners and variation-aware read aligners. A series of benchmarks on simulated and HTS datasets demonstrated that it can simultaneously achieve good throughput, sensitivity and accuracy in various kinds of read alignment tasks. DeBGA-VARA is open source and free for non-commercial use.

deBGA-VARA is mainly designed and developed by Hongzhe Guo in Center for Bioinformatics, Harbin Institute of Technology, China.

### Memory requirement
The memory usage of deBGA-VARA can fit the configurations of most modern servers and workstations. Its peak memory footprint depends on the length of reference genome, i.e., 40.32 Gigabytes for the real H. Sapiens, on a server with  with 2 Intel E5-2630v3 CPUs at 2.4 GHz (12 cores in total), 512 GB RAM and 48 TB hard disk space.

The wall time and memory footprint of the index construction for the references are almost as same as deBGA (https://github.com/HongzheGuo/deBGA). e.g, the memory footprints for GRCh37/hg19 is about 40 Gigabytes and the time is in about 4.5 hours (k-mer size is 22).

The memory footprint of deBGA-VARA when aligning the reads to genome is almost as same as deBGA.

### Installation

Current version of deBGA-VARA needs to be run on Linux operating system.  
The source code is written in C, and can be directly download from: https://github.com/HongzheGuo/deBGA-VARA  
The makefile is attached. Use the make command for generating the executable file.  

### Synopsis

deBGA index [options] reference.fasta \<index_route\>  
Index reference in RdBG-Index format  

deBGA aln [options] \<index_route\> \<single_end_read.fastq [pair_end_read1.fastq pair_end_read2.fastq]\> \<result_file.sam\>  
Align read to its primitive location in Reference  

### Parameters (could be updated in the future for adding new functions)
```
deBGA index   
--ext-alt 				STR	(default: not set) the index construction option for variation aware reference. When –ext-alt option is set, deBGA builds the index of reference and variation knowledge in source .vcf file. 

-k,                     INT the k-mer length of the vertices of RdBG. This is a basic parameter for building theRdBG-index. For the current version of deBGA, 
							the range of -k parameter is restricted to 21-28 bp, considering both of the effectiveness of the seeds and memory footprint[22]. 
 
 
deBGA aln 
--ext-alt-aln 				(default: not set) the local alignment option for VARA algorithm. When –ext-alt-aln option is set, deBGA carries out the VARA algorithm with variation index during extension. 

-k,                     INT the minimum length of a valid Uni-MEM seed. For current version of deBGA, this setting should be equal to the k-mer length of the RdBG-index[22].    

-s,                     INT the number of iterations of re-seeding. deBGA iteratively aligns a read in at most (-s + 1) iterations with various set of seeds. 
							This parameter works combining with the minimum interval of seeding (the -i option) and the maximum allowed number of hits per seed (the -n option). 
							That is, in the r-th iteration (r = 1 ,…, -s), deBGA tries to generate seeds at every ((-s)– r +1)*(-i) bp along the read. 
							If the read still cannot be successfully aligned after -s iterations, deBGA would ignore -n option to handle very repetitive reads in the (-s+1)-th iteration[4].    

-i,                     INT the minimum interval of seeding. This parameter determines the density of seeds, which is related to the sensitivity and efficiency of alignment. 
							Configuring this parameter with lower value will make deBGA generate seed more densely, which could improve the sensitivity, but at the expense of throughput[5].   

-n,                     INT the maximum allowed number of hits per seed. In the first -s iterations of the alignment process, the seeds with more than -n hits would be discarded for achieving faster speed.
							DeBGA ignores this restriction to introduce repetitive seeds if the read still cannot be successfully aligned after -s iterations[300].  
							
-c,                     NUM the threshold on the edit distance for early stop. In each iteration, deBGA checks the edit distance of the obtained best alignment. 
							If the ratio ED_best/RL <(-c), where ED_best and RL are respectively the edit distance of the best alignment and the read length, 
							deBGA considers that the read is confidently aligned and early-stops the alignment[0.05].    

--cl,                   NUM the adjusted threshold on the edit distance for early stop. When --cl option is set, in any given iteration, if there is at least one Uni-MEM seed available for extension, 
							but no successful alignment is obtained, the threshold on the edit distance(-c) can be dynamically adjusted to the value of --cl in next iterations. 
							This is a heuristic may acceleratethe alignment ofdivergent reads, e.g., reads having many low quality bases. 
							If --cl is not set, there will be no change on the -c option during the process (default: following the setting of -c).    

--local,                    the local alignment option for confident alignment. When --local option is set, in any given iteration, if there is at least one Uni-MEM seed available for extension, 
							but no successful alignment is obtained,deBGA perform local alignment instead of end-to-end alignment in following iterations.The best obtained local alignment will be output as the result. 
							It is also worthnoting that the --cl option and --local option should not be simultaneously set (default: not set).    

--local-match,          INT the score for a matched base in the local alignment.This option will take effect only if --local option is set[1].     

--local-mismatch,       INT the penalty for a mismatched base in the local alignment. This option will take effect only if --local option is set[4].    

--local-gap-open,       INT the penalty for a gap open in the local alignment.This option will take effect only if --local option is set[6].    

--local-gap-extension,  INT the penalty for gap extension in the local alignment.This option will take effect only if --local option is set[1].     

--stdout,					(default: not set) output alignments by stdout. This option will let deBGA directly output alignments by stdout instead of user defined file.

-u,                     INT the upper limit of insert size. For a pair-end read, deBGA pairs the alignments of the two ends according to the upper (-u option) and lower (-f option) limits of the insert size.
							deBGA will consider it as a suitable pair-end alignment only if the inferred insert size is within the range [-f, -u][1000].

-f,                     INT the lower limit of insert size. For a pair-end read, deBGA pairs the alignments of the two ends according to the upper (-u option) and lower (-f option) limits of the insert size. 
							deBGA will consider it as a suitable pair-end alignment only if the inferred insert size is within the range [-f, -u][50].        

-o,                     INT the maximum number of alignment output.deBGA outputs at most -o alignments for the read. This is except for the pair-end reads which are handled with the anchoring alignment strategy. 
							For thosereads, the number of outputsis determined by the -x option[20].  

-x,                     INT the maximum number of alignment output for anchoring alignment. For the pair-end reads aligned with the anchoring alignment strategy, deBGA will output at most -x alignments[150].    

-l,                     INT the maximum allowed read length. For the current version of deBGA, reads shorter than -l bp will be normally processed, andfor reads longer than -l bp, only the first -l bp will be aligned, 
							and the other parts will be trimmed. Set -l option with a larger number may slightly increasethememory footprint. For most nowadays next generation sequencing reads, e.g., reads from Illumina platforms, 
							the default setting is long enough to work withoutthe trimming. Moreover, the current version ofdeBGA can support reads upto 4096 bp (setting -l to 4096)[512].     

-e,                     INT the budget for single-end alignment. In single-end read alignment, deBGA sets a budget on the computation resource in advance for balancing the efficiency and the sensitivity. More precisely, in the extension phase, 
							deBGA subsequently extend the candidate seeds in order of their coverage lengths, until more than -e extension operations have been totally executed after handling some of the seeds, or all the seeds are extended[100].

-p,                     INT the number of threads. The current version of deBGA supports upto 32 threads in read alignment[1].    
```

### Quick start

Genome indexing:
deBGA index --ext-alt variation.vcf Reference Index_Dir

Read alignment:
deBGA aln --ext-alt-aln Index_Dir Fastq_File Sam_file

### Simulation benchmarking

We simulated two datasets from the human reference genome (GRCh37/hg19) and 200 MB VCF dataset (NA12878) through Mason Simulator (version0.1.2). We simulated Illumina-like pair-end reads with lengths of 100 bp and 250 bp and the mean and standard deviation of the insert size are respectively 500 bp and 25 bp. These datasets helped us to evaluate the performance of deBGA. The datasets have been uploaded to Google Drive, and can be downloaded through the following link:


### Reference

Fast variation-aware read alignment with deBGA-VARA.

### Contact

For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn; hzguo@hit.edu.cn

