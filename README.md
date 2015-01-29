ScaffMatch README page -- Georgia State University -- December 2014

ScaffMatch v0.9 Igor Mandric - Alex Zelikovsky

email: mandric.igor@gmail.com


1) Software Description:

ScaffMatch is a novel scaffolding tool based on Maximum-Weight Matching able to produce high-quality scaffolds from NGS data (reads and contigs). The tool is written in Python 2.7. It also includes a bash script wrapper that calls aligner in case one needs to first map reads to contigs (instead of providing .sam files).

The arguments accepted by ScaffMatch are:

  -w) Working directory -- this is the directory where ScaffMatch files are stored. These are .sam files produced after mapping reads to contigs and the resulting scaffolds file `scaffolds.fa` fasta file;

  -c) Contig fasta file;

  -m) Command line argument with no options. It is used when .sam files are used instead of reads .fastq files. Do not use this option if you provide reads files;

  -1) (Comma separated list of) either .fastq or .sam file(s) corresponding to the first read of the read pair;

  -2) (Comma separated list of) either .fastq or .sam file(s) corresponding to the second read of the read pair;

  -i) (Comma separated list of) insert size(s) of the library(-ies);

  -s) (Comma separated list of) library(-ies) standard deviation(s) of insert size(s);

  -t) Bundle threshold. Pairs of contigs supported by number of read pairs less than the value of this argument are discarded. Optional argument, by default it is equal to 5;

  -g) Matching heuristics: use `max_weight` for Maximum Weight Matching heuristics with the Insertion step, use `backbone` for Maximum Weight Matching heuristics without the Insertion step, use `greedy` for Greedy Matching heuristics;

  -l) Log file - where to store the logs. Optional argument. By default, stdout is used.


One can use directly scaffmatch.py Python script when using .sam files.


2) Requirements: 

    * Python 2.7.x
    * Bash >= 4
    * Networkx >= 1.7
    * Numpy >= 1.6.2
    * Bowtie2

3) Scaffolding algorithm:

    3.1) Algorithm Overview:

        ScaffMatch algorithm consists of the following main steps:

            1. Mapping reads to contigs - optional.
            2. Constructing the scaffolding graph.
            3. Maximum Weight Matching step - producing the backbone scaffolds.
            4. Insertion step - inserting singletone contigs into the backbone.
            5. Writing the final scaffolds.fa file.


    3.2) Algorithm step by step:

        1. We use bowtie2 to map reads to contigs.

        2. The scaffolding graph G = (V, E) is constructed as follows: each vertex of the scaffolding graph G corresponds to one of the contig strands and each inter-contig edge corresponds to a bundle of read pairs connecting two strands of different contigs. The weight of an inter-contig edge is equal to the size of the corresponding bundle. Also for each contig we have a dummy edge connecting its two strands.

        3. In our interpretation, the Scaffolding Problem is reduced to the problem of finding the Maximum Weight Matching in the scaffolding graph G. We use either the well-known blossom algorithm (implemented in Networkx library) or a greedy O(N * log N) heuristic. After the matching is found, we obtain the so-called backbone scaffolds. 
        
        4. After the matching step, an insertion of singletones into the backbone is performed. It helps to increase the number of correct contig joins. The usefulness of this step is demonstrated in the corresponding publications of the authors.

        5. We write the scaffolds as a .fasta file. The gaps are filled with 'N's.
