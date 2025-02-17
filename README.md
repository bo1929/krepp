# krepp
A k-mer-based maximum likelihood method for estimating distances of reads to genomes and phylogenetic placement.

- [Installing krepp or compiling from source.](#installing-krepp-or-compiling-from-source)
- [Building an index from multiple reference genomes](#building-an-index-from-multiple-reference-genomes)
- [Estimating distances of reads to all reference genomes](#estimating-distances-of-reads-to-all-reference-genomes)
- [Phylogenetic placement of reads on a backbone tree](#phylogenetic-placement-of-reads-on-a-backbone-tree)
- [Practical distance estimation by sketching a file](#practical-distance-estimation-by-sketching-a-file)
- [A toy example for testing](#a-toy-example-for-testing)

## Quick start
### Installing krepp or compiling from source
Pre-compiled binaries are only available for Linux and x64 architecture, see the [latest release](https://github.com/bo1929/krepp/releases/tag/v0.0.4).

To compile from the source, simply clone the repository with its submodules and compile with
```bash
git clone --recurse-submodules -j8 https://github.com/bo1929/krepp.git
cd krepp && make
```
and run `./krepp --help`. Then, perhaps, copy it to a directory you have in your `$PATH` (e.g., `cp ./krepp ~/.local/bin`).

You may not have `libcurl`  on your system, and in that case, compilation would fail. You can install it (e.g., `sudo apt install curl`), and try running make again.
Otherwise, you could simply compile without libcurl by running `make WLCURL=0`.
Note that, without `libcurl`, you cannot use URLs to retrieve FASTA/FASTQ files directly from FTP servers, and all query/reference files will have to be stored locally.

### Building an index from multiple reference genomes
Given a set of reference genomes and a backbone tree, krepp can build an LSH index with colored k-mers by running
```bash
krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS 
```
where
* `$NUM_THREADS` is simply the number of threads,
* `-o` is the directory where the index will be stored,
* `-i` is a tab-separated file where the first column is for the IDs of references, i.e., leaves of the tree, and the second column is for the paths (or the URLs) to corresponding reference genomes (e.g., FASTA),
* `-t` is the optional backbone tree in Newick format matching the genome IDs (i.e., leaf labels) given (must be rooted and with labeled internal nodes).

Run `krepp build --help` to see more options and details.
If you don't have a tree and if you are mainly interested in estimating distances, run `krepp index` without the `-t` option (with the side effect of potentially making the color index slightly larger).
An index constructed without a backbone cannot be used for phylogenetic placement, but distance estimates stay the same.

Instead of building an index from scratch, you could use one of the public ones.
Currently, only two such indexes are available for microbial and archaeal genomes; both allow phylogenetic placement.
Note that the smaller index (v1) was built using a reference set which is a subset of the larger one.
Therefore, these indexes overlap, and you could just pick the one that you can afford in terms of memory available on your machine.

* Web of Life - v2 (15,493 archaeal and bacterial genomes): [index](https://ter-trees.ucsd.edu/data/krepp/index_WoLv2-k29w35-h14.tar.gz), [tree](https://ter-trees.ucsd.edu/data/krepp/misc/backbone_tree-WoLv2.nwk.gz), [metadata](https://ter-trees.ucsd.edu/data/krepp/misc/metadata-WoLv2.tsv.gz)
* Web of Life - v1 (10,576 archaeal and bacterial genomes): [index](https://ter-trees.ucsd.edu/data/krepp/index_WoLv1-k29w35-h14.tar.gz), [tree](https://ter-trees.ucsd.edu/data/krepp/misc/backbone_tree-WoLv1.nwk.gz), [metadata](https://ter-trees.ucsd.edu/data/krepp/misc/metadata-WoLv1.tsv.gz)

The genome IDs that will be reported with these indexes themselves may not be informative.
You can use the provided metadata files to analyze your distance estimates further; perhaps for taxonomic classification or abundance profiling.
Similarly, the backbone tree could be used for UniFrac computation.

We note that memory use increases almost linearly with `$NUM_THREADS` for the `index` subcommand, and hence, you may want to decrease it if you run out of memory.
This is not the case for `dist`, `place`, and `seek` commands, feel free to use all cores available during the query time.

Another way of making the index smaller is increasing the minimizer window size (`-w`) with respect to `-k`.
This may cost you some accuracy as fewer *k*-mers will be indexed, but shouldn't be an issue as long as it is not too aggressive (e.g., `w-k>9`).

The peak memory usage during the index construction could be also reduced to the memory level available by partitioning the index into smaller pieces.
This is done by a variation of FracMinHash, controlled by options `-m` and `-r`.
`krepp` partitions the index into `-m` (more or less) equally sized pieces, and these partitions could be built independently, but one can query them together.
The `-r` option determines the partition that is going to be constructed: if `--no-frac` is given only the `-r`th partition, otherwise all partitions from 0th to `-r`th, will be constructed and saved (see `krank index --help` for details).
You don't need to construct all partitions, `krepp` will search in whatever is available, and these partitions can be distributed independently.
The default is `-m 5 -r 1 --frac`, so 40% (0th and 1st of 5 partitions) of the minimized *k*-mers will be indexed.
The only requirement is keeping the `-m` value (and of course `-i` and `-t`) fixed across all partitions.
For instance, one can index 3% percent of the reference *k*-mers and construct a lightweight index by running:
```bash
krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS -m 100 -r 99 --no-frac
krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS -m 100 -r 98 --no-frac
krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS -m 100 -r 97 --no-frac
```
or, alternatively
```bash
krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS -m 100 -r 2 --frac
```
This would be faster, but with increased memory use during the index construction, despite resulting in an index of the same size.
Perhaps later, if you don't think this works well for your task, you can build the remaining 47% of the index and complete it to 50% by running
```bash
krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS -m 100 -r 46 --frac
```
or maybe, if you have a small memory machine, you can iterate over partitions
```bash
seq 1 46 | xargs -I{} -P2 bash -c "krepp index -o $INDEX_DIR -i $INPUT_FILE -t $BACKBONE_NEWICK --num-threads $NUM_THREADS -m 100 -r {} --no-frac"
```

Just stick to the defaults if everything works or if you don't run out of memory.

### Estimating distances of reads to all reference genomes
Once you have the index built, query reads against it to get distance estimates is quite simple:
```bash
krepp dist -i $INDEX_DIR -q $QUERY_FILE --num-threads $NUM_THREADS -o ${QUERY_NAME}.tsv
```
where `-q` is the path (or URL) of a FASTA/Q file containing query sequences, and `krepp` simply outputs everything to stdout and writes log messages to stderr.
The output is in a tab-separated format where the first column stands for the sequence ID, the second column is the ID of the reference matching, and the third column is the distance estimate.

### Phylogenetic placement of reads on a backbone tree
In addition to distance estimation, one could place reads on the backbone tree given as input while building the index
```bash
krepp place -i $INDEX_DIR -q $QUERY_FILE --num-threads $NUM_THREADS -o ${QUERY_NAME}.jplace
```
where the output is in `jplace` format (version 3) (see the description [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009)).
We leave the `p` field empty for sequences without any *k*-mer match.
Some down-stream analysis tools that takes a `jplace` file as the input may not be compatible with empty placement fields (such as [`gappa`](https://github.com/lczech/gappa)), in this case you can simply filter those lines (e.g., by parsing `jplace` as a JSON file in Python or via bash scripting `grep -v "\[ \]" ${QUERY_NAME}.jplace | sed -z "s/},\n\t]/}\n\t]/g" > ${QUERY_NAME}-filtered.jplace`).

### Practical distance estimation by sketching a file

In addition to indexing multiple references together, `krepp` can also create a sketch from a FASTA/Q file for practical analysis with respect to a single reference.
Simply run
```bash
krepp sketch -i $INPUT_FILE -o $SKETCH_PATH --num-threads $NUM_THREADS
```
where `-i` is an URL or a filepath containing reference sequences from which *k*-mers will be extracted, and `-o` is the path that the resulting skecth will be saved in.
A sketch can not be used for phylogenetic placement, but you can efficiently seek query sequences in a sketch to get distance estimates by running
```bash
krepp seek -i $SKETCH_PATH -q $QUERY_FILE --num-threads $NUM_THREADS -o ${QUERY_NAME}.tsv
```

The output is in a tab-separated format with two columns: *i)* the sequence ID and *ii)* the distance estimate.

### A toy example for testing

#### Constructing a small index
You can build it from scratch consisting of only 25 genomes provided in `test/` to make yourself familiar with `krepp`.
```bash
cd test/
tar -xvf references_toy.tar.gz && xz -d references_toy/*
../krepp index -h 11 -k 27 -w 35 -o index_toy -i input_map.tsv -t tree_toy.nwk --num-threads 8
```
This command took only a couple of seconds and used $<$1.5GB memory for 6,975,500 indexed *k*-mers on a machine with Intel Xeon Silver 4110 CPUs.
The resulting index will be stored in `index_toy`.
Alternatively, you could download one of the larger public libraries to make it more realistic and use it also for your novel query sequences.

#### Querying novel sequences against the reference index
Once you have your index (e.g., the one we built above: `index_toy`), you can estimate distance by running:
```bash
../krepp --num-threads 8 dist -i index_toy -q query_toy.fq | tee distances_toy.tsv
```
The first five lines of `distances_toy.tsv` are going to look like:
```
#software: krepp	#version: v0.0.4	#invocation :../krepp --num-threads 8 dist -l index_toy -q query_toy.fq
SEQ_ID	REFERENCE_NAME	DIST
||61435-4122	G000341695	0.0898062
||61435-4949	G000830905	0.147048
||61435-4949	G000341695	0.0740587
||61435-4949	G000025025	0.131182
||61435-4949	G000741845	0.0395985
```

Quite similarly, you can place reads by running:
```bash
../krepp --num-threads 8 place -i index_toy -q query_toy.fq | tee placements_toy.jplace
```

The resulting placement file is a JSON file in a special format called `jplace`:
```bash
head -n15 placements_toy.jplace
```
```
{
	"version" : 3,
        "fields" : ["edge_num", "likelihood", "like_weight_ratio", "placement_distance", "pendant_length", "distal_length"],
	"metadata" : {
		"software" : "krepp",
		"version" : "v0.0.4",
		"repository" : "https://github.com/bo1929/krepp",
		"invocation" : "../krepp --num-threads 8 place -l index_toy -q query_toy.fq"
	},
	"placements" :
		[
			{"n" : ["||61435-4122"], "p" : [[39, 0, 0.00108799, -12.2542, 0, 0.0257952]]},
			{"n" : ["||61435-4949"], "p" : [[38, 0, 0.000709545, -35.0081, 0.96283, 0.0362437]]},
			{"n" : ["||61435-317"], "p" : [[35, 0, 0.00585016, -28.701, 6.56765, 0.0363296]]},
			{"n" : ["||61435-2985"], "p" : [[38, 0, 0.000709545, -25.4031, 0.0492549, 0.0114681]]},
```

You can proceed with your down-stream analysis using other tools, such as [`gappa`](https://github.com/lczech/gappa):
```bash
# filtering sequences without any match and hence an empty placement field
grep -v "\[ \]" placements_toy.jplace | sed -z "s/},\n\t]/}\n\t]/g" > placements_toy-filtered.jplace
# generating a colored tree based on placement densities across the backbone tree
gappa examine heat-tree --jplace-path placements_toy-filtered.jplace --write-svg-tree
```
