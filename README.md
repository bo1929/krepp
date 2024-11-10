# krepp
A k-mer-based maximum likelihood method for estimating distances of reads to genomes and phylogenetic placement.

## Quick start
### Installing and compiling krepp from source
Pre-compiled binaries are not available yet but hopefully will be soon.

To compile from the source, simply clone the repository with its submodules and compile with
```
git clone --recurse-submodules -j8 https://github.com/bo1929/krepp.git
cd krepp && make
```
and run `./krepp --help`. Then, perhaps, copy it to a directory you have in your `$PATH` (e.g., `cp ./krepp ~/.local/bin`).

### Building a krepp index from reference genomes
Given a set of reference genomes and a backbone tree, krepp can build an LSH index of colored k-mers by simply running
```bash
krepp --num-threads $NUM_THREADS build -l $INDEX_DIR -i $PATHS_MAPPING -t $BACKBONE_NEWICK
```
where
* `$NUM_THREADS` is simply the number of threads,
* `-l` is the directory where the index will be stored (preferably a unique name),
* `-i` is a tab-separated file where the first column is IDs of references, i.e., leaves of the tree, and the second column is for paths to reference genomes (e.g., FASTA),
* `-t` is the backbone tree in Newick format matching the genome IDs (i.e., leaf labels) given.

Run `krepp build --help` to see more options.
If you don't have a tree and are mainly interested in estimating distances, an arbitrary (or random) tree might just work as a workaround (with the side effect of potentially making the color index a little larger).

Instead of building an index from scratch, you could use one of the public ones.
Currently, there are only three such indexes available for microbial and archaeal genomes.
Note that smaller indexes were built using a reference subset of the larger one.
Therefore, these indexes overlap, and you could just pick the one that you can afford memory-wise.

* Web of Life - v2 (15,493 archaeal and bacterial genomes): [index](https://ter-trees.ucsd.edu/data/krepp/index_WoLv2-k29w35-h14.tar.gz), [tree](https://ter-trees.ucsd.edu/data/krepp/misc/backbone_tree-WoLv2.nwk.gz), [metadata](https://ter-trees.ucsd.edu/data/krepp/misc/metadata-WoLv2.tsv.gz)
* Web of Life - v1 (10,576 archaeal and bacterial genomes): [index](https://ter-trees.ucsd.edu/data/krepp/index_WoLv1-k29w35-h14.tar.gz), [tree](https://ter-trees.ucsd.edu/data/krepp/misc/backbone_tree-WoLv1.nwk.gz), [metadata](https://ter-trees.ucsd.edu/data/krepp/misc/metadata-WoLv1.tsv.gz)

Only the index is required to query novel sequences, but genome IDs are not informative.
You can use the provided metadata to analyze the distance estimates further for taxonomic classification or abundance profiling.
Similarly, the backbone tree could be used for UniFrac computation.

Lastly, we note that memory use increases almost linearly with `$NUM_THREADS` for the `build` subcommand, and hence, you may want to decrease it if you run out of memory.
The peak memory usage can always be reduced to the memory level available by splitting the LSH index into batches and building each batch separately by using `-m`, `-r` and `--no-frac` options (see relevant documentation and `krank --help`).
This is not the case for `dist` and `place` commands, feel free to use all cores available during the query-time.

### Estimating distances of reads to reference genomes
Once you have the index built, query reads against it to get distance estimates is quite simple:
```bash
krepp --num-threads $NUM_THREADS dist -l $INDEX_DIR -q $QUERY_PATH
```
where `-q` is the path of a FASTA/Q file containing query reads, and `krepp` simply outputs everything to stdout and writes log messages to stderr.
The output is in a tab-separated format in which the first column stands for the read ID, the second column is the ID of the reference matching, and the third column is the distance estimate of krepp.

### Phylogenetic placement of reads on the backbone tree
In addition to distance estimation, one could place reads on the backbone tree given as an input while building the index
```bash
krepp --num-threads $NUM_THREADS place -l $INDEX_DIR -q $QUERY_PATH
```
where the output is in `jplace` format (version 3) (see the description [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009)).
We allow the `p` field to be empty for sequences without any *k*-mer match.
Some down-stream analysis tools that takes a `jplace` file as the input may not be compatible with empty placement fields (such as [`gappa`](github.com/lczech/gappa)), in this case you can simply filter those lines (e.g., by parsing `jplace` as a JSON file in Python or via bash scripting `grep -v "\[ \]" $QUERY_NAME.jplace | sed -z "s/},\n\t]/}\n\t]/g" > $QUERY_NAME-filtered.jplace`).

### A toy example for testing

#### Constructing a small index
You can build it from scratch consisting of only 25 genomes provided in `test/` to make yourself familiar with `krepp`.
```bash
cd test/
tar -xvf references_toy.tar.gz && xz -d references_toy/*
../krepp --num-threads 8 build -h 11 -k 27 -w 35 -l index_toy -i input_map.tsv -t tree_toy.nwk
```
This command took less than 10 seconds and used 1.5GB memory for 6,975,500 indexed *k*-mers on a machine with Intel Xeon Silver 4110 CPUs.
The resulting index will be stored in `index_toy`.
Alternatively, you could download one of the larger public libraries to make it more realistic and use it also for your novel query sequences.

#### Querying novel sequences against the reference index
Once you have your index (e.g., the one we built above: `index_toy`), you can estimate distance by running:
```bash
../krepp --num-threads 8 dist -l index_toy -q query_toy.fq | tee distances_toy.tsv
```
The first five lines of `distances_toy.tsv` are going to look like:
```
#software: krepp	#version: v0.0.2	#invocation :../krepp --num-threads 8 dist -l index_toy -q query_toy.fq
SEQ_ID	REFERENCE_NAME	DIST
||61435-4122	G000341695	0.0898062
||61435-4949	G000830905	0.147048
||61435-4949	G000341695	0.0740587
||61435-4949	G000025025	0.131182
||61435-4949	G000741845	0.0395985
```

Quite similarly, you can place reads by running:
```bash
../krepp --num-threads 8 place -l index_toy -q query_toy.fq | tee placements_toy.jplace
```

The resulting placement file is a JSON file in a special format called `jplace`:
```bash
head -n15 placements_toy.jplace
```
```
{
	"version" : 3,
	"fields" : ["edge_num", "like_weight_ratio", "likelihood", "pendant_length", "distal_length"],
	"metadata" : {
		"software" : "krepp",
		"version" : "v0.0.1",
		"repository" : "https://github.com/bo1929/krepp",
		"invocation" : "../krepp --num-threads 8 place -l index_toy -q query_toy.fq"
	},
	"placements" :
		[
			{"p" : [[39, 0.000000, 20.671067, 0.000010, 0.001088]], "n" : ["||61435-4122"]},
			{"p" : [[40, 2.308570, 45.363721, 0.000010, 0.008363]], "n" : ["||61435-4949"]},
			{"p" : [[41, 1.938949, 37.197833, 0.000010, 0.149723]], "n" : ["||61435-317"]},
			{"p" : [[40, 0.043721, 38.058347, 0.000010, 0.008363]], "n" : ["||61435-2985"]},

```

You can proceed with your down-stream analysis using other tools, such as [`gappa`](github.com/lczech/gappa):
```bash
grep -v "\[ \]" placements_toy.jplace | sed -z "s/},\n\t]/}\n\t]/g" > placements_toy-filtered.jplace # filtering sequences without any match and hence an empty placement field
gappa examine heat-tree --jplace-path placements_toy-filtered.jplace --write-svg-tree # generating a colored tree based on placement densities across the backbone tree
```
