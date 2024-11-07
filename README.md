# krepp
A k-mer-based maximum likelihood method for estimating distances of reads to genomes and phylogenetic placement.

## Quickstart
### Installing and compiling krepp
Pre-compiled binaries are not available yet, but hopefully will be soon. Clone the repository and compile with
```
git clone --recurse-submodules -j8 https://github.com/bo1929/krepp.git
cd krepp && make
```
and run `./krepp --help`. Perhaps copy it to a directory you have in your `$PATH`: `cp ./krepp ~/.local/bin`.

### Building a krepp index from reference genomes
Given a set of reference genomes and a backbone tree, krepp can build an LSH index of colored k-mers by simply running
```bash
krepp --num-threads $NUM_THREADS build -l $LIBRARY_DIR -i $PATHS_MAPPING -t $BACKBONE_NEWICK
```
where
* `$NUM_THREADS` is simply the number of threads (note that memory use increase as more threads used),
* `-l` is the directory where the index will be stored (preferably a unique name),
* `-i` is a tab-separated file where the first column is IDs of references, i.e., leaves of the tree, and the second column is paths to genomes (e.g., FASTA),
* `-t` is the backbone tree in Newick format matching the genome IDs given.

Run `krepp build --help` to see more options.
If you don't have a tree and mostly interested in estimating distances, an arbitrary (or random) tree might just work as a workaround (with the side effect of making the color index potentially a little larger).

### Estimating distances of reads to reference genomes
Once  you have the index built, query reads against it to get distance estimates is quite simple:
```bash
krepp dist -l $LIBRARY_DIR -q $QUERY_PATH
```
where `-q` is the path of a FASTA/Q file containing query reads, and `krepp` simply outputs everything to stdout and write log messages to stderr. The output is in a tab-separated format in which the first column stands for the read ID, second column is the ID of the reference matching, and the third column is the distance estimate of krepp.

### Phylogenetic placement of reads on the backbone tree
In addition to distance estimation, one could place reads on the backbone tree given as an input while building the index
```bash
krepp place -l $LIBRARY_DIR -q $QUERY_PATH
```
where output is again in a tab-separated format and each row corresponds to a read.
We will make the standard placement format `jplace` available very soon with the next release.
