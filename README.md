# krepp
krepp is a k-mer-based maximum likelihood tool for estimating distances of reads to genomes and phylogenetic placement.

For the descripton of the method, please refer to the preprint [here](https://www.biorxiv.org/content/10.1101/2025.01.20.633730v1).

See the [Wiki](https://github.com/bo1929/krepp/wiki) for a detailed documentation, a list of available databases, and various tutorials.

## Installation
### Using `conda` (recommended)
The easiest way to install krepp is by using [conda](https://anaconda.org/bioconda/krepp)
```bash
conda install bioconda::krepp 
```
This will install the latest available version, simply run `krepp --help` to test.

### Compiling from the source
To compile from the source, clone the repository with its submodules (might take a while) and compile with
```bash
git clone --recurse-submodules -j8 https://github.com/bo1929/krepp.git
cd krepp && make
```
and run `./krepp --help`. Then, perhaps, copy it to a directory you have in your `$PATH` (e.g., `cp ./krepp ~/.local/bin`).


## Quickstart with a toy example
### Building a small index
You can build an index from scratch, consisting of only 25 genomes provided in `test/`, to make yourself familiar with `krepp`.
```bash
cd test/
tar -xvf references_toy.tar.gz && xz -d references_toy/*
krepp index -h 11 -k 27 -w 35 -o index_toy -i input_map.tsv -t tree_toy.nwk --num-threads 8
```
This command took only a couple of seconds and used $<$1.5GB memory for 6,863,411 indexed *k*-mers.
The resulting index will be stored in `index_toy`.
Alternatively, you could download one of the larger [public libraries](Available-reference-indexes) to make it more realistic and use it also for your novel query sequences.

### Querying sequences against the reference index
Once you have your index (e.g., the one we built above: `index_toy`), you can estimate distance by running:
```bash
krepp dist -i index_toy -q query_toy.fq --num-threads 4 | tee distances_toy.tsv
```
The first five lines of `distances_toy.tsv` are going to look like:
```
#software: krepp	#version: v0.6.0	#invocation :krepp dist -i index_toy -q query_toy.fq --num-threads 4
SEQ_ID	REFERENCE_NAME	DIST
||61435-4122	G000341695	0.0898062
||61435-4949	G000830905	0.147048
||61435-4949	G000341695	0.0740587
||61435-4949	G000025025	0.131182
||61435-4949	G000741845	0.0395985
```

Quite similarly, you can place reads by running:
```bash
krepp --num-threads 8 place -i index_toy -q query_toy.fq | tee placements_toy.jplace
```

The resulting placement file is a JSON file in a special format called [jplace](https://doi.org/10.1371/journal.pone.0031009):
```bash
head -n20 placements_toy.jplace
```
```
{
        "version" : 3,
        "fields" : ["edge_num", "pendant_length", "distal_length", "likelihood", "like_weight_ratio", "distance"],
        "placements" : [
                        {"n" : ["||61435-4122"], "p" : [[39, 0.0986, 0.0011, -18.9251, 1.0000, 0.0933]]},
                        {"n" : ["||61435-4949"], "p" : [
                                [38, 0.0433, 0.0007, -31.8820, 0.2240, 0.0427],
                                [37, 0.0523, 0.0000, -38.5286, 0.1582, 0.0505],
                                [40, 0.0400, 0.0084, -35.2375, 0.1900, 0.0468],
                                [36, 0.0314, 0.0020, -24.9208, 0.2696, 0.0326],
                                [39, 0.0512, 0.0011, -38.5314, 0.1582, 0.0505]]
                        },
                        {"n" : ["||61435-317"], "p" : [
                                [38, 0.0058, 0.0007, -18.0543, 0.1837, 0.0065],
                                [37, 0.0060, 0.0000, -18.0878, 0.1844, 0.0060],
                                [40, -0.0021, 0.0084, -18.0759, 0.1842, 0.0062],
                                [36, 0.0052, 0.0020, -18.0116, 0.1812, 0.0071],
                                [39, 0.0049, 0.0011, -18.0924, 0.1844, 0.0060],
                                [41, -0.1333, 0.1497, -22.2532, 0.0821, 0.0162]]
                        },
```
Here, `n` field is for the read ID as it appeared in `query_toy.fq`, and `p` is for the placement information, with the following fields:
```
["edge_num", "pendant_length", "distal_length", "likelihood", "like_weight_ratio", "distance"]
```
At the end of the jplace file, you can find the phylogeny decorated with edge numbers, which corresponds to the first field of `p`.

You can proceed with your downstream analysis using other tools, such as [`gappa`](https://github.com/lczech/gappa). e.g., by generating a heat tree, colored based on placement densities across the backbone tree:
```bash
gappa examine heat-tree --jplace-path placements_toy.jplace --write-svg-tree
```

Alternatively, for a simpler format, give the `--tabular` flag. Then, the first 20 lines would look like:
```
# software: krepp       version: v0.6.0 invocation :krepp --num-threads 8 place -i index_toy -q query_toy.fq --tabular
# (G001917855:0.4290{0},(G001918235:0.5280{1},(((((G000526415:0.1764{2},G001306135:0.2276{3})N1779:0.0365{4},(G000016665:0.1683{5},(G000735195:0.0276{6},G000018865:0.0229{7})N2640:0.1610{8})N1780:0.0609{9})N1532:0.1058{10},G000021685:0.3879{11})N1303:0.0337{12},(G002010545:0.3919{13},((G001050235:0.1711{14},G001306055:0.1635{15})N5461:0.1783{16},G001567105:0.2932{17})N2355:0.1933{18})N1304:0.0560{19})N1099:0.0288{20},(G000702505:0.1445{21},G001914715:0.1983{22})N1305:0.2353{23},(G001796575:0.3977{24},G001795015:0.4213{25},((G001796415:0.2756{26},(G002010445:0.3607{27},G001795205:0.3152{28})N3984:0.0458{29})N3292:0.0233{30},(((G000025025:0.0115{31},G001889305:0.0109{32})N4334:0.0068{33},G000830905:0.0305{34})N3987:0.0117{35},((G000741845:0.0039{36},G001610775:0.0001{37})N5905:0.0014{38},G000341695:0.0022{39})N4337:0.0167{40})N3634:0.2994{41})N2954:0.1317{42})N1788:0.1622{43})N916:0.0295{44})N736:0.0348{45})N432:0.0257{46};
SEQ_ID  DISTAL_NODE     EDGE_NUM        LWR     DIST
||61435-4122    G000341695      39      1.0000  0.0933
||61435-4949    G000741845      36      0.2696  0.0326
||61435-4949    N5905   38      0.2240  0.0427
||61435-4949    G000341695      39      0.1582  0.0505
||61435-4949    G001610775      37      0.1582  0.0505
||61435-4949    N4337   40      0.1900  0.0468
||61435-317     N3634   41      0.0821  0.0162
||61435-317     G000741845      36      0.1812  0.0071
||61435-317     N5905   38      0.1837  0.0065
||61435-317     G000341695      39      0.1844  0.0060
||61435-317     G001610775      37      0.1844  0.0060
||61435-317     N4337   40      0.1842  0.0062
||61435-2985    G000741845      36      0.2048  0.0158
||61435-2985    N5905   38      0.2023  0.0175
||61435-2985    G000341695      39      0.1966  0.0189
||61435-2985    G001610775      37      0.1966  0.0189
||61435-2985    N4337   40      0.1997  0.0183
```

## Citation
```bibtex
@misc{sapci_k-mer-based_2025,
	title = {A k-mer-based maximum likelihood method for estimating distances of reads to genomes enables genome-wide phylogenetic placement.},
	copyright = {2025, Posted by Cold Spring Harbor Laboratory. The copyright holder for this pre-print is the author. All rights reserved. The material may not be redistributed, re-used or adapted without the author's permission.},
	url = {https://www.biorxiv.org/content/10.1101/2025.01.20.633730v1},
	doi = {10.1101/2025.01.20.633730},
	language = {en},
	urldate = {2025-01-27},
	publisher = {bioRxiv},
	author = {Sapci, Ali Osman Berk and Mirarab, Siavash},
	month = jan,
	year = {2025},
	note = {Pages: 2025.01.20.633730
Section: New Results},
}

```
