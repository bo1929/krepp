# krepp
krepp is a *k*-mer-based maximum likelihood tool for estimating distances of reads to genomes and phylogenetic placement.

For the description of the method, refer to the main paper [here](https://link.springer.com/article/10.1186/s13059-026-03999-y).

See the [Wiki](https://github.com/bo1929/krepp/wiki) for detailed documentation, a list of available databases, and various tutorials.

See [this tutorial for metagenomic analysis of a toy marine sample](https://alisapci.xyz/documents/krepp-tutorial-ismb2026/index.html) (from ISMB 2026).

## Installation
### Using `conda` (recommended)
The easiest way to install krepp is by using [conda](https://anaconda.org/bioconda/krepp)
```bash
conda install bioconda::krepp
```
This will install the latest available version. Simply run `krepp --help` to test.

### Compiling from the source
To compile from the source, clone the repository with its submodules and compile with
```bash
git clone --recurse-submodules -j8 https://github.com/bo1929/krepp.git
cd krepp && make
```
and run `./krepp --help`. Then, perhaps, install it to a directory you have in your `$PATH` with `sudo make install` (or `make install PREFIX=~/.local`, or simply `cp ./krepp ~/.local/bin`).


## Quickstart with a toy example
### Building a small index
You can build an index from scratch, using only the 25 genomes provided in `test/`, to familiarize yourself with `krepp`.
```bash
cd test/
tar -xvf references_toy.tar.gz && xz -d references_toy/*
krepp index -h 11 -k 27 -w 35 -o index_toy -i input_map.tsv -t tree_toy.nwk --num-threads 8
```
This command took only a couple of seconds and used ~1.5GB of memory to index 6,932,567 *k*-mers.
The resulting index will be stored in `index_toy`.
Alternatively, you could download one of the larger [public libraries](https://github.com/bo1929/krepp/wiki/Available-reference-indexes) to make it more realistic and use it also for your novel query sequences.

### Querying sequences against the reference index
Once you have your index (e.g., the one we built above: `index_toy`), you can estimate distance by running:
```bash
krepp dist -i index_toy -q query_toy.fq --num-threads 4 | tee distances_toy.tsv
```
The first few lines of `distances_toy.tsv` are going to look like:
```
# software: krepp	version: v0.11.1	invocation :krepp dist -i index_toy -q query_toy.fq --num-threads 4
SEQ_ID	REFERENCE_NAME	DIST
||61435-4122	G000341695	0.04617
||61435-4949	G000741845	0.02924
||61435-4949	G000341695	0.03017
||61435-4949	G001610775	0.03022
```
One row is reported per query and reference, up to the `--dist-max` cutoff.
The `dist` output is byte-for-byte reproducible; with `--num-threads` above 1, only the
order in which queries appear can vary between runs.
Add `--p-value` to append a column that says how distinguishable a hit is from the best hit of the same query:
```
SEQ_ID	REFERENCE_NAME	DIST	P_VALUE
||61435-4122	G000341695	0.04617	0.00000e+00
||61435-4949	G000741845	0.02924	0.00000e+00
||61435-4949	G000341695	0.03017	6.23556e-02
||61435-4949	G001610775	0.03022	6.61218e-02
```
It is 0 for anything tied with the best hit (there is no evidence against it) and
grows towards 1 as the hit becomes clearly worse, so smaller is better.
The option `--filter` keeps the hits whose statistic is below `--chisq`, i.e.
the ones with *p*-value is below 0.9 for the default `--chisq 2.706`.

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
			{"n" : ["||61435-4122"], "p" : [[39, 0.04738, 0.00027, -10.94198, 1.00000, 0.04617]]},
			{"n" : ["||61435-4949"], "p" : [
				[36, 0.02872, 0.00193, -34.17883, 0.18834, 0.02924],
				[41, 0.02982, 0.10049, -17.82452, 0.06163, 0.04962],
				[37, 0.03082, 0.00003, -20.30400, 0.18769, 0.03022],
				[39, 0.03035, 0.00108, -20.28802, 0.18776, 0.03017],
				[38, 0.03073, 0.00070, -27.47240, 0.18785, 0.03009],
				[40, 0.04097, 0.00635, -24.11682, 0.18673, 0.03080]]
			},
			{"n" : ["||61435-317"], "p" : [
				[36, 0.01090, 0.00136, -18.07846, 0.30264, 0.00713],
				[37, 0.01858, 0.00004, -26.26133, 0.14447, 0.01839],
				[39, 0.01744, 0.00116, -26.24350, 0.14479, 0.01836],
				[38, 0.01563, 0.00064, -22.66043, 0.22905, 0.01334],
				[40, 0.01925, 0.00590, -24.76044, 0.17905, 0.01623]]
			},
```
Here, the `n` field is for the read ID as it appeared in `query_toy.fq`, and `p` is for the placement information, with the following fields:
```
["edge_num", "pendant_length", "distal_length", "likelihood", "like_weight_ratio", "distance"]
```
At the end of the jplace file, you can find the phylogeny decorated with edge numbers, which corresponds to the first field of `p`.

You can proceed with your downstream analysis using other tools, such as [`gappa`](https://github.com/lczech/gappa), e.g., by generating a heat tree, colored based on placement densities across the backbone tree:
```bash
gappa examine heat-tree --jplace-path placements_toy.jplace --write-svg-tree
```

Alternatively, for a simpler format, give the `--tabular` flag. Then, the first 20 lines would look like:
```
# software: krepp	version: v0.11.1	invocation :krepp --num-threads 8 place -i index_toy -q query_toy.fq --tabular
# (G001917855:0.42900{0},(G001918235:0.52803{1},(((((G000526415:0.17636{2},G001306135:0.22758{3})N1779:0.03648{4},(G000016665:0.16835{5},(G000735195:0.02760{6},G000018865:0.02290{7})N2640:0.16098{8})N1780:0.06088{9})N1532:0.10584{10},G000021685:0.38795{11})N1303:0.03371{12},(G002010545:0.39187{13},((G001050235:0.17105{14},G001306055:0.16351{15})N5461:0.17835{16},G001567105:0.29316{17})N2355:0.19330{18})N1304:0.05605{19})N1099:0.02879{20},(G000702505:0.14447{21},G001914715:0.19828{22})N1305:0.23530{23},(G001796575:0.39771{24},G001795015:0.42134{25},((G001796415:0.27559{26},(G002010445:0.36066{27},G001795205:0.31521{28})N3984:0.04579{29})N3292:0.02335{30},(((G000025025:0.01149{31},G001889305:0.01088{32})N4334:0.00679{33},G000830905:0.03055{34})N3987:0.01170{35},((G000741845:0.00392{36},G001610775:0.00006{37})N5905:0.00142{38},G000341695:0.00218{39})N4337:0.01673{40})N3634:0.29945{41})N2954:0.13171{42})N1788:0.16221{43})N916:0.02948{44})N736:0.03481{45})N432:0.02566{46};
SEQ_ID	DISTAL_NODE	EDGE_NUM	LWR	DIST
||61435-4122	G000341695	39	1.00000	0.04617
||61435-4949	N5905	38	0.18785	0.03009
||61435-4949	G000341695	39	0.18776	0.03017
||61435-4949	N3634	41	0.06163	0.04962
||61435-4949	G000741845	36	0.18834	0.02924
||61435-4949	G001610775	37	0.18769	0.03022
||61435-4949	N4337	40	0.18673	0.03080
||61435-317	N5905	38	0.22905	0.01334
||61435-317	G000341695	39	0.14479	0.01836
||61435-317	G001610775	37	0.14447	0.01839
||61435-317	G000741845	36	0.30264	0.00713
||61435-317	N4337	40	0.17905	0.01623
||61435-2985	N5905	38	0.18785	0.00881
||61435-2985	G000341695	39	0.19048	0.00803
||61435-2985	N3634	41	0.06517	0.01653
||61435-2985	G001610775	37	0.19048	0.00804
||61435-2985	G000741845	36	0.17634	0.00995
```
Placement values are deterministic up to the last printed digit of `like_weight_ratio`
(a ratio of sums accumulated in hash-map order); with `--num-threads` above 1, the order
of queries in the file can also vary between runs.

## Python bindings

`make python` builds an extension module (using `pybind11`) into `python/krepp/`.
Distances, placements and per-query *k*-mer match summaries come back as Python objects

```python
import krepp                          # with PYTHONPATH=python (repo root)

s = krepp.Session("index_toy", krepp.Options(num_threads=4, p_value=True))

for qd in s.dist("query_toy.fq"):     # like `krepp dist`
    for row in qd.rows:
        print(qd.query, row.ref_name, row.dist, row.p_value)

for qp in s.place("query_toy.fq"):    # like `krepp place`
    for p in qp.placements:
        print(qp.query, p.edge_num, p.distal_node, p.like_weight_ratio)

for qm in s.matches("query_toy.fq"):  # k-mer match summaries per reference
    print(qm.query, qm.onmers, qm.enmers, qm.hdist_hist)   # aggregate histogram by Hamming distance
    for r in qm.rows:
        print("   ", r.ref_name, r.match_count, r.hdist_min, r.hdisthist)  # per reference
```

Use `make python PYTHON=/path/to/python3` to target a specific interpreter
(default: `python3`).

## Citation
```bibtex
@article{sapci_krepp_2026,
	title = {krepp: a k-mer-based maximum pseudo-likelihood method for estimating read distances and genome-wide phylogenetic placement},
	volume = {27},
	issn = {1474-760X},
	shorttitle = {krepp},
	url = {https://doi.org/10.1186/s13059-026-03999-y},
	doi = {10.1186/s13059-026-03999-y},
	number = {1},
	urldate = {2026-03-28},
	journal = {Genome Biology},
	author = {Şapcı, Ali Osman Berk and Mirarab, Siavash},
	month = feb,
	year = {2026},
	keywords = {Average nucleotide identity, k-mer-based sequence comparison, Metagenomics, Phylogenetic placement},
	pages = {108},
}
```
