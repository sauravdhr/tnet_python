# TNet
TNet: Phylogeny-Based Inference of Disease Transmission Networks UsingWithin-Host Strain Diversity

Input: TNet takes a rooted phylogeny as input. The file should be in [Newick](https://en.wikipedia.org/wiki/Newick_format) format. To get the phylogeny from mapped reads we used [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html). We used RAxML Version 8 for the current analysis. The phylogeny has to be rooted. We used [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) again for rooting the phylogeny. Also make sure the input tree is bifurcating. The name format for the leaves is `<hostID>_<sequenceID>`. Only `<hostID>` will be used in the tool for computation and the output will be based on the `<hostID>` too.

Output: The output file is a file specified by the user. Each line of the output file is an edge in the transmission network. The edge format is `<hostID> <hostID>` separated by a tab. The first line contains the source of the outbreak in `None <hostID>` format. TNet also outputs the minimum parsimony cost of the phylogeny.

Running TNet: TNet is implemented in [Python3](https://www.python.org/download/releases/3.0/). You can just use the `tnet.py` file to run TNet from command line.
```bash
./tnet.py input_file output_file
```
or
```bash
python3 tnet.py input_file output_file
```

### Example
The `input` directory on the git repository contains some rooted phylogeny files. There is a `output` directory for the TNet outputs. Use the following command line to run TNet on the test files.

```bash
./tnet.py input/RAxML_rootedTree.1 output/1.out
```
