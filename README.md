<p align="center"><img src="media/logo.svg" alt="BKTyper - v.0.1" width="450"></p>

BKTyper is a human polyoma BK (BKPyV) typing tool. It can automatically classify transcriptional blocks from the Non-Coding Control Region (NCCR) and the Viral Protein 1 (VP1). It can take a nucleotide sequence in fasta format containing the complete BKPyV genome or subgenomic region containing the VP1 gene or the NCCR region.

## Requirments

BKTyper depends on 4 different external programs that should be installed in your operational system and be on path:

- [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [EMBOSS toolkit](ftp://emboss.open-bio.org/pub/EMBOSS/)
- [IQTree](http://www.iqtree.org/)

## Installation

BKTyper is written in Python3 and depends on several Python libraries. Clone the repository and install the requiered Python libraries:

```
git clone https://github.com/joanmarticarreras/BKTyper.git
cd BKTyper/
pip3 install -r requirements.txt
```

## Run

BKTyper requires a fasta file as input, which should contain a BKPyV nucleotide sequence (full genome, VP1 or NCCR regions). Additionally it requires as extra arguments the type of sequence is being input (complete for the complete genome, VP1 or NCCR as what the sequence contains). Finally, it is possible to have phylogenetic placing throuigh a Maxmimum-Likelihood tree with IQTree, type "Tree" or "noTree" as to enable the options (for "complete" or "VP1" mode only).

`python3 BKTyper.py <input-sequences> <mode: VP1, NCCR, complete> <phylogeny: Tree noTree>`
