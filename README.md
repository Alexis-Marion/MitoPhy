# MitoPhy 1.0
![banner.png]("https://github.com/Alexis-Marion/MitoPhy/blob/main/banner.png")
## Contents

- [Contents](#Contents)
- [1 Introduction](#1-Introduction)
- [2 Installation](#2-Installation)
	- [2.1 MitoPhy for Mac](#21-MitoPhy-for-Mac)
	- [2.2 MitoPhy for Linux](#22-MitoPhy-for-Linux)
- [3 Download manager](#3-Download-manager)
	- [3.1 Sequence selection](#31-Sequence-selection)
    - [3.1.1 Write taxa name](#311-Write-taxa-name)
    - [3.1.2 Select taxa list](#312-Select-taxa-list)
	- [3.2 Save, load, and convert](#32-Save,-load,-and-convert)
		- [3.2.1 Save in directory](#321-Save-in-directory)
		- [3.2.2 Convert alignment file](#322-Convert-alignment-file)
	- [3.3 Sequence option](#33-Sequence-option)
		- [3.3.1 Only RefSeq sequences](#331-Only-RefSeq-sequences)
		- [3.3.2 Date range](#332-Date-range) 
    		- [3.3.3 Sequence range](#333-Sequence-range)
  - [3.4 Research bar](#34-Research-bar)
- [4 Sequence manipulation](#4-Sequence-manipulation)
	- [4.1 Sequence conversion](#41-Sequence-conversion)
    	 - [4.1.1 Select the operating directory](#411-Select-the-operating-directory)
    	 - [4.1.2 Select the output directory](#412-Select-the-output-directory)
    	 - [4.1.3 Turn into fasta](#413-Turn-into-fasta)
    	 - [4.1.4 Consensus](#414-Consensus)
	- [4.2 Alignment parameters](#42-Alignment-parameters)
    	 - [4.2.1 Select the operating directory](#421-Select-the-operating-directory)
    	 - [4.2.2 Align](#422-Align)
    	 - [4.2.3 Merge](#423-Merge)
    	 - [4.2.4 Alignement info](#424-Alignement-info)
- [5 Phylogeny](#5-Phylogeny)
	- [5.1 General options](#51-General-options)
		- [5.1.1 Select input sequence](#511-Select-input-sequence)
		- [5.1.2 Select outgroup](#512-Select-outgroup)
		- [5.1.3 Select bootstrap value](#513-Select-bootstrap-value)
		- [5.1.4 Select a constraint tree](#514-Select-a-constraint-tree)
		- [5.1.5 Save in directory](#515-Save-in-directory)
  	- [5.2 Distance and parsimony](#52-Distance-and-parsimony)
		- [5.2.1 Select method](#521-Select-method)
		- [5.2.2 Select output format](#522-Select-output-format)
		- [5.2.3 Select consensus method](#523-Select-consensus-method)
		- [5.2.4 Build distance or parsimony tree](#524-Build-distance-or-parsimony-tree)
  	- [5.3 Maximum likelihood](#53-Maximum-likelihood)
		- [5.3.1 Select substitution model](#531-Select-substitution-model)
		- [5.3.2 Use gamma distribution](#532-Use-gamma-distribution)
		- [5.3.3 Use invariable sites](#533-Use-invariable-sites)
		- [5.3.4 Select the partition file](#534-Select-the-partition-file)
		- [5.3.5 Select bootstrap type](#535-Select-bootstrap-type)
		- [5.3.6 Build maximum likelihood tree using IQTREE](#536-Build-maximum-likelihood-tree-using-IQTREE)
- [6 Known issues](#6-Known-issues)
- [7 Future projects](#7-Future-projects)
- [Reference](#Reference)

## 1 Introduction

<p align="justify"> MitoPhy is a visual interface written in python 3 for Unix computers dedicated to sequence management and molecular phylogenetics. The purpose of MitoPhy is to provide a quick and efficient way to build phylogenies de novo. I created the general idea for MitoPhy during my Master's internship in 2022. I was asked to build a large phylogeny using mainly mitochondrial genes for more than 350 species. No sequences were already available in my lab, so I relied on using NCBI's GENBANK. However the task was a bit too overwhelming for me and before this internship, I had no genuine idea of how Genbank works. I asked myself several questions such as: what is a "good" sequence: is it a long one, a well-attributed one, or a contamination-free one? To make it simple, all of the above. But since my time was limited, I could not simply examine each possible gene for each species more than 5200 times and look for these criteria. Moreover, with such a large number, human errors can happen easily. So I decided to build a pipeline with several programs and scripts, making this task easier. Six months after and the first version is complete, with numerous options included and an in-built visual interface. </p>

<p align="justify"> Mitophy mainly rely on taking actions on directory recursively, meaning that when proposed to select a direcotry, the user should create a new one. Otherwise the built-in loop in the programm can go out of hand, and cause several issues. In doubt work only with pristine clear directories or kill the process in the terminal. </p>

<p align="justify"> Many programs used in MitoPhy were not developed by me such as muscle (Edgar, 2004), IQ-TREE (Minh et al., 2020) as well as python pacakges : module Bio (Cock et al., 2009), pandas (McKinney, W. et al., 2010), so please cite each of them when using the associated function. If you have any questions/requests feel free to ask me at my Github. </p>

## 2 Installation

MitoPhy is a standalone program, as such, no prerequisites are needed.
Mac and Linux (Ubuntu/Debian) are distributed and can be respectively be obtained at my figshare repository or here.

### 2.1 MitoPhy for Mac

MitoPhy for Mac OS can be obtained at "". MitoPhy.app is a bundled app that can be used like any other app from the App Store. However this distribution do not display the terminal, as such bugs are hardly tracked down. To counteract this, the user can enter manually inside the MitoPhy.app and locate the MitoPhy executable file (MitoPhy/Contents/MacOS/MitoPhy). This procedure is quite handly as it can help you to keep tracks on long procedure, such as phylogenetic reconstruction.

### 2.2 MitoPhy for Linux

MitoPhy for Linux can be downloaded right here. I recommend to install MitoPhy on the $HOME directory. (ADDING MITOPHY TO PATH)

```diff
- Mandatory step
```

## 3 Download manager

In this chapter, we will address the different download options for sequences, using the E-utilities packages in python.
### 3.1 Sequence selection

#### 3.1.1 Write taxa name

```diff
- Mandatory step
```

<p align="justify"> This panel allows the user to enter a taxa name (at the moment the binomial species name, ie: Genus species format). This panel is quite straightforward as you have only to write the species name and hit the button confirm next to hit. Keep in mind that while this method is useful for few species, you will have to hit the search button each time you change the species name to download it, which can be not pratical for large datasets. At the moment this panel will allow you to donwnload the sequences corresponding to the name entered. As such, in theory, any order/family/genus name will be considered as correct an be downloaded, but not separated, meaning that the user should always enter the lowest taxonomical rank for an organism. As an example, if one write "Canis lupus" (grey wolf), all subspecies will also be downloaded </p>

#### 3.1.2 Select taxa list

```diff
- Mandatory step
```

<p align="justify"> This panel allows the user to enter a .txt file as input. This file can have more than one species, each written on a different line. This download method is far superior to the previous, especially with numerous species. </p>

### 3.2 Save, load, and convert

#### 3.2.1 Save in directory

```diff
- Mandatory step
```

This panel allows the user to enter an output directory for the downloaded file. 

#### 3.2.2 Convert alignment file

This panel allows the user to enter any sequence file and convert it from one of the following formats to another. 

• **clustal**: Format used by the clustal X alignment tool, files must start with CLUSTAL X or CLUSTALX. The input and output files must and will end with .aln.

• **fastq**: Format including nucleotide sequence along with a quality score. Works only for input file, must end with .fastq or .fq.

• **fasta**: Standard format used by many programs. The input file must and will end with .aln.

• **Genbank**: Format used by NCBI's Genbank, including sequences along with detailed descriptions The input and output files must and will end with .genbank or .gb.

• **nexus**: Format used by many programs, such as mrbayes or BEAST. The input and output files must and will end with .nex.

• **phylip**: Format used by many programs, such as IQ-TREE or RAxML. The input and output files must and will end with .phylip.

### 3.3 Sequence option

#### 3.3.1 Only RefSeq sequences

This option allows the user to be looking only for RefSeq sequences. Refseq sequences are of high quality but are less represented than any non-RefSeq sequence, usually, each RefSeq sequence contains the whole mitogenome.

#### 3.3.2 Date range

<p align="justify"> This option allows the user to enter a date interval for the downloaded sequence. As such if one writes 2014/08/03 and 2022/03/01, the downloaded sequence cannot be older than the third of august 2014 but are only younger than the first of march 2022. The format must be YYYY/MM/DD and both of the panels must be completed in a valid format. This step is optional, keep the panel blank if you do not want any date constraints. </p>

#### 3.3.3 Sequence range

<p align="justify"> This option allows the user to enter a sequence length interval in which the downloaded sequence cannot be smaller or bigger. Same as "Date range”, you must complete both of the panels if you want a constraint. This step is optional, keep the panel blank if you do not want any sequence length constraint. </p>

### 3.4 Research bar

<p align="justify"> At the moment, the research bar submenu only contains one button: the research button. Taking into account all of the previous options and input, the research button will perform research on NCBI's Genbank database. For each species provided, a folder will be created along with a .gb (GenBank format) containing many sequences for the species (capped at 100 to not overload the NCBI's servers). </p>

## 4 Sequence manipulation

In this chapter, I will show how to convert and align the sequence obtained in [chapter 3](#3-Download-manager) and other sequence-related operations.

### 4.1 Sequence conversion

#### 4.1.1 Select the operating directory

```diff
- Mandatory step
```

This panel allows the user to enter an input directory for the following operation. The structure of this directory is with two levels, the first is the species level, and the second is inside the first one. If none is written here, the previous from section [3.2.1](#321-Save-in-directory) will be entered.

#### 4.1.2 Select the output directory

```diff
- Mandatory step
```

Similar to [3.2.1](#321-Save-in-directory), the user can enter a path to an output directory.

#### 4.1.3 Turn into fasta

```diff
- Mandatory step
```

In this panel, all .gb files in the operating directory will be transformed into .fasta files and transferred to new subdirectories for each gene. The read direction is *5*’ - *3*’ (+).

#### 4.1.4 Consensus

```diff
- Mandatory step
```

In this panel, in each subdirectory defined in [4.1.3](#413-Turn-into-fasta), all sequences will be transferred to the new output directory by genes. Here three options are present :

• You have **no** sequences in the directory, as such a placeholder sequence, containing only "-" is created and moved.

• You have only **one** sequence, no consensus is made and the sequence is moved to the output directory unchanged.

• You have more than **two** sequences and a majority of the sequence is made and moved to the output directory. Ambiguous sites are "N”.

### 4.2 Alignment parameters

#### 4.2.1 Select the operating directory

```diff
- Mandatory step
```

Similar to [3.1.1](#311-Write-taxa-name), this panel allows the user to select an operating directory for the rest of his data manipulation. The structure here consists of the main directory with a gene-related subdirectory. 

#### 4.2.2 Align

```diff
- Mandatory step
```

With the directory selected in [4.1.1](#411-Select-the-operating-directory), align all the files in each subdirectory with muscle and create output in the main root directory.

#### 4.2.3 Merge

```diff
- Mandatory step
```

Merge all the .fasta files from the main directory by species names. Also creates two partition files in .txt format, which could be used in the [last chapter](#5-Phylogeny).

#### 4.2.4 Alignement info

Select a multi-sequence file and give information about : the number of sequence, the number of invariable site, variable sites and only indel sites. 

## 5 Phylogeny

### 5.1 General options

#### 5.1.1 Select input sequence

```diff
- Mandatory step
```
Select an input multi-sequence of valid format to build the phylogeny. 

#### 5.1.2 Select outgroup

```diff
- Mandatory step
```
Select an outgroup among the sequence referenced in your dataset.

#### 5.1.3 Select bootstrap value

Select a non-parametric bootstrap value. 

#### 5.1.4 Select a constraint tree

Select a tree in newick format to constraint certains nodes in the phylogenetic analysis.

#### 5.1.5 Save in directory

```diff
- Mandatory step
```

Similar to [3.2.1](#321-Save-in-directory), this panel allows the user to enter an output directory for the reconstructed tree and additional files. 

### 5.2 Distance and parsimony

All options used for distance-based and maximum parsiomny methods.

#### 5.2.1 Select method

```diff
- Mandatory step
```

Select the algortihm best suited for your data.

• **UPGMA** (Distance based algortihm, assumes a constant substitution rate, over time and phylogenetic lineages (known as the molecular clock hypothesis)

• **Neighbour Joining (NJ)** (Distance based algorithm, reconstructs phylogenetic trees from evolutionary distance data)

• **Parsimony** (Cladistic,  optimality criterion under which the phylogenetic tree that minimizes the total number of character-state changes)

#### 5.2.2 Select output format

Select output format for the tree, default is Newick.

• **Newick** (or New Hampshire tree format, way of representing trees with parentheses and commas, standard format for phylogenetic trees)

• **nexus** (tree format, starting with "#Nexus" and followed by a  "TREES block" with trees written in Newick format)

• **nexml** (similar to nexus format, but using XML for richer phylogenetic data)

• **phyloxml** (XML format for rich node and branch annotation)

#### 5.2.3 Select consensus method


• **Strict consensus** methods only show relationships that are unambiguously supported by the data for all the trees

• **Majority consensus** methods only show relationships that are unambiguously supported by the data the majority of the trees

• **Adams consensus** methods, as described in Adams (1972)

#### 5.2.4 Build distance or parsimony tree

```diff
- Mandatory step
```

This panel allows the user to build a distance-based tree with maximum parsimony using the options in sections [5.1](#51-General-options) and [5.2](#52-Distance-and-parsimony).

### 5.3 Maximum likelihood

All options used for building maximum likelihood trees using IQ-TREE. 

#### 5.3.1 Select substitution model

```diff
- Mandatory step
```

Select a substitution model for your sequence, this option is overwritten when using partition files.

• **Jukes and Cantor 69 (JC69)** (Equal substitution rates and base frequencies, + 0 degree of freedom, Jukes and Cantor, 1969)

• **Felsenstein 81 (F81)** (Equal substitution rates and unequal base frequencies, + 3 degree of freedom, Felsentein 1981)

• **Kimura two parameters (K2P)** (Unequal transition/transversion rates an equal base frequencies, + 1 degree of freedom, Kimura 1980)

• **Hasegana-Kishino-Yana 85 (HKY85)** (Unequal transition/transversion rates an unequal base frequencies, + 4 degree of freedom, Hasegawa, Kishino and Yano 1985)

• **Tamura Nei 93 (TN93)** (Similar to HKY85 but unequal purine/pyrimidine substitution rates, + 5  degree of freedom, Tamura and Nei, 1993)

• **Kimura three parameters (K3P)** (Two transversion and one transition rates and equal base frequencies, + 3  degree of freedom, Kimura 1981)

• **Transition model (TIM)** (Unequal AC and GT *vs* AT and CG substitution rates and unequal base frequencies, + 6 degree of freedom)

• **Transversion model (TVM)** (Unequal substitution rates except AG and CT and unequal base frequencies, + 6 degree of freedom)

• **Symmetric model (SYM)** (Unequal substitution rates and equal base frequencies, + 5 degree of freedom, Zharkikh 1994)

• **General time revesible (GTR)** (Unequal substitution rates and unequal base frequencies, + 8 degree of freedom, Tavare 1986)

• **Model finder (MPF)** (Algorithm designed to search for the likeliest model suited for the data, time-consuming, Kalyaanamoorthy et al., 2017)

#### 5.3.2 Use gamma distribution

Allows for a proportion of invariable sites.

#### 5.3.3 Use invariable sites

Use discrete Gamma model for rate heterogeneity accross sites with 4 categories (Yang, 1994)

#### 5.3.4 Select the partition file

Select a partition file for the analysis. Partition are written in order to fit multiple substitution model for each of them. Partition files can be written by hand, but are also autogenerated by [step 4.2.3](#423-Merge) and ends with .nex.

#### 5.3.5 Select bootstrap type

```diff
- Mandatory step
```

• **Bootstrap** Non-parametric bootstrap is a random sampling with replacement creating pseudo-matrices for building alternative trees to compare the originnal to. Node with >= 70 % support are robust. Final tree ends with *.contree*.

• **Ultrafastbootstrap** Ultra fast bootstrap approximation for 1000 replicates, only nodes with  >= 95 % support are robust (D.T. Hoang et al., 2018). Final tree ends with *.contree*.

#### 5.3.6 Build maximum likelihood tree using IQTREE

```diff
- Mandatory step
```
Build a maximum likelihood tree with IQ-TREE, several other files will also be created.

## 6 Known issues

# Reference

Cock, P. J. Antao, T. Chang, J. T. Chapman, B. A. Cox, C. J. Dalke, A. De Hoon, M. J. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 25(11):1422-1423.

D.T. Hoang, O. Chernomor, A. von Haeseler, B.Q. Minh, and L.S. Vinh (2018) UFBoot2: Improving the ultrafast bootstrap approximation. Mol. Biol. Evol., 35:518–522. https://doi.org/10.1093/molbev/msx281

Edgar, R.C. 2004. MUSCLE: a multiple sequence alignment method with reduced time and space complexity. *BMC Bioinformatics*. 5:113 https://doi.org/10.1186/1471-2105-5-113

Felsenstein, J. Evolutionary trees from DNA sequences: A maximum likelihood approach. J Mol Evol 17, 368–376 (1981). https://doi.org/10.1007/BF01734359

Hasegawa, M., Kishino, H. & Yano, Ta. Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. J Mol Evol 22, 160–174 (1985). https://doi.org/10.1007/BF02101694

Jukes, Thomas H., and Charles R. Cantor. "Evolution of Protein Molecules." Mammalian Protein Metabolism, (1969): 21-132. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7.

S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, L.S. Jermiin (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. Nat. Methods, 14:587-589. https://doi.org/10.1038/nmeth.4285

Kimura, M. A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. J Mol Evol 16, 111–120 (1980). https://doi.org/10.1007/BF01731581

Kimura M. Estimation of evolutionary distances between homologous nucleotide sequences. Proc Natl Acad Sci U S A. 1981 Jan;78(1):454-8. doi: 10.1073/pnas.78.1.454. PMID: 6165991; PMCID: PMC319072.

Minh, B. Q. Schmidt, H. A. Chernomor, O. Schrempf, D. Woodhams, M. D. von Haeseler, A. Lanfear, R. 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular biology and evolution*. 37(5):1530–1534. https://doi.org/10.1093/molbev/msaa015  

McKinney, W., 2010. Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference. pp. 51–56.

Tamura, K., Nei, M., Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees., Molecular Biology and Evolution, Volume 10, Issue 3, May 1993, Pages 512–526, https://doi.org/10.1093/oxfordjournals.molbev.a040023

Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. Lectures on mathematics in the life sciences, 17(2), 57-86.

Yang, Z. Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: Approximate methods. J Mol Evol 39, 306–314 (1994). https://doi.org/10.1007/BF00160154

Zharkikh, A. (1994). Estimation of evolutionary distances between nucleotide sequences. Journal of molecular evolution, 39(3), 315-329.
