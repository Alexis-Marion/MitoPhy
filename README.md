# MitoPhy

## Contents

- [Contents](#Contents)
- [1 Introduction](#1-Introduction)
- [2 Installation](#2-Installation)
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
- [Reference](#Reference)

## 1 Introduction

MitoPhy is a visual interface written in python 3 for Unix computers dedicated to sequence management and molecular phylogenetics. The purpose of MitoPhy is to provide a quick and efficient way to build phylogenies de novo. I created the general idea for MitoPhy during my Master's internship in 2022. I was asked to build a large phylogeny using mainly mitochondrial genes for more than 350 species. No sequences were already available in my lab, so I relied on using NCBI's GENBANK. However the task was a bit too overwhelming for me and before this internship, I had no genuine idea of how Genbank works. I asked myself several questions such as: what is a "good" sequence: is it a long one, a well-attributed one, or a contamination-free one? To make it simple, all of the above. But since my time was limited, I could not simply examine each possible gene for each species more than 5200 times and look for these criteria. Moreover, with such a large number, human errors can happen easily. So I decided to build a pipeline with several programs and scripts, making this task easier. Six months after and the first version is complete, with numerous options included and an in-built visual interface. 

Many programs used in MitoPhy were not developed by me, so please cite each of them when used. If you have any questions/requests feel free to ask me at my Github.

## 2 Installation

MitoPhy is a standalone program, as such, no prerequisites are needed.


## Download manager

In this chapter, we will address the different download options for sequences, using the E-utilities packages in python.
### 3.1 Sequence selection

#### 3.1.1 Write taxa name

This panel allows the user to enter a taxa name (at the moment the binomial species name, ie: Genus species format). This panel is quite straightforward as you have only to write the species name and hit the
button confirm next to hit. Keep in mind that this method is useful when you have only a few species as you will have to hit the search button each time you change the species name to download it. This and "Select taxa list" are mutually exclusive, either way, one of these options is mandatory for the research.

#### 3.1.2 Select taxa list

This panel allows the user to enter a .txt file as input. This file can have more than one species, each written on a different line. This download method is far superior to the previous, especially with numerous species. This and "Write taxa name" are mutually exclusive, in both cases, one of these options is mandatory for the research option.

### 3.2 Save, load, and convert

#### 3.2.1 Save in directory

This panel allows the user to enter an output directory for the downloaded file. This step is mandatory for
the research option.

#### 3.2.2 Convert alignment file

This panel allows the user to enter any sequence file and convert it from one of the following formats to another. This step is optional.
• clustal: Format used by the clustal X alignment tool, files must start with CLUSTAL X or CLUSTALX. The input and output files must and will end with .aln.
• fastq: Format including nucleotide sequence along with a quality score. The input and output files must and will end with .fastq or .fq.
• fasta: Standard format used by many programs. The input file must and will end with .aln.
• Genbank: Format used by NCBI's Genbank, including sequences along with detailed descriptions The input and output files must and will end with .genbank or .gb.
• nexus: Format used by many programs, such as mrbayes or BEAST. The input and output files must and will end with .nex.
• phylip: Format used by many programs, such as IQ-TREE or RAxML. The input and output files must and will end with .phylip.

### 3.3 Sequence option

#### 3.3.1 Only RefSeq sequences

This option allows the user to be looking only for RefSeq sequences. Refseq sequences are of high quality but are less represented than any non-RefSeq sequence, usually, each RefSeq sequence contains the whole
mitogenome. This step is optional.

#### 3.3.2 Date range

This option allows the user to enter a date interval for the downloaded sequence. As such if one writes 2014/08/03 and 2022/03/01, the downloaded sequence cannot be older than the third of august 2014 but
are only younger than the first of march 2022. The format must be YYYY/MM/DD and both of the panels must be completed in a valid format. This step is optional, keep the panel blank if you do not want any date
constraints.

#### 3.3.3 Sequence range

This option allows the user to enter a sequence length interval in which the downloaded sequence cannot be smaller or bigger. Same as "Date range”, you must complete both of the panels if you want a constraint. This step is optional, keep the panel blank if you do not want any sequence length constraint.

### 3.4 Research bar

At the moment, the research bar submenu only contains one button: the research button. Taking into account all of the previous options and input, the research button will perform research on NCBI's Genbank database. For each species provided, a folder will be created along with a .gb (GenBank format) containing many sequences for the species (capped at 100 to not overload the NCBI's servers).

## 4 Sequence manipulation

In this chapter, I will show how to convert and align the sequence obtained in chapter 3 and other sequence-related operations.

### 4.1 Sequence conversion


#### 4.1.1 Select the operating directory

This panel allows the user to enter an input directory for the following operation. The structure of this directory is with two levels, the first is the species level, and the second is inside the first one. If none is written here, the previous from section 3.2.1 will be entered. Either way, this step stays mandatory.

#### 4.1.2 Select the output directory

Similar to 3.2.1, the user can enter a path to an output directory. This step is mandatory.

#### 4.1.3 Turn into fasta

In this panel, all .gb files in the operating directory will be transformed into .fasta files and transferred to new subdirectories for each gene. The read direction is 5’ - 3’ (+). This step is mandatory for step 4.1.4.

#### 4.1.4 Consensus

In this panel, in each subdirectory defined in 4.1.3, all sequences will be transferred to the new output directory by genes. Here three options are present :
• 1 You have no sequences in the directory, as such a placeholder sequence, containing only "-" is created and moved.
• 2 You have only one sequence, no consensus is made and the sequence is moved to the output
directory unchanged.
• 3 You have more than two sequences and a majority of the sequence is made and moved to the output directory. Ambiguous sites are "X”.

### 4.2 Alignment parameters

#### 4.2.1 Select the operating directory

Similar to 3.1.1, this panel allows the user to select an operating directory for the rest of his data manipulation. The structure here consists of the main directory with a gene-related subdirectory. This step is mandatory.

#### 4.2.2 Align

With the directory selected in 5.1.1, align all the files in each subdirectory with muscle and create output in the main root directory.

#### 4.2.3 Merge

Merge all the .fasta files from the main directory by species names. Also creates two partition files in .txt format, which could be used in the last chapter.

#### 4.2.4 Alignement info


## 5 Phylogeny

### 5.1 General options

#### 5.1.1 Select input sequence

#### 5.1.2 Select outgroup

#### 5.1.3 Select bootstrap value

#### 5.1.4 Select a constraint tree

#### 5.1.5 Save in directory

### 5.2 Distance and parsimony

#### 5.2.1 Select method

• UPGMA
• Neighbour Joining (NJ)
• Parsimony
#### 5.2.2 Select output format

• Newick
• nexus
• nexml
• phyloxml

#### 5.2.3 Select consensus method

• Strict consensus methods only show relationships that are unambiguously supported by the data for all the trees
• Majority consensus methods only show relationships that are unambiguously supported by the data the majority of the trees
• Adams consensus methods, as described in Adams (1972)

#### 5.2.4 Build distance or parsimony tree

This panel allows the user to build a distance-based tree with maximum parsimony using all the options in sections 5.1 and 5.2.

### 5.3 Maximum likelihood

#### 5.3.1 Select substitution model

Select a substitution model for your sequence, this option is overwritten when using partition files.
• Jukes and Cantor 69 (JC69)
• Felsenstein 81 (F81)
• Kimura two parameters (K2P)
• Hasegana-Kishino-Yana 85 (HKY85)
• Tamura Nei 93 (TN93)
• Kimura three parameters (K3P)
• Transition model (TIM)
• Transversion model (TVM)
• Symmetric model (SYM)
• General time revesible (GTR)
• Model finder (MPF)

#### 5.3.2 Use gamma distribution

#### 5.3.3 Use invariable sites

#### 5.3.4 Select the partition file

#### 5.3.5 Select bootstrap type

• Bootstrap
• Ultrafastbootstrap
#### 5.3.6 Build maximum likelihood tree using IQTREE

# Reference
