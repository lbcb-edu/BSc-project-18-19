# BSc project (computer science - 2018/2019)

Software Design Project is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the fifth semester of the undergraduate study. The main focus is to promote cooperation between students while solving specific problems. Under the supervision of prof. Mile Šikić, students will get familiar with C++, basics of compilation methods, version control and continuous integration, and will be introduced to algorithms used in bioinformatics. This project will be executed in several steps, each with defined outcomes which are required for succeeding steps. Instructions and guidelines for the project will be written in this README file which will be updated through the semester.

## Preliminaries

Students are required to get through the following tutorials: [C++](http://www.cplusplus.com/doc/tutorial/), [GitHub](http://rogerdudler.github.io/git-guide/), [CMake](https://cmake.org/cmake-tutorial/) and [TravisCI](https://docs.travis-ci.com/user/getting-started/). While writing C++ code, it is advised to follow the [Google C++ style guide](https://google.github.io/styleguide/cppguide.html).

Students will be assigned to one of five teams which are **blue**, **brown**, **orange**, **pink** and **white**. Each team will have a separate branch and only team members will have write permission. Although, additional branches can be created if needed, but should have names starting with the team name (e.g. `blue_feature_one`).

## Objective

At the end of the project, students will have implemented several libraries which will enable alignment of a large amount of substrings (of various sizes) to a much larger string from which they originate. The objective is to join these libraries into a single program, often called mapper, in order to map long erroneous fragments from third generation of sequencing technologies to a reference genome, which has various use cases in bioinformatics. A visual example can be seen bellow.

![](misc/sample_mappings.png)

## Setup

Each team's main branch should be up to date with this README. This can be achieved by merging the master branch or rebasing onto it. The setup of the projects consists of creating a program and naming it in form of `<team name>_mapper` (e.g. `blue_mapper`). The program has to accept two files as floating arguments and enable options `-h` (`--help`) and `-v` (`--version`), which are used for the help and version messages, respectively. Suggested argument parser to include is `optarg`, but this feature can also be implemented independently.

The first file will contain a set of fragments in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) or [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format, while the second one will contain a corresponding reference genome in FASTA format. The files need to be parsed and stored in memory, and some statistics have to be outputted to `stderr` which includes number of sequences, average length, minimal and maximal length etc. There is no need to implement a parser. Students can  add [bioparser](https://github.com/rvaser/bioparser) to the project as a submodule via `git` and integrate it with `cmake`. It supports several bioinformatics formats where files can also be compressed with `gzip`. Therefore, before passing file paths with `bioparser` check whether file extension is one in the following list: `.fasta`, `.fa`, `.fastq`, `.fq`, `.fasta.gz`, `.fa.gz`, `.fastq.gz`, `.fq.gz`.

Sample program runs, after the setup step is completed, can be seen bellow:

```bash
blue_mapper escherichia_coli_r7_reads.fastq escherichia_coli_reference.fasta
<basic statistics of input files>
```

```bash
blue_mapper -h
<appropriate message describing supported arguments>
```

```bash
blue_mapper -v
v0.1.0
```

## Disclaimer

Laboratory for Bioinformatics and Computational Biology cannot be held responsible for any copyright infringement caused by actions of students contributing to any of its repositories. Any case of copyright infringement will be promptly removed from the affected repositories and reported to appropriate faculty organs.
