# Bayesembler

The Bayesembler is a Bayesian method for doing transcriptome assembly from RNA-seq data. A manuscript describing the method has been accepted for publication (a reference will follow here upon publication).

## Installing the Bayesembler
The Bayesembler runs on both Linux and OS X and is freely available under the [MIT License](http://opensource.org/licenses/MIT).
### Static builds 
* Linux x86_64 ([latest release](https://github.com/bioinformatics-centre/bayesembler/releases/download/v1.1.1/Bayesembler_v1.1.1_Linux_x86_64.tar.gz))
* OS X x86_64 (will be available soon)

### Building the Bayesembler from source using CMake
1. Install [CMake](http://www.cmake.org/)
1. Install the [Boost](http://www.boost.org/), [Eigen](http://eigen.tuxfamily.org/) and [Bamtools](https://github.com/pezmaster31/bamtools) libraries
1. Install the [CEM](http://alumni.cs.ucr.edu/~liw/cem.html) assembler and [samtools](http://www.htslib.org/) and set the environment variables CEM_PROCESSSAM_PATH and SAMTOOLS_PATH to the locations of the processsam (part of CEM) and samtools binaries, respectively
1. Download the Bayesembler source code ([latest release](https://github.com/bioinformatics-centre/bayesembler/archive/v1.1.1.tar.gz)) or `git clone https://github.com/bioinformatics-centre/bayesembler.git`
1. Add a `CMakeLists.txt` file to your `src` directory ([example Bayesembler CMakeLists.txt](https://github.com/bioinformatics-centre/bayesembler/wiki/CMakeLists.txt-example))
1. `mkdir build` next to your `src` directory
1. In `build`, `cmake ../src` and `make`

## Running the Bayesembler

###Input
The Bayesembler requires *paired-end* RNA-seq data mapped using [TopHat2](http://ccb.jhu.edu/software/tophat/index.shtml).

###Assembly
The Bayesembler can be invoked with default settings using
```
bayesembler -b <tophat2_map.bam>
```
Important options include:
* `-s`: Data is strand-specific. Use `-s first` to indicate mate orientation as in the dUTP protocol or `-s second` to indicate opposite orientation. If not provided, the program assumes unstranded data.
* `-p`: Number of threads used for assembly.

###Output
Assembly in [Gene Annotation Format](http://genome.ucsc.edu/FAQ/FAQformat.html#format4) (GTF) as `assembly.gtf`.

Please consult the [wiki](https://github.com/bioinformatics-centre/bayesembler/wiki) for additional options and information about the output.

## Bug reports, help and feature requests
Please use the [GitHub issue tracking system](https://github.com/bioinformatics-centre/bayesembler/issues) to report bugs and request help or new features. Please use [bayesembler@binf.ku.dk](mailto:bayesembler@binf.ku.dk) for private communication only.

## Citation and credits
A manuscript describing the method has been accepted for publication; a reference to the paper will follow here as soon as it is published. The Bayesembler is being developed by Jonas Andreas Sibbesen, Lasse Maretty and Anders Krogh at the [Bioinformatics Centre](http://www.binf.ku.dk/), a part of the [Section for Computational and RNA Biology](http://www1.bio.ku.dk/binf/) at the [University of Copenhagen](http://www.ku.dk/english).
