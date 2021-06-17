nopilesum
=========

This tool functions as an alternative to GATK4's GetPileupSummaries which has very high memory usage that may be incompatible with some setups. It attempts to output results that are the same as GATK's tool, though in my testing I ahve noticed some discrepencies. Though my validation leads me to believe there may be bugs in the GATK program. In any case, `nopilesum` is a *D* program that relies on [dhtslib](https://github.com/blachlylab/dhtslib) and [htslib](https://github.com/samtools/htslib). `nopilesum`'s name is derivative of the fact that we avoid pileup and perform allele counting to achieve a `nopile` summary that is similar to GATK's result with minimal memory usage.

## Installation
### Dependencies
#### htslib
To intall htslib dependencies:
```
Debian / Ubuntu
---------------

sudo apt-get update  # Ensure the package list is up to date
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

Note: libcurl4-openssl-dev can be used as an alternative to libcurl4-gnutls-dev.

RedHat / CentOS / Amazon Linux
---------------

sudo yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel

Alpine Linux
------------

sudo apk update  # Ensure the package list is up to date
sudo apk add autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev libressl-dev

OpenSUSE
--------

sudo zypper install autoconf automake make gcc perl zlib-devel libbz2-devel xz-devel libcurl-devel libopenssl-devel

MacOS
-----

brew install xz autoconf automake
```
You will then need to download [htslib](http://www.htslib.org/download/).
As of now we support versions >=1.10. To install htslib:
```
curl https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -o htslib-1.12.tar.bz2
tar -xjf htslib-1.12.tar.bz2
cd htslib-1.12
./configure
make 
sudo make install
```
### Grab a release
You can find a linux binary [here]().

Or you can build nopilesum yourself. This will require a *D* compiler.

### D compiler
To install *D* run the following command to install into your home folder. 
```
curl https://dlang.org/install.sh | bash -s ldc
```
To activate the *D* environment: `source ~/dlang/ldc-*/activate`

You could use `dmd` instead of `ldc`.
```
curl https://dlang.org/install.sh | bash -s dmd
```
We prefer `ldc` for its better performance and it is the compiler we actively test `recontig` with.
Your results may vary.

If you already have htslib 1.12 installed on your system and on your path then activate your *D* environment
```
# LD_LIBRARY_PATH and LIBRARY_PATH may not be required if 
# you already have them defined
LD_LIBRARY_PATH=/usr/local/lib/ LIBRARY_PATH=/usr/local/lib/ dub build
```

If your htslib install is in a non-standard location:
```
LD_LIBRARY_PATH=path/to/htslib/ LIBRARY_PATH=path/to/htslib/ dub build
```

## Usage
`nopilesum` inputs can be BAM/SAM or VCF/BCF respectively. `nopilesum` can automatically resolve and open remote files, gzipped files, and bgzipped files.
```
./nopilesum <in.bam> <in.vcf> > summary.txt
```

### Cmdline
```
nopilesum: get pileup summaries
usage: nopilesum <in.bam> <in.vcf> > summary.txt
    --max-af VCF records with INFO/AF fields above this threshold won't be used (default 0.2)
    --min-af VCF records with INFO/AF fields below this threshold won't be used (default 0.01)
-v --verbose see warnings
     --debug see extra info
-h    --help This help information.
```
