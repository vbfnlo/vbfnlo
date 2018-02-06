# VBFNLO

* Version: 3.0.0 beta 5
* Release date:  06 February 2018

---

* [arXiv:1107.4038], [arXiv:1207.4975], [arXiv:1404.3940]
* http://www.itp.kit.edu/vbfnlo
* vbfnlo@itp.kit.edu

---

  VBFNLO is a fully flexible parton level Monte Carlo program for the
simulation of vector boson fusion, double and triple vector boson production in
hadronic collisions at next-to-leading order in the strong coupling constant. 
Electroweak corrections to Higgs production via vector boson fusion are also 
included in the SM and MSSM.  VBFNLO includes Higgs and vector boson decays
with full spin correlations and all off-shell effects.  In addition, VBFNLO 
implements CP-even and CP-odd Higgs boson production via gluon fusion, associated 
with two jets, at the leading-order one-loop level.

  A variety of effects arising from beyond the Standard Model physics are 
implemented for selected processes. Processes can be run in the MSSM, anomalous
couplings of Higgs and vector bosons can be used, a Three-Site Higgsless model 
and a Warped Higgsless extra dimension model are included.  The program offers the 
possibility of generating Les Houches Accord event files for most processes 
available at leading order.

  Complete documentation, including instructions for building, installing
and running VBFNLO, is available in the doc/Manual.pdf file.


## Installation

* `autoreconf -vi` (only needed for git checkouts)
* `./configure --OPTIONS` see below
* `make`
* `make install`

## Running VBFNLO

To run VBFNLO:
`./bin/vbfnlo --input=PATH/TO/INPUT_FILES`

If the input -- e.g. vbfnlo.dat -- files are in the current directory, the
--input can of course be omitted.


### Configuration Options

VBFNLO can be configured with several options, a complete list of which can be
found in the [manual](https://www.itp.kit.edu/vbfnlo/wiki/doku.php?id=documentation:manual)
or by running `./configure --help`. 

The main options are:

* `--prefix=/path/to/install/dir`: Install VBFNLO in the location given by path.

* `--enable-processes=LIST`: Comma-separated list of processes to enable. Options
                          are: vbf, ggf, diboson, dibosonjet, triboson,
                          tribosonjet, hjjj, top, qcdvvjj, qcdvjj,
                          all_except_hexagons, all. 
                          Default: all_except_hexagons

* `--enable-kk`: enable simulation of Kaluza-Klein resonances
* `--enable-MPI`: use MPI parallelization
* `--disable-NLO`: disable next-to-leading order QCD
* `--enable-madgraph`: include code for MadGraph comparisons
* `--enable-quad`: enable quad precision for difficult phase space points

Optional Packages:
* `--with-LHAPDF=DIR`: location of LHAPDF installation
* `--with-LOOPTOOLS=/path/to/LOOPTOOLS/`: needed for the calculation of electroweak corrections
* `--with-FEYNHIGGS=/path/to/FEYNHIGGS/`: needed for the calculation of MSSM Higgs sector by FeynHiggs
* `--with-root=DIR`: location of ROOT installation
* `--with-hepmc=DIR`: location of HepMC installation
* `--with-gsl=DIR`: location of gsl installation, default=system lib path

## Bug reports
Please report any problems to vbfnlo@itp.kit.edu
with a short report including the `configure` options used to build
VBFNLO, as well as the versions of compilers and external libraries
used.

You can also open issues in github: https://github.com/vbfnlo/vbfnlo/issues


[arXiv:1107.4038]: https://arxiv.org/abs/1107.4038
[arXiv:1207.4975]: https://arxiv.org/abs/1207.4975
[arXiv:1404.3940]: https://arxiv.org/abs/1404.3940
