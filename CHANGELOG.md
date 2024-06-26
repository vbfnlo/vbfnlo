# Changelog

## [3.0.0] - 2024-04-XX

### Major Enhancements
* VBFNLO website is now at https://ific.uv.es/vbfnlo/index.html (the old KIT website 
at https://www.itp.kit.edu/vbfnlo/wiki/doku.php?id=overview is no longer updated).
* VBFNLO can now be linked to Herwig 7 via the BLHA interface, for all VBS, double and triple vector
boson production processes.

### Updates and new features
* Added additional form factor types (2018-07).

### Bug fixes
* A bug has been fixed in the calculation of histogram-bin errors (2020-07). 

## [3.0.0 beta 5] - 2018-02-06

### Major Enhancements
* Large speed improvement in real-emission calculation of Hjjj production (processes 110-117)

### Bug fixes
* Fix two bugs in NLO calculation of Hjjj production (processes 110-117)
* Fix number of active flavours in BLHA interface

## [3.0.0 beta 4] - 2017-04-19

### Bug fixes
* Fix bug in NLO calculation of ZAjj production (processes 290,291) when including anomalous couplings
* Fix bug in event output for ZAjj production

## [3.0.0 beta 3] - 2017-01-18

### Major Enhancements
* VBFNLO is now published also at https://github.com/vbfnlo/vbfnlo
* Processes formerly generated by running `ggflo` now use the binary `vbfnlo` as well
* There is now a quick low-statistics version of the regress tests. They require the
pytest package on python3. Run them using `cd regress; py.test`

### Updates and new features
* Add dimension-8 operator O<sub>S,2</sub> and adjust K-matrix unitarization accordingly
* Support for parallelization using MPI via `--enable-MPI`. See manual for details.

### Other additions and bug fixes
* Fix bug in BLHA interface for same-flavour boson decays
* Fix PDF convolution in real-emission part of VBF-Zjj
* Fix bug in event output for QCD-ZZjj and -ZAjj
* Show timing information by default
* Add git hash of source version to runtime version output
* Fix compilation issue when flag `--disable-nlo` is used
* Fix warnings for automake version >= 1.15
* Tested with gfortran 5.4 and 6.2
* Redirect LoopTools output for VBF-H EW corrections to LoopTools.out
* Removed compilation flag `--enable-spin2`. The spin2 code is now always included at
compile time and can be switched on/off at runtime as before
* Reduce memory consumption of Histogram routine. Better scaling for many histograms.

## [3.0.0 beta 2] - 2015-10-30

### Updates and new features
* Add possibility for Breit-Wigner or Passarino CPS scheme distribution of Higgs boson 
(VBF-H only, process IDs 100-117,1010)

### Other additions and bug fixes
* Fix bug in K-matrix unitarisation of W^+ W^- jj

## [3.0.0 beta 1] - 2015-08-05

### Major enhancements (for details see manual)

* Several new processes are available:
WWj, ZZj production including the loop-induced gluon-fusion contribution
* Interface according to the Binoth Les Houches Accord (BLHA) added
for all VBF processes with fully leptonic decays
* K-matrix unitarisation procedure implemented for the two 
dimension 8 operators O_{S,0} and O_{S,1}.


## [2.7.1] - 2015-08-05 

* New processes:
QCD-ZZjj, QCD-Zgammajj, VBF-Zgammajj, Z(->l+l-)Z(->nu nubar)gamma
* Higgs decays for VBF-HHjj process
* Enable linking with LHAPDF v6
* Fix PDF convolution in QCD-(W-/W-Z/W-gamma)jj
* Fixed missing entry of VBF-HHjj as VBF process in cut routine
* Fixed missing axis labels on 2D histograms


2014-04-15 *** Changes in Version 2.7.0

  *** Major enhancements (for details see manual): ***

        * Several new processes are available:
          Wgammajj in VBF, HHjj in VBF, W, Wj, WH, WHj, pp->Spin2jj in VBF (with 
          Spin2->WW/ZZ->leptons), several QCD induced processes (WZjj, Wgammajj, same sign 
          WWjj, Wjj)
          Anomalous gauge boson couplings are supported in VBF Wgammajj, WH, WHj.
        * Semileptonic decay modes are available in the following diboson, triboson and 
          VBF diboson plus two jets processes: VV, VVV, VVgamma, VVjj with V=W,Z.
          Additionally, the Higgs production processes pp -> H -> WW/ZZ in VBF
          with semileptonic decays are available as well.
          The processes with semileptonic decays support anomalous gauge boson couplings.
        * The VBF VVjj processes can be run taking a second Higgs resonance into account
          (simplified 2HDM).
        * Anomalous gauge boson couplings are now available in all triboson (plus jet)
          processes and all VVjj production processes via VBF.
        * A larger set of operators leading to anomalous gauge boson couplings 
          is included in VBFNLO.
        * Added ability to use anomalous HVV couplings in Higgs production via
          gluon fusion, where the Higgs decays into WW or ZZ (prodIDs 4105-7)

  ***  Further updates and new features: ***

        * It is now possible to specify a desired number of unweighted events for event
          output in vbfnlo.dat
          In case VBFNLO is not able to achieve the desired number of unweighted
          events due to single events with very large weights it's also possible to
          opt for "partially unweighted events", where most of the weights will have
          weight 1, and a few will have larger weights.
        * Added 'any leptons' switch (LEPTONS = 98/99), allowing for all possible lepton
          final state combinations of 2 or 3 generations.  See manual for details.
        * m_ll cut can be switched between cutting only opposite-sign lepton combinations
          or all combinations
        * H->bb decay mode for Higgs boson production processes:
          jet cuts are applied to final state b-quarks
        * All default values of anomalous couplings are set to zero
          to avoid confusion
        * implemented staggered pT cuts for leptons and photons

  *** Other additions and bug fixes: ***

        * Phase space improvements to increase unweighting efficiency for event output in
          several VBF, di- and triboson processes
        * Phase space optimizations for Z decays in VV and VVjj

        * Path to procinfo.dat is no longer hardcoded, now it's also searched in the input
          directory
        * Increased allowed file name length for slha file
        * Print out values of g, g' used in calculation
        * in raw histogram data output: include gnuplot command files
        * New command line option --version : provides version information of VBFNLO
        * add workaround for enormous memory consumption during compilation with gfortran 4.8

        * event output for W-W-jj in VBF: color information fixed
        * Fixed bug in WWZ virtuals implementation -- new results are about a percent lower
        * Several dynamical renormalization and factorization scale choices corrected which 
          gave wrong results for the real emission part of NLO calculations:
          min(pT_jet) (ID 2), min(ET_bosons) (ID 7)
        * fix particle-antiparticle-assignment in W+jj, Zjj
        * event output for W-W+W-: fixed momenta assignment
        * Fixed line length for ifort
        * default value for cut MDIJ_MIN is changed to 0 GeV
        * works with LHAPDF version 6
        * output of per-bin error in histograms now available in the ROOT output
          of histograms defined in histograms.F


2013-06-19 *** Changes in Version 2.6.3

        * Fixed particle ID in LHE output for process 191 (Spin2->gamgam in VBF)
        * Fixed bug in couplings for spin-2 VBF->ZZ->llnunu (process 211)
        * Fixed bug with pdf sum for LO+j VBF processes
        * Fixed several bugs in les houches output concerning tau mass output and
          helicity output#
        * Added tau mass output for LHE to all processes
        * Fix linking HepMC on some newer distributions
        * In pp->Hjjj: fix PDF sum if nfl<4
        * In ggf Higgs processes: remove debug output when using event output
        * Process 191 is now treated as VBF process in JetVeto-cuts
        * LHE output for process W-W-jj: parton and beam particle IDs fixed
        * Add higher-order corrections to H->gg partial width
        * Fixed particle ID in LHE output for gg->VV processes
        * Fixed bugs in electroweak corrections to VBF-H
        * For event output, the default choice is now to generate unweighted events.
        * Consistent use of nflVBF for Hjjj processes


2012-12-07 *** Changes in Version 2.6.2

        * Fixed lepton assignment in H->(WW,ZZ)->4l decays.
        * Fixed event output for processes with more than one phase space.
        * For the Ajj production in VBF (procID 150) a different form factor is used with
          respect to the other processes.
        * Phase space optimizations in WWW, ZZW and ZZZ production
        * Fixed bug in dipole subtraction for ZZZ production
        * Phase space optimization for VBF->4l including spin-2 case
        * Improved output of event numbers written to event file
        * Faster code for process 191 (spin2 -> gamma gamma)


2012-08-24 *** Changes in Version 2.6.1

	* Fixed problem with TREEFAC for anomalous Higgs couplings
        * Fixed problem with root-only histograms (only histograms only in 
             rootusershist.cpp were affected)
        * Fixed coefficient of FB_ODD in a3_HZZ coupling
        * Fixed symmetry factor in 2106 (pp -> HAjj -> ZZAjj -> llllAjj) and 
             4106 (gluon fusion pp -> Hjj -> ZZjj -> lllljj) for identical final leptons
        * Increased allowed path length for input/output files


2012-07-20 *** Changes in Version 2.6.0

        *** New processes added (see manual for details):

        * W+W+jj and W-W-jj via VBF: process IDs 250-260
        * diboson processes WZ, ZZ, WA, ZA, AA: process IDs 310-370
        * triboson + jet processes WAAj: process IDs 800-810
        * gluon-induced diboson production WW, ZZ, ZA, AA: process IDs 4300-4370


	*** Updates and new features:

        * Anomalous couplings included in vector boson production via VBF: process IDs 120-150
        * Spin-2 model included in diboson + 2 jet production via VBF: process IDs 200-240
        * The switch FERMIONLOOP controls the inclusion of the gluon-induced higher order 
             processes to diboson production: process IDs 300-370, 4300-4370
        * Anomalous couplings (VVV and HVV) included in diboson production: process IDs 300-370,
             4300-4370


	*** Other additions and bug fixes:

        * Added input HWIDTH: the Higgs width can now be fixed by the user
        * When anomalous VVV couplings are used, the Higgs width is now calculated with the
            appropriate anomalous HVV couplings in relevant processes
        * Altered histograms to correct the errors and allow the user to control the 
            smearing between bins
        * Altered formfactor momentum dependence for anomalous gauge boson couplings 
        * Increased allowed width of virtuality of resonant bosons in phasespace generation 
            affected processes: 106,117,116,117,2107,4107,211
        * Fixed compilation error to allow both HepMC and ROOT to be linked
        * Streamlined calls to LHAPDF (no change to results)
        * Fixed compilation error with LHAPDF for newer versions of gcc
        * Fixed bug in print out of alfa_s for scale choice 6 (no change to results)
        * Fixed problem in TREEFAC for anomalous HVV couplings, and separated into 
            TREEFACZ and TREEFACW
        * Fixed problem in MASS_SCALE_KAPPA for anomalous VVV couplings - separated into
            MASS_SCALE_AKAPPA and MASS_SCALE_ZKAPPA
        * Fixed bug in HVV anomalous coupling formfactor F2
        * Input grid files are now taken from the --input=INPUT directory: output files still
            written to the working directory
        * Added ability to include electroweak loop squared contributions in SUSY 
            (switch LOOPSQR_SWITCH)
        * Added ability to switch between tree-level and corrected Higgs masses in electroweak
            loop corrections (switch MH_LOOPS)
        * Bug fixed in virtual contributions to H+3j processes: process IDs 110-117


	*** STRUCTURAL CHANGE:
       
        * The main files, vbfnlo_main.F and ggflo_main.F, have been moved from src/ to lib/
            This is to enable compilation with newer versions of gfortran

	
2011-09-08 *** Changes in Version 2.5.3
	
        * Fixed implementation of anomalous couplings in diboson plus jet 
            processes (610 - 640)
        * Fixed problem in root histogram for invisible particles (neutrinos)

2011-08-26 *** Changes in Version 2.5.2
	
	* Fixed write-out of alfa_strong
	* Fixed bug in HepMC output
	* Fixed bug in H -> ZZ -> 4l
	* Fixed bug in anomalous WWgamma coupling for VVV production
        * Stopped calculation of diboson plus jet processes with anomalous couplings 
            while we fix a bug in their implementation

2011-07-21 *** Changes in Version 2.5.1

	* Fixed bug in real emission part of VBF processes.

2011-07-20 *** New release: version 2.5.0 ***

	*** New processes added (see manual for details):

	* H + photon + jj: process IDs 2100-2107
	* photon + jj: process ID 150
	* W + photon + j, W + Z + j: process IDs 610-640
	* WW + photon, ZZ + photon: process IDs 460,470
	* WZ + photon: process IDs 480,490
	* W + 2 photons, Z + 2 photons: process IDs 500-521
	* 3 photons: process ID 530

	*** Updates to existing processes:
	* Anomalous quartic gauge boson couplings for VVV production
	* Extension to the MSSM (via FeynHiggs, SLHA input, user input)
	* Electroweak corrections to VBF Higgs production in the SM and MSSM (requires
	    link to LoopTools)
	* Triple vector production in Higgsless KK


	*** Other additions:
	* Ability to vary histogram range via input histogram.dat file
	* More cut options implemented
	* New drawing option


	*** Bug fixes
	* Bug fixed in readinput routine
	* Correct evolution of bottom mass for gluon fusion
	* Allow low pTj for Higgs production via VBF
	* Updates to LHAPDF interface
	* ROOT interface updated


2009-11-25 *** Changes in Version 2.0.3

	* Fixed bug and improved mapping in WWW phase space

2009-04-07 *** Changes in Version 2.0.2

	* Bug fixed in the real emission part of the WWjj process (VBF process)

2008-11-27 *** Release of Version 2.0.0

[Unreleased]: https://github.com/vbfnlo/vbfnlo/compare/v3.0.0beta4...HEAD
[3.0.0 beta 4]: https://github.com/vbfnlo/vbfnlo/compare/v3.0.0beta3...v3.0.0beta4
[3.0.0 beta 3]: https://github.com/vbfnlo/vbfnlo/compare/v3.0.0beta2...v3.0.0beta3
[3.0.0 beta 2]: https://github.com/vbfnlo/vbfnlo/compare/v3.0.0beta1...v3.0.0beta2
[3.0.0 beta 1]: https://github.com/vbfnlo/vbfnlo/compare/v2.7.1...v3.0.0beta1
[2.7.1]: https://github.com/vbfnlo/vbfnlo/compare/v2.7.0...v2.7.1
