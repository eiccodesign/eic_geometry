# EIC geometry for eiccodesign

## Installation instructions

This repository relies on the ip6 repository (commit 45ed8e8a77862dd727aabe5aad943e4adccec020). Make sure to clone and install both repositories. This is done automatically in the `get_frameworks.sh` script in `generate_data`.

## Editing the simulation
By default, the simulation includes just the beampipe and the hadron endcap with HCal components: HCal insert (W/Sc + Steel/Sc) & HCal (ATHENA: 20/3 mm Steel/Sc). Files are included for a homogeneous ECal (`compact/ecal/ecal_forward_homogeneous.xml`), a homogeneous ECal insert (`compact/ecal/ecal_forward_insert_homogeneous.xml`), and the ECCE LFHCAL (`compact/hcal/hcal_forward_ECCE.xml`). 

To change what detectors are simulated, add/remove the desired components in `hadron_endcap.xml`.

Some simple parameters for the geometry are contained in `compact/configuration_default.xml`. If you want to simulate without a beampipe and hole, replace `<include ref="compact/configuration_default.xml"/>` with `<include ref="compact/configuration_nohole.xml"/>` and comment out `<include ref="ip6/central_beampipe.xml"/>` in `hadron_endcap.xml`.

#### NOTE: If you adjust any compact or src files, you need to `make install` in your build directory before running the simulation.

## Running the simulation
`scripts/run_sim_hepmc.sh` generates a HepMC file and feeds it to npsim and DD4hep. The resulting sim file is then sent through Juggler for digitization, reconstruction, and clustering. The sim and reco files are saved.

To run the simulation, use `$DETECTOR_PATH/scripts/run_sim_hepmc.sh` after sourcing `setup_env.sh` from `generate_data`.

Some basic, adjustable paramaters are listed at the top of `run_sim_hepmc.sh` and the particle type, momentum, and number of events can be fed in with options `-part`, `-p`, and `-n`, respectively.

The simulation can be run from any directory and the output data will be stored in your current working directory. `scripts/hadron_endcap_reco.py` controls the digitization, reconstruction, and clustering. The clustering for the HCal insert isn't completely correct curently and needs adjusting (it will still run and generate output files though).

#### NOTE: If you remove/add any detector components from the simulation, you must also remove/add the related lines in `scripts/hadron_endcap_reco.py`

## Output
There will be two output files: a sim file and a reco file. The sim file contains the Geant4 level information while the reco file contains the digitized and reconstructed information. Both contain information about the MCParticles. Use the `XHitsReco` and `XClusters` branches in the reco files, where `X` is a detector name.  

<!--
## Addressing Homogenous ECal
To avoid high memory usage, the ECal and ECal insert use a mixed material of W/Polystyrene where the weight percentages are calculated based on the empirical weight and density of a prototype W/ScFi ECal. To incorporate the sampling fraction of the W/ScFi, a smearing procedure is needed: $E_{tower} = \mathrm{gRandom->Gaus}\left(E_{tower}*.03, \sigma\right)$, where $\sigma = E_{tower}\sqrt{a^2/E_{tower} + b^2}$, $a = 0.1$, and $b = 0.0015$. $\mathrm{gRandom}$ here is ROOT's random generator. This maintains the mean of W/ScFi and reproduces fluctuations with a random Gaussian. See [this presentation](https://github.com/rymilton/eic_endcap_insert/files/9172710/Smearing.of.mixture.structure.for.EMCal.pdf) (especially slide 3) by Zhiwan Xu, et al. for more details. -->

## Images of simulation geometry
<!--
Whole endcap:

<img src="https://user-images.githubusercontent.com/87345122/180581581-c85ece7d-7137-4392-ada6-f52dbaaec1ff.png" width="450"> <img src="https://user-images.githubusercontent.com/87345122/180581647-cc728b5b-3e67-48fc-bd74-570ddc7933d0.png" width="404"> -->

Front of HCal and insert:

<img src="https://user-images.githubusercontent.com/87345122/180581758-455f3b03-5d8b-4bf1-83a3-18e5a30dafdc.png" width="450">

W/Steel insert:

<img src="https://user-images.githubusercontent.com/87345122/180581614-37ec62c5-132e-4979-8864-8f8996bd86f7.png" width="800">


<!-- ## GPS documentation
---------------------------------
If you want to adjust the particle gun, here's the documentation for the general particle source (GPS):

   Manual: https://www.fe.infn.it/u/paterno/Geant4_tutorial/slides_further/GPS/GPS_manual.pdf

   Examples: https://hurel.hanyang.ac.kr/Geant4/Geant4_GPS/reat.space.qinetiq.com/gps/examples/examples.html -->
