from Gaudi.Configuration import *
import json
import os
import ROOT

from Configurables import ApplicationMgr, EICDataSvc, PodioInput, PodioOutput, GeoSvc
from GaudiKernel.SystemOfUnits import MeV, GeV, mm, cm, mrad, ns

detector_name = str(os.environ.get("JUGGLER_DETECTOR", "${DETECTOR}"))
detector_path = str(os.environ.get("DETECTOR_PATH", "${DETECTOR_PATH}"))
compact_path = os.path.join(detector_path, detector_name)

# input arguments from calibration file

# get sampling fractions from system environment variable, 1.0 by default
# ci_hcal_sf = float(os.environ.get("CI_HCAL_SAMP_FRAC", 1.))
# ci_hcal_insert_sf = float(os.environ.get("CI_HCAL_INSERT_SAMP_FRAC", 1.))

ci_hcal_sf = "1."
ci_hcal_insert_sf = "1."
ci_ecal_sf = "0.03"
ci_ecal_insert_sf = "0.03"

# input and output
input_sims = [f.strip() for f in str.split(os.environ["JUGGLER_SIM_FILE"], ",") if f.strip()]
output_rec = str(os.environ["JUGGLER_REC_FILE"])
n_events = int(os.environ["JUGGLER_N_EVENTS"])

# geometry service
geo_service = GeoSvc("GeoSvc", detectors=["{}.xml".format(compact_path)], OutputLevel=INFO)
# data service
podioevent = EICDataSvc("EventDataSvc", inputs=input_sims)


# juggler components
from Configurables import Jug__Digi__CalorimeterHitDigi as CalHitDigi
from Configurables import Jug__Reco__CalorimeterHitReco as CalHitReco

# from Configurables import Jug__Fast__InclusiveKinematicsTruth as InclusiveKinematicsTruth

# branches needed from simulation root file
sim_coll = [
    "MCParticles",
    "HcalEndcapPHits",
    "HcalEndcapPInsertHits",
#     "EcalEndcapPHits",
#     "EcalEndcapPInsertHits"
]

# input and output
podin = PodioInput("PodioReader", collections=sim_coll)
podout = PodioOutput("out", filename=output_rec)

# Hadron endcap HCal
ci_hcal_daq = dict(
         dynamicRangeADC=800.*MeV,
         capacityADC=32768,
         pedestalMean=400,
         pedestalSigma=10)
ci_hcal_digi = CalHitDigi("ci_hcal_digi",
         inputHitCollection="HcalEndcapPHits",
         outputHitCollection="HcalEndcapHitsDigi",
         timeResolution=1.*ns,
         **ci_hcal_daq)
ci_hcal_reco = CalHitReco("ci_hcal_reco",
        inputHitCollection=ci_hcal_digi.outputHitCollection,
        outputHitCollection="HcalEndcapPHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_hcal_sf,
        **ci_hcal_daq)

# Hadron Endcap HCal Insert
ci_hcal_insert_daq = dict(
         dynamicRangeADC=800.*MeV,
         capacityADC=32768,
         pedestalMean=400,
         pedestalSigma=10)
ci_hcal_insert_digi = CalHitDigi("ci_hcal_insert_digi",
         inputHitCollection="HcalEndcapPInsertHits",
         outputHitCollection="HcalEndcapPInsertHitsDigi",
         timeResolution=1.*ns,
         **ci_hcal_insert_daq)
ci_hcal_insert_reco = CalHitReco("ci_hcal_insert_reco",
        inputHitCollection=ci_hcal_insert_digi.outputHitCollection,
        outputHitCollection="HcalEndcapPInsertHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_hcal_insert_sf,
        **ci_hcal_insert_daq)

# Hadron Endcap ECal
# ci_ecal_daq = dict(
#          dynamicRangeADC=3.*GeV,
#          capacityADC=65536,
#          pedestalMean=100,
#          pedestalSigma=0.7)
# ci_ecal_digi = CalHitDigi("ci_ecal_digi",
#          inputHitCollection="EcalEndcapPHits",
#          outputHitCollection="EcalEndcapHitsDigi",
#          scaleResponse=ci_ecal_sf,
#          energyResolutions=[0.00340,0.0009,0.0],
#          timeResolution=1.*ns,
#          **ci_ecal_daq)
# ci_ecal_reco = CalHitReco("ci_ecal_reco",
#         inputHitCollection=ci_ecal_digi.outputHitCollection,
#         outputHitCollection="EcalEndcapPHitsReco",
#         thresholdFactor=5.0,
#         thresholdValue=2.0,
#         samplingFraction=ci_ecal_sf,
#         **ci_ecal_daq)

# Hadron Endcap ECal insert
# Uncomment this section when the ecal insert is simulated
# ci_ecal_insert_daq = dict(
#          dynamicRangeADC=3.*GeV,
#          capacityADC=65536,
#          pedestalMean=100,
#          pedestalSigma=0.7)
# ci_ecal_insert_digi = CalHitDigi("ci_ecal_insert_daq",
#          inputHitCollection="EcalEndcapPInsertHits",
#          outputHitCollection="EcalEndcapInsertHitsDigi",
#          scaleResponse=ci_ecal_insert_sf,
#          energyResolutions=[0.00340,0.0009,0.0],
#          timeResolution=1.*ns,
#          **ci_ecal_insert_daq)
# ci_ecal_insert_reco = CalHitReco("ci_ecal_insert_reco",
#         inputHitCollection=ci_ecal_insert_digi.outputHitCollection,
#         outputHitCollection="EcalEndcapPInsertHitsReco",
#         thresholdFactor=5.0,
#         thresholdValue=2.0,
#         samplingFraction=ci_ecal_insert_sf,
#         **ci_ecal_insert_daq)

# Output
podout.outputCommands = ['drop *',
        'keep MCParticles',
        'keep *Digi',
        'keep *Reco*']

ApplicationMgr(
    TopAlg = [podin,
            ci_hcal_digi, ci_hcal_reco, 
            ci_hcal_insert_digi, ci_hcal_insert_reco,
        #     ci_ecal_digi, ci_ecal_reco,
        # Uncomment this when using the ecal insert
        #     ci_ecal_insert_digi, ci_ecal_insert_reco, 
	    podout],
    EvtSel = 'NONE',
    EvtMax = n_events,
    ExtSvc = [podioevent]
)
        