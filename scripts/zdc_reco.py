from Gaudi.Configuration import *
import json
import os
import ROOT

from Configurables import ApplicationMgr, EICDataSvc, PodioInput, PodioOutput, GeoSvc
from GaudiKernel.SystemOfUnits import MeV, GeV, mm, cm, mrad, ns, ps

detector_name = str(os.environ.get("JUGGLER_DETECTOR", "${DETECTOR}"))
detector_path = str(os.environ.get("DETECTOR_PATH", "${DETECTOR_PATH}"))
compact_path = os.path.join(detector_path, detector_name)

# input arguments from calibration file

# get sampling fractions from system environment variable, 1.0 by default
# ci_hcal_sf = float(os.environ.get("CI_HCAL_SAMP_FRAC", 1.))
# ci_hcal_insert_sf = float(os.environ.get("CI_HCAL_INSERT_SAMP_FRAC", 1.))

ci_zdc_hcal_sf = "1."
ci_zdc_ecal_sf = "1."

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
    "ZDCHcalHits",
#     "ZDCEcalHits"
]

# input and output
podin = PodioInput("PodioReader", collections=sim_coll)
podout = PodioOutput("out", filename=output_rec)

# Hadron endcap HCal
ci_zdc_hcal_daq = dict(
        dynamicRangeADC=800.*MeV,
        capacityADC=32768,
        pedestalMean=400,
        pedestalSigma=10)
ci_zdc_hcal_digi = CalHitDigi("ci_zdc_hcal_digi",
        inputHitCollection="ZDCHcalHits",
        outputHitCollection="ZDCHcalHitsDigi",
        timeResolution=1.*ns,
        **ci_zdc_hcal_daq)
ci_zdc_hcal_reco = CalHitReco("ci_zdc_hcal_reco",
        inputHitCollection=ci_zdc_hcal_digi.outputHitCollection,
        outputHitCollection="ZDCHcalHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_zdc_hcal_sf,
        **ci_zdc_hcal_daq)

# ci_zdc_ecal_daq = dict(
#         dynamicRangeADC=2000.*MeV,
#         capacityADC=32768,
#         pedestalMean=400,
#         pedestalSigma=3.2)
# ci_zdc_ecal_digi = CalHitDigi("ci_zdc_ecal_digi",
#         inputHitCollection="ZDCEcalHits",
#         outputHitCollection="ZDCEcalHitsDigi",
#         timeResolution=0.*ns,
#         resolutionTDC = 10*ps,
#         scaleResponse = 1.0,
#         **ci_zdc_ecal_daq)
# ci_zdc_ecal_reco = CalHitReco("ci_zdc_ecal_reco",
#         inputHitCollection=ci_zdc_ecal_digi.outputHitCollection,
#         outputHitCollection="ZDCEcalHitsReco",
#         thresholdFactor=4.0,
#         thresholdValue = 0.0,
#         samplingFraction=ci_zdc_ecal_sf,
#         **ci_zdc_ecal_daq)

# Output
podout.outputCommands = ['drop *',
        'keep MCParticles',
        'keep *Digi',
        'keep *Reco*']

ApplicationMgr(
    TopAlg = [podin,
            ci_zdc_hcal_digi, ci_zdc_hcal_reco,
        #     ci_zdc_ecal_digi, ci_zdc_ecal_reco,
	    podout],
    EvtSel = 'NONE',
    EvtMax = n_events,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
)
        