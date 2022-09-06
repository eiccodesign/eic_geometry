'''
    An example option file to digitize/reconstruct/clustering calorimeter hits
'''
from Gaudi.Configuration import *
import json
import os
import ROOT

from Configurables import ApplicationMgr, EICDataSvc, PodioInput, PodioOutput, GeoSvc
from GaudiKernel.SystemOfUnits import MeV, GeV, mm, cm, mrad

detector_name = str(os.environ.get("JUGGLER_DETECTOR", "endcapP_insert"))
detector_path = str(os.environ.get("DETECTOR_PATH", "${INSERT_PATH}"))
compact_path = os.path.join(detector_path, detector_name)

# input arguments from calibration file

# get sampling fractions from system environment variable, 1.0 by default
ci_hcal_sf = float(os.environ.get("CI_HCAL_SAMP_FRAC", 1.))
ci_hcal_insert_sf = float(os.environ.get("CI_HCAL_INSERT_SAMP_FRAC", 1.))

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
from Configurables import Jug__Reco__CalorimeterHitsMerger as CalHitsMerger
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG

from Configurables import Jug__Fast__InclusiveKinematicsTruth as InclusiveKinematicsTruth

# branches needed from simulation root file
sim_coll = [
    "MCParticles",
    "HcalEndcapPHits",
    "HcalEndcapPInsertHits",
]

# input and output
podin = PodioInput("PodioReader", collections=sim_coll)
podout = PodioOutput("out", filename=output_rec)

# Hadron endcap HCal
ci_hcal_daq = dict(
         dynamicRangeADC=200.*MeV,
         capacityADC=32768,
         pedestalMean=400,
         pedestalSigma=10)
ci_hcal_digi = CalHitDigi("ci_hcal_digi",
         inputHitCollection="HcalEndcapPHits",
         outputHitCollection="HcalEndcapHitsDigi",
         **ci_hcal_daq)
ci_hcal_reco = CalHitReco("ci_hcal_reco",
        inputHitCollection=ci_hcal_digi.outputHitCollection,
        outputHitCollection="HcalEndcapPHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_hcal_sf,
        **ci_hcal_daq)

ci_hcal_merger = CalHitsMerger("ci_hcal_merger",
        inputHitCollection=ci_hcal_reco.outputHitCollection,
        outputHitCollection="HcalEndcapPHitsRecoXY",
        readoutClass="HcalEndcapPHits",
        fields=["layer", "slice"],
        fieldRefNumbers=[1, 0])

ci_hcal_cl = IslandCluster("ci_hcal_cl",
        inputHitCollection=ci_hcal_merger.outputHitCollection,
        outputProtoClusterCollection="HcalEndcapPProtoClusters",
        splitCluster=False,
        minClusterCenterEdep=30.*MeV,
        localDistXY=[15.*cm, 15.*cm])

ci_hcal_clreco = RecoCoG("ci_hcal_clreco",
        inputProtoClusterCollection=ci_hcal_cl.outputProtoClusterCollection,
        outputClusterCollection="HcalEndcapPClusters",
        logWeightBase=6.2)

# Hadron Endcap HCal Insert
ci_hcal_insert_daq = dict(
         dynamicRangeADC=200.*MeV,
         capacityADC=32768,
         pedestalMean=400,
         pedestalSigma=10)
ci_hcal_insert_digi = CalHitDigi("ci_hcal_insert_digi",
         inputHitCollection="HcalEndcapPInsertHits",
         outputHitCollection="HcalEndcapPInsertHitsDigi",
         **ci_hcal_insert_daq)
ci_hcal_insert_reco = CalHitReco("ci_hcal_insert_reco",
        inputHitCollection=ci_hcal_insert_digi.outputHitCollection,
        outputHitCollection="HcalEndcapPInsertHitsReco",
        thresholdFactor=0.0,
        samplingFraction=ci_hcal_insert_sf,
        **ci_hcal_insert_daq)
# ci_hcal_insert_merger = CalHitsMerger("ci_hcal_insert_merger",
#         inputHitCollection=ci_hcal_insert_reco.outputHitCollection,
#         outputHitCollection="HcalEndcapPInsertHitsRecoXY",
#         readoutClass="HcalEndcapPInsertHits",
#         fields=["layer", "slice"],
#         fieldRefNumbers=[1, 0])

ci_hcal_insert_cl = IslandCluster("ci_hcal_insert_cl",
        inputHitCollection=ci_hcal_insert_reco.outputHitCollection,
        outputProtoClusterCollection="HcalEndcapPInsertProtoClusters",
        splitCluster=False,
        minClusterCenterEdep=30.*MeV,
        localDistXY=[15.*cm, 15.*cm])

ci_hcal_insert_clreco = RecoCoG("ci_hcal_insert_clreco",
        inputProtoClusterCollection=ci_hcal_insert_cl.outputProtoClusterCollection,
        outputClusterCollection="HcalEndcapPInsertClusters",
        logWeightBase=3.6)

# Truth level kinematics
truth_incl_kin = InclusiveKinematicsTruth("truth_incl_kin",
        inputMCParticles = "MCParticles",
        outputInclusiveKinematics = "InclusiveKinematicsTruth"
)

# Output
podout.outputCommands = ['drop *',
        'keep MCParticles',
        'keep *Digi',
        'keep *Reco*',
        'keep *Cluster*',
	'keep Inclusive*']

ApplicationMgr(
    TopAlg = [podin,
            ci_hcal_digi, ci_hcal_reco, ci_hcal_merger, ci_hcal_cl, ci_hcal_clreco,
            ci_hcal_insert_digi, ci_hcal_insert_reco, ci_hcal_insert_cl, ci_hcal_insert_clreco,
	    truth_incl_kin, podout],
    EvtSel = 'NONE',
    EvtMax = n_events,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
)
        