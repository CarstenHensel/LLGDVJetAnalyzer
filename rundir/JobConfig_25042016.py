import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('inputfile',
          '',
          VarParsing.multiplicity.list,
          VarParsing.varType.string,
          "Input File")

options.parseArguments()


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12', '')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:./LLGDV_MiniAOD_withNonCHSJets_1.root',
    )
)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates') 

process.demo = cms.EDAnalyzer('LLGDVJetAnalyzer',
    IsData = cms.bool(False),
    StoreGenData = cms.bool(False),
    RunLeptonTriggers = cms.bool(True),
    ignoreTriggers = cms.bool(False),
    bits = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjects = cms.InputTag("selectedPatTrigger"),
    METFilters = cms.InputTag("TriggerResults", "","PAT"),
    pupInfo = cms.InputTag("slimmedAddPileupInfo"),
    conversions = cms.InputTag('allConversions'),
    GenEventInfo = cms.InputTag("generator"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secVertices = cms.InputTag("slimmedSecondaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    jetsnochs = cms.InputTag("ak4PFJets", "", "Demo"),
    genJets = cms.InputTag("slimmedGenJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    pfCands = cms.InputTag("packedPFCandidates" ),
)


process.p = cms.Path(process.ak4PFJets * process.demo )
