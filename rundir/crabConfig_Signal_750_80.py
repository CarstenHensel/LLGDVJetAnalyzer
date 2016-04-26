from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'LLGDV_Analysis_Signal_750_80' 
config.General.workArea = 'crabJobs'
config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_25042016.py'
config.JobType.outputFiles = ['RecoOutput.root']
config.Data.inputDataset = '/CRAB_PrivateMC/mhamer-LLG_750_80_Fall15_stepPAT-17d438ff51ec6b3cada9e499a5a978e0/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'LLGDV_Analysis_AnalysisNTuple_Signal_750_80_Fall15'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
