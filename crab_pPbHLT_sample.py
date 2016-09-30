from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'YOUR_REQUEST_NAME'
config.General.workArea = 'YOUR_WORKAREA_NAME'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'YOUR_CONFIG'
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.useParent = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.inputDataset = 'RECO OR AOD INPUT'
config.Data.unitsPerJob = 5000
config.Data.outLFNDirBase = 'YOUR_output_dir'
config.Data.outputDatasetTag = 'YOUR_output_name'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
