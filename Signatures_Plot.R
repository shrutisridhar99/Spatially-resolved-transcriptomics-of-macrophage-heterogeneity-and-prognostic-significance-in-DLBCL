#Signature plot 2023
#Kevin Mulder

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(scCustomize)

MoMac = readRDS("MoMac_reference.rds")
DimPlot(MoMac)

DefaultAssay(MoMac) = "RNA"

# custom gene signature score of gene list from Signature_Plot_List.R
MoMac <- AddModuleScore(
  object = MoMac,
  features = list(Interfol_up2),
  name = 'Interfol_upV'
)

#Plot signatures

FeaturePlot_scCustom(MoMac, features = "Interfol_upV1", order = F)

FeaturePlot_scCustom(MoMac, features = "Relapse1", order = F)

FeaturePlot_scCustom(MoMac, features = "Dark_ZoneV1", order = F)

FeaturePlot_scCustom(MoMac, features = "Light_ZoneV1", order = F)
FeaturePlot_scCustom(MoMac, features = "Light_Zone1", order = F)


