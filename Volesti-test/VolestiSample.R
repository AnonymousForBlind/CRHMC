# % 1 Abiotrophia_defectiva_ATCC_49176.mat: (952, 1069) -> 157
# % 2 Acidaminococcus_fermentans_DSM_20731.mat: (1009, 1090) -> 164
# % 3 Acidaminococcus_intestini_RyC_MR95.mat: (917, 994) -> 123
# % 4 Acidaminococcus_sp_D21.mat: (856, 851) -> 103
# % 5 Acinetobacter_calcoaceticus_PHEA_2.mat: (1319, 1561) -> 328
#   % 6 AntCore.mat: (77, 90) - Empty A
# % 15 cardiac_mit_glcuptake_atpmax.mat: (230, 220) -> 12
#   % 7 Ec_iAF1260_flux1.mat: (1668, 2382) -> 524 - duplicate
# % 16 ecoli_core_model.mat: (72, 95) -> 24
# % 17 iAF1260.mat: (1668, 2382) -> 524
# % 18 iJO1366.mat: (1805, 2583) -> 582
# % 19 modelReg.mat: (72, 95) -> 24
# % 9 Recon1.0model.mat: (2766, 3742) -> 932
# % 10 Recon2.0model.mat: (5063, 7440) -> 2430
#   % 11 Recon2.v04.mat: (5063, 7440) - duplicate
#   % 12 Recon2.v05.mat: (5063, 7440) - duplicate
#   % 13 Recon3DModel_301.mat: (5835, 10600) - duplicate 
#   % 14 Recon3D_301.mat: (8399, 13543) - too large

# pass: 1~7, 9, 11, 19
# session aborted: 8, 10, 12~15, 17~18, 20
# infinite loop?: 16
# not supported: 21

library(R.matlab)
library(volesti)
name_list = list.files("../Instances/2cdhr")

for (modelname in name_list[c(16)]){
  print(modelname)
  print(Sys.time())
  modelmat = readMat(paste("../Instances/2cdhr/",modelname,sep=''))
  polytope = modelmat$poly
  
  A = polytope[1];
  A = A[[1]];
  A = A[[1]];
  
  b = polytope[2]
  b = b[[1]]
  b = b[[1]]
  
  P = Hpolytope(A=A, b=as.numeric(b))
  dim = length(A[1,])
  print(dim)
  
  # Estimate how many samples can be drawn for 24hrs
  pre = system.time({points =  sample_points(P, n=1, random_walk = list("walk" = "CDHR", 
                                                                        "walk_length"=dim^2))})
  sample_num <- min(round(24*3600/pre[3]), 1000)
  print(sample_num)
  # Sample 
  tim = system.time({points =  sample_points(P, n=sample_num, random_walk = list("walk" = "CDHR", 
                                                                                 "walk_length"=dim^2))})
  filetosave <- list("dim" = length(A[1,]), "n" = sample_num, "time" = tim[3], "step" = dim^2*sample_num, "points" = points)
  writeMat(paste(getwd(), "/SavedPoints/Sampled",modelname,sep=''), saved = filetosave)
}
