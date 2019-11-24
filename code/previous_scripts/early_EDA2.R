meta_joint = readLines(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/AD_joint_meta.txt"))
v3_joint = readLines(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/AD_joint_v3.txt"))

meta_joint = readLines(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_joint_meta.txt"))
v3_joint = readLines(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_joint_v3.txt"))

meta_single = list()
fl = list.files(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_separate/meta"))
setwd(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_separate/meta"))
for(i in 1:length(fl)){
  meta_single[[i]] = readLines(fl[i])
}
v3_single = list()
fl = list.files(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_separate/v3"))
setwd(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_separate/v3"))
for(i in 1:length(fl)){
  v3_single[[i]] = readLines(fl[i])
}

#meta_single = list()
#fl = list.files(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/AD_separate/meta"))
#setwd(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/AD_separate/meta"))
#for(i in 1:length(fl)){
#  meta_single[[i]] = readLines(fl[i])
#}
v3_single = list()
fl = list.files(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/AD_separate/v3"))
setwd(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/AD_separate/v3"))
for(i in 1:length(fl)){
  v3_single[[i]] = readLines(fl[i])
}

max(unlist(lapply(v3_single,length)))
length(unique(unlist(v3_single[7:16])))
sum(unlist(lapply(v3_single,length))[c(7:16)])
sum(unlist(lapply(v3_single,length))[c(27)])
"ZNF227"   "CEACAM19" "CLU"      "APOC4"    "CASP2"    "NRG2"    
#for(k in length(v3_single[[27]])){
#for(j in (1:44)[-27]){
#  if(!v3_single[[27]][k]%in%v3_single[[j]]){
#    break
#  }
#}
#}
meta_single = unique(unlist(meta_single))
v3_single = unique(unlist(v3_single))
length(meta_joint); length(meta_single); length(intersect(meta_joint, meta_single))
length(v3_joint); length(v3_single); length(intersect(v3_joint, v3_single))
v3_joint[!v3_joint%in%v3_single]
[1] "RAB43"    "RAPSN"    "RPN1"     "CEACAM16" "DGKZ"     "FNBP4"    "HBEGF"    "MDK"  

[1] "ATXN7"      "PPP1R16B"   "PRRG2"      "RBMS3"      "RCN3"       "RPL3"       "RPS17"      "SATB2"      "SETD6"      "SFR1"      
[11] "SP4"        "SREBF1"     "BRD8"       "BTN1A1"     "TH"         "TMEM81"     "TMTC1"      "TNFRSF13C"  "TRAPPC3"    "TSKS"      
[21] "TSNAXIP1"   "TTC14"      "TWF2"       "UBXN2B"     "VPS29"      "ZDHHC5"     "ZNF703"     "C6orf15"    "ADAM10"     "CDC25C"    
[31] "CORO7"      "DBF4"       "DNPH1"      "FAM178B"    "FOXN2"      "FSD2"       "GATAD2A"    "GLT8D1"     "GPN3"       "HDDC3"     
[41] "HIRIP3"     "HIST1H2AJ"  "HIST1H2BG"  "HIST1H3C"   "HIST1H4B"   "HIST2H2AA3" "ITIH3"      "JKAMP"      "KCNG2"      "ARL14EP"   
[51] "MZB1"       "NES"        "NGF"        "PBRM1"      "PCDHA1"     "PCDHA4"     "PCDHA5"    
results_v3[results_v3[,2]%in%v3_joint[!v3_joint%in%v3_single], ][,2:4]
als_single = list()
fl = list.files(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/50traits/sig_gene/ALS/v3/local"))
setwd(paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/50traits/sig_gene/ALS/v3/local"))
for(i in 1:length(fl)){
  als_single[[i]] = readLines(fl[i])
}
als_single = unique(unlist(als_single))

7*6+15
