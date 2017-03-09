rm(list = ls())
gc()
#path = file.choose()
setwd("/Users/jackbaska/Documents/JackBaska Interesting/Group/Hackthon/Git/hackthon/Water.gov/practice/Test/")
library(dplyr)

raw.air = read.csv("air.csv")
raw.rain = read.csv("rain.csv")
raw.water = read.csv("new_opendata_waterquality.csv")
##### clean raw.water
colnames(raw.water) = raw.water[1,] %>% as.matrix() %>% as.vector()
raw.water = raw.water[,-55]
raw.water = raw.water[2:nrow(raw.water),]
colnames(raw.water) = raw.water[2,] %>% as.matrix() %>% as.vector() %>%
  paste(colnames(raw.water),"(",.,")",sep="")
colnames(raw.water)[1:3] = c("區處別","系統代號","單位")
std.drink = raw.water[1,]
raw.water = raw.water[-c(1,2,(nrow(raw.water)-1),nrow(raw.water)),]
#################################################################
################### clean for analyzing data
#################################################################
data.water = raw.water %>% as.matrix() %>% as.data.frame()
data.analysis = data.water[,c(1,2,3,6,10)]
data.analysis = data.analysis %>% 
  sapply(function(x){gsub("=","",x)}) %>% 
  as.data.frame()
col.buf = data.analysis$單位 %>% as.character() %>% 
  strsplit("[()<>{}]|\\[|\\]") %>% lapply(function(x){
    x[which(x!="")]
  }) %>% do.call(rbind.data.frame,.)
colnames(col.buf) = c("淨水廠","淨水廠(英)","地址","代碼1","代碼2","縣市","縣市(英)")
col.buf = col.buf[c(1,3,6)]
data.buf = data.analysis
data.analysis[,(ncol(data.analysis)+1):(ncol(data.analysis)+ncol(col.buf))] = col.buf
data.analysis = data.analysis[,c(1:3,6,7,8,4,5)]
#################################################################
################### analyzing
#################################################################



##################################################################
##################################################################

# raw.rain$Rainfall24hr
# raw.rain %>% filter(Rainfall24hr!=0)
