#initial content

rm(list= ls())

load.lib <- c('PerformanceAnalytics','xtable','stargazer', 'ggpubr','xts',"readxl","data.table","lubridate",'tidyr','rddensity','estimatr','rdrobust',
              'dplyr','ggplot2','forecast','RColorBrewer','quadprog','NMOF','fBonds','cvar','mFilter',"x13binary","seasonal","lmtest","RCurl") 
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib)
sapply(load.lib,require,character=TRUE)

url_robust <- "https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R"
eval(parse(text = getURL(url_robust, ssl.verifypeer = FALSE)),
     envir=.GlobalEnv
     
     options(scipen=10000)
     col<-brewer.pal(6,"Blues")
     
     dir <- 'C:/Users/57319/Documents/GitHub/RDD/Data'
     setwd(dir)
     
     datos=as.data.table(read.table("hansen_dwi.csv",header = T,sep=","))
     
     #3d point: Dummy variables only with BAC 0.08 
     
     datos$D <- as.numeric(datos$bac1 >= 0.08)
     datos$bac1 = datos$bac1
     
     datos1<- datos[bac1>=0.03 & bac1<=0.13]
     datos2<- datos[bac1>=0.055 & bac1<=0.105]
     
     #4 point: density test 
     
     MacCrary=rddensity(datos$bac1, c = 0.08,all=T)
     rdplotdensity(MacCrary, datos$bac1)
     
     
     histMC <-ggplot(datos,aes(x=datos$bac1))+
       geom_histogram(aes(y=..density..),color="darkorchid4", fill="darkorchid",bins = 100)+
       geom_density(col=col[6],size=2)+
       geom_vline(xintercept =c(0.08) ,
                  colour = c('red1'), linetype ="longdash", size = .6)+
       xlab('BAC') +ylab('Density')+
       theme(axis.title.x =element_text(size=13,  family="serif"),
             panel.background = element_blank(),
             panel.border = element_rect( fill=NA,colour = "black"),
             axis.text.x = element_text(angle = 90))
     histMC
     
     
     #Point 5: Checking the covariate balance 
     
     
     Balance_1 <- rdrobust(y = datos1$male,
                           x = datos1$bac1, c = 0.08,kernel = 'uniform',all = T)
     Balance_2 <- rdrobust(y = datos1$white,
                           x = datos1$bac1, c = 0.08,kernel = 'uniform',all = T)
     Balance_3 <- rdrobust(y = datos1$aged,
                           x = datos1$bac1, c = 0.08,kernel = 'uniform',all = T)
     Balance_4 <- rdrobust(y = datos1$acc,
                           x = datos1$bac1, c = 0.08,kernel = 'uniform',all= T)
     
     
     Balance_1_1 <- lm_robust(male ~ D + bac1 + D*bac1+ white+ aged+ acc, 
                              data = datos)
     Balance_2_1 <- lm_robust(white ~ D + bac1 + D*bac1+ male+ aged+ acc, 
                              data = datos)
     Balance_3_1 <- lm_robust(aged ~ D + bac1 + D*bac1+ white+ male +acc, 
                              data = datos)
     Balance_4_1 <- lm_robust(acc ~ D + bac1 + D*bac1+ white+ aged+ male, 
                              data = datos)
     
     Balance_1_2 <- lm(male ~ D + bac1 + D*bac1+ white+ aged+ acc, 
                       data = datos)
     Balance_2_2 <- lm(white ~ D + bac1 + D*bac1+ male+ aged+ acc, 
                       data = datos)
     Balance_3_2 <- lm(aged ~ D + bac1 + D*bac1+ white+ male +acc, 
                       data = datos)
     Balance_4_2 <- lm(acc ~ D + bac1 + D*bac1+ white+ aged+ male, 
                       data = datos)
     
     