####Traits paper analysis

library(reshape2)
library(vegan)
library(MASS) 
library(plyr)
library(ggfortify)
library(pls)
library(Hmisc)
#test3
#test

{#import the plant datafile
plants<- read.csv('data/sare_plants_traits.csv')
head(plants)
summary(plants)

#import other plant triat files to eventually merge
corolla<-read.csv('data/corolla_width_traits.csv')
head(corolla)

#merge corolla width data with plants file so that corolla width is represented for each plant
#note corolla width was only collected once at each site as an average of 5 samples.
plants<- merge(plants, corolla)
head(plants)
colnames(plants)

#also need to merge the color traits
traits<-read.csv('data/flw_col_traits.csv')
head(traits)

plants<- merge(plants, traits)
head(plants)
colnames(plants)

#also need to merge the color traits
nectar<-read.csv('data/nectar.csv')
head(nectar)

plants<- merge(plants, nectar)
#rid yourself of NAs
plants<-na.omit(plants)
summary(plants)
head(plants)
write.csv(plants, 'data/plant_traits.csv')
}#creating the traits data file. *do not need to use moving forward. just use plant_traits.csv*


plants<-read.csv('data/plant_traits.csv')

head(plants)

#create pollen per plot trait
plants$pol_plot<-plants$pol.unit*plants$tot_flw
plants<-na.omit(plants)
plants$per_tar_cov<-as.numeric(plants$per_tar_cov)

#create nectar per plot trait
plants$nec_plot<-plants$nectar*plants$tot_flw


#now it's time to look for correlations between plant trait data

{
  
  plants<-read.csv("data/plant_traits.csv")
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  plants$nec_plot<-plants$nectar*plants$tot_flw
  
#restructure data
library(plyr)
plants_cor<- ddply(plants, c("site","plant_sp"), summarise,
             Week_Bloom = mean(week_num),
             Number_Flowers = mean(tot_flw),
             avg_flw_area= mean(avg_flor_area),
            Floral_Area= mean(tot_flor_area),
             Flower_Height= mean(avg_tall_flw),
             avg_per_cov= mean(per_tar_cov),
             Corolla_Width=mean(avg_fw),
             Chroma=mean(chroma),
             Hue=mean(hue),
             #avg_pol_flw= mean(pol.unit),
             Pollen_Plot=mean(pol_plot),
            Nectar_Plot=mean(nec_plot))
head(plants_cor)
colnames(plants_cor)
plants_corr<-plants_cor[3:13]


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
#this prints the matrix 
library(Hmisc)
res2<-rcorr(as.matrix(plants_corr), type = c("pearson"))
flattenCorrMatrix(res2$r, res2$P)
#performance metrics
library("PerformanceAnalytics")
chart.Correlation(plants_corr, histogram=TRUE, pch=19)


}# plant trait correlation code
{
  plants<-read.csv('data/plant_traits.csv')
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  
  library(plyr)
  plants<- ddply(plants, c("site","plant_sp"), summarise,
                     Week_Bloom = mean(week_num),
                     Number_Flowers = mean(tot_flw),
                     avg_flw_area= mean(avg_flor_area),
                     Floral_Area= mean(tot_flor_area),
                     Flower_Height= mean(avg_tall_flw),
                     avg_per_cov= mean(per_tar_cov),
                     Corolla_Width=mean(avg_fw),
                     Chroma=mean(chroma),
                     Hue=mean(hue),
                     #avg_pol_flw= mean(pol.unit),
                     Pollen_Plot=mean(pol_plot))
  head(plants)
  
  #number of flowers vs floral area
  #week bloom vs flower height
  #floral area vs flower height
  #number flowers vs corolla width
  #floral area vs corolla width
  #flower height vs corolla width
  #floral area vs hue
  #flower height vs hue
  
}#simple linear regression models for relating plant traits

####OKAY so based off of the correlation information we have that 
#Floral area and Pollen per plot are correlated r2=0.53
#Chroma and Hue are Correlated r2=-0.5


# to avoid issues with autocorrelation, I removed Pollen and Hue variables from MLR

{
###combine plant data with bee data
#import the bee dataset
specimen <- read.csv('data/all_bees3.csv') 
head(specimen)
colnames(specimen) #make sure it looks correct

#reshape bee file to merge with plants_cor file
bees_sum<-dcast(specimen,plant_sp ~family, length)   
head(bees_sum) 
colnames(bees_sum)

#combine these two data files
plants_bees <- merge(plants_cor, bees_sum)
head(plants_bees)
}#combining bee data with plant data to determine best models

{
  #Apidae
  hist(plants_bees$Apidae)
  plants_bees$logApidae<-log(plants_bees$Apidae +1)
  hist(plants_bees$logApidae)
  
  #Andrenidae
  hist(plants_bees$Andrenidae)
  plants_bees$logAndrenidae<-log(plants_bees$Andrenidae +1)
  hist(plants_bees$logAndrenidae)
  
  #Colletidae
  hist(plants_bees$Colletidae)
  plants_bees$logColletidae<-log(plants_bees$Colletidae +1)
  hist(plants_bees$logColletidae)
  
  #Megachilidae
  hist(plants_bees$Megachilidae)
  plants_bees$logMegachilidae<-log(plants_bees$Megachilidae +1)
  hist(plants_bees$logMegachilidae)
  
  #Halictidae
  hist(plants_bees$Halictidae)
  plants_bees$logHalictidae<-log(plants_bees$Halictidae +1)
  hist(plants_bees$logHalictidae)
  
  
}# log transformations of wild bee abundance data

{
#Apidae model
model.Ap<- lm(logApidae~ Week_Bloom + Number_Flowers + Floral_Area + Flower_Height + Corolla_Width + Hue,
                    data = plants_bees)
  
  summary(model.Ap)
  #ANOVA F tests
  modela<-aov(model.Ap)
  summary(modela)
  
  #Using the dredge function to determine best model
  library(MuMIn)
  options(na.action = "na.fail")
  pos.model = dredge(model.Ap)
  best.mods = subset(pos.model, delta < 3)
  best.mods
  avg.mod = model.avg(best.mods)
  summary(avg.mod) 
  
#Andrenidae model
model.An<- lm(logAndrenidae~ Week_Bloom + Number_Flowers + Floral_Area + Flower_Height + Corolla_Width + Hue,
           data = plants_bees)

  summary(model.An)
  modela<-aov(model.An)
  summary(modela)

  #Using the dredge function to determine best model
  pos.model = dredge(model.An)
  best.mods = subset(pos.model, delta < 3)
  best.mods
  avg.mod = model.avg(best.mods)
  summary(avg.mod)


#Colletidae model
model.Co<- lm(logColletidae~ Week_Bloom + Number_Flowers + Floral_Area + Flower_Height + Corolla_Width + Hue,
           data = plants_bees)

  summary(model.Co)
  modela<-aov(model.Co)
  summary(modela)

  #Using the dredge function to determine best model
  pos.model = dredge(model.Co)
  best.mods = subset(pos.model, delta < 3)
  best.mods
  avg.mod = model.avg(best.mods)
  summary(avg.mod)


#Halictidae model
model.Ha<- lm(logHalictidae~ Week_Bloom + Number_Flowers + Floral_Area + Flower_Height + Corolla_Width + Hue,
           data = plants_bees)

  summary(model.Ha)
  modela<-aov(model.Ha)
  summary(modela)

  #Using the dredge function to determine best model
  pos.model = dredge(model.Ha)
  best.mods = subset(pos.model, delta < 3)
  best.mods
  avg.mod = model.avg(best.mods)
  summary(avg.mod)

#Megachilidae model
model.Me<- lm(logMegachilidae~ Week_Bloom + Number_Flowers + Floral_Area + Flower_Height + Corolla_Width + Hue,
           data = plants_bees)

  summary(model.Me)
  modela<-aov(model.Me)
  summary(modela)

  #Using the dredge function to determine best model
  pos.model = dredge(model.Me)
  best.mods = subset(pos.model, delta < 3)
  best.mods
  avg.mod = model.avg(best.mods)
  summary(avg.mod)


}#using dredging functions to determine best models and an average best model

{# now lets do some PCA because PCA can utilize trait variables that are correlated.
  #Multivariate analysis at the family level- pollen excluded from analysis
  plants<-read.csv('data/plant_traits.csv')
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants$nec_plot<-plants$nectar*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  head(plants)
  
  plants_t<-ddply(plants, c("plant_sp", "site"), summarise,
                  Week_Bloom = mean(week_num),
                  Number_Flowers = mean(tot_flw),
                  #avg_flw_area= mean(avg_flor_area),
                  Floral_Area= mean(tot_flor_area),
                  Flower_Height= mean(avg_tall_flw),
                  Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot))
  
  plants_t #this the plants trait file
  head(plants_t)
  
  #import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #this creates the specimen matrix with site, year, bloom period attached
  bee_genus<-dcast(specimen,plant_sp+site ~family, length)   
  head(bee_genus) 
  colnames(bee_genus)
  
  #merge the plant_t and specimen data files
  plants_bees <- merge(plants_t, bee_genus)
  head(plants_bees)
  
  #remove NAs because NMDS can't hang
  plants_bees<- na.omit(plants_bees)
  head(plants_bees)
  colnames(plants_bees)
  
  #look at genus sums- will remove genera with very low numbers
  colSums(plants_bees[12:16])
  
  #remove genera with low numbers of bees 
  #plants_bees$Anthidiellum<-NULL #0 total
  colnames(plants_bees)
  
  
  tot_bees<-rowSums(plants_bees[12:16])
  plants_bees_condensed<-plants_bees[tot_bees>1,] #removing interactions zeros
  colnames(plants_bees_condensed)
  
  
  #creating bee and plant matrix data for PCA analysis
  bee_matrix<- plants_bees_condensed[12:16]
  colnames(bee_matrix)
  head(bee_matrix)
  
  plant_data<- plants_bees_condensed[1:11]
  colnames(plant_data)
  
  #remove unneccesary information
  plant_data$plant_sp<-NULL
  plant_data$site<-NULL
  plant_data$year<-NULL
  plant_data$block<-NULL
  plant_data$Pollen_Plot<-NULL
  
  ##PCA #I like this the most
  ord<-rda(bee_matrix, distance="bray", autotransform=T, scale=T)
  ord
  plot(ord, display="species", type= "t")
  
  ef <- envfit(ord, plant_data, permu = 1000)
  ef
  plot(ef)
}##PCA analysis family level does not include pollen data

{
  
  #Multivariate analysis at the genus level- pollen excluded from analysis
plants<-read.csv('data/plant_traits.csv')
head(plants)
plants$pol_plot<-plants$pol.unit*plants$tot_flw
plants$nec_plot<-plants$nectar*plants$tot_flw
plants<-na.omit(plants)
plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
head(plants)

plants_t<-ddply(plants, c("plant_sp", "site"), summarise,
                Week_Bloom = mean(week_num),
                Number_Flowers = mean(tot_flw),
                #avg_flw_area= mean(avg_flor_area),
                Floral_Area= mean(tot_flor_area),
                Flower_Height= mean(avg_tall_flw),
                Percent_Tar_Cov= mean(per_tar_cov),
                Corolla_Width=mean(avg_fw),
                Chroma=mean(chroma),
                Hue=mean(hue),
                #avg_pol_flw= mean(pol.unit),
                Pollen_Plot=mean(pol_plot))

plants_t #this the plants trait file
head(plants_t)

#import specimen dataset
specimen <- read.csv('data/all_bees3.csv') 
head(specimen)
colnames(specimen) #make sure it looks correct

#this creates the specimen matrix with site, year, bloom period attached
bee_genus<-dcast(specimen,plant_sp+site ~genus, length)   
head(bee_genus) 
colnames(bee_genus)

#merge the plant_t and specimen data files
plants_bees <- merge(plants_t, bee_genus)
head(plants_bees)

#remove NAs because NMDS can't hang
plants_bees<- na.omit(plants_bees)
head(plants_bees)
colnames(plants_bees)

#look at genus sums- will remove genera with very low numbers
colSums(plants_bees[12:42])

#remove genera with low numbers of bees 
plants_bees$Anthidiellum<-NULL #0 total
colnames(plants_bees)


tot_bees<-rowSums(plants_bees[12:41])
plants_bees_condensed<-plants_bees[tot_bees>1,] #removing interactions zeros
colnames(plants_bees_condensed)


#creating bee and plant matrix data for PCA analysis
bee_matrix<- plants_bees_condensed[12:41]
colnames(bee_matrix)
head(bee_matrix)

plant_data<- plants_bees_condensed[1:11]
colnames(plant_data)

#remove unneccesary information
plant_data$plant_sp<-NULL
plant_data$site<-NULL
plant_data$year<-NULL
plant_data$block<-NULL
plant_data$Pollen_Plot<-NULL

##PCA #I like this the most
ord<-rda(bee_matrix, distance="bray", autotransform=T, scale=T)
ord
plot(ord, display="species", type= "t")

ef <- envfit(ord, plant_data, permu = 1000)
ef
plot(ef)
}##PCA analysis genus level but does not include any pollen data

{
  #Multivariate analysis at the genus level- pollen included in this analysis
  plants<-read.csv('data/plant_traits.csv')
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants$nec_plot<-plants$nectar*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  head(plants)
  
  plants_t<-ddply(plants, c("plant_sp", "site"), summarise,
                  Week_Bloom = mean(week_num),
                  Number_Flowers = mean(tot_flw),
                  #avg_flw_area= mean(avg_flor_area),
                  Floral_Area= mean(tot_flor_area),
                  Flower_Height= mean(avg_tall_flw),
                  Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot))
  
  plants_t #this the plants trait file
  head(plants_t)
  
  #import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #this creates the specimen matrix with site, year, bloom period attached
  bee_genus<-dcast(specimen,plant_sp+site ~genus, length)   
  head(bee_genus) 
  colnames(bee_genus)
  
  #merge the plant_t and specimen data files
  plants_bees <- merge(plants_t, bee_genus)
  head(plants_bees)
  
  #remove NAs because NMDS can't hang
  plants_bees<- na.omit(plants_bees)
  head(plants_bees)
  colnames(plants_bees)
  
  #look at genus sums- will remove genera with very low numbers
  colSums(plants_bees[12:42])
  
  #remove genera with low numbers of bees 
  plants_bees$Anthidiellum<-NULL #0 total
  colnames(plants_bees)
  
  
  tot_bees<-rowSums(plants_bees[12:41])
  plants_bees_condensed<-plants_bees[tot_bees>1,] #removing interactions zeros
  colnames(plants_bees_condensed)
  
  
  #creating bee and plant matrix data for PCA analysis
  bee_matrix<- plants_bees_condensed[12:41]
  colnames(bee_matrix)
  head(bee_matrix)
  
  plant_data<- plants_bees_condensed[1:11]
  colnames(plant_data)
  
  #remove unneccesary information
  plant_data$plant_sp<-NULL
  plant_data$site<-NULL
  plant_data$year<-NULL
  plant_data$block<-NULL
  #plant_data$Pollen_Plot<-NULL
  
  ##PCA #I like this the most
  ord<-rda(bee_matrix, distance="bray", autotransform=T, scale=T)
  ord
  plot(ord, display="species", type= "t")
  
  ef <- envfit(ord, plant_data, permu = 1000)
  ef
  plot(ef)
 }##PCA analysis does include pollen data

{
  plants<- read.csv('data/sare_plants_traits2.csv')
  head(plants)
  summary(plants) 
  
  #import other plant trait files to eventually merge
  corolla<-read.csv('data/corolla_width_traits.csv')
  head(corolla)
  
  #merge corolla width data with plants file so that corolla width is represented for each plant
  #note corolla width was only collected once at each site as an average of 5 samples.
  plants<- merge(plants, corolla)
  head(plants)
  colnames(plants)
  
  #also need to merge the color traits
  traits<-read.csv('data/flw_col_traits.csv')
  head(traits)
  
  plants<- merge(plants, traits)
  head(plants)
  colnames(plants)
  
  plants<-na.omit(plants) 
  
  #restructure this data to merge with bee data
  library(plyr)
  plants_mlr<- ddply(plants, c("site", "plant_sp"), summarise,
                     Week_Bloom = mean(week_num),
                     Number_Flowers = mean(tot_flw),
                     avg_flw_area= mean(avg_flor_area),
                     Floral_Area= mean(tot_flor_area),
                     Flower_Height= mean(avg_tall_flw),
                     #avg_per_cov= mean(per_tar_cov),
                     Corolla_Width=mean(avg_fw),
                     Chroma=mean(chroma),
                     Honey_Bees=mean(tot_hb),
                     Syrphid=mean(tol_syr),
                     Hue=mean(hue),
                     avg_pol_flw= mean(pol.unit))
                     #Pollen_Plot=mean(pol_plot))
  head(plants_mlr)
  colnames(plants_mlr)
  
  #bring in bee data
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #this creates the specimen matrix with site, year, bloom period attached
  bee_genus<-dcast(specimen,plant_sp+site ~family, length)   
  head(bee_genus) 
  colnames(bee_genus)
  
  #merge the plant_t and specimen data files
  plants_bees <- merge(plants_mlr, bee_genus)
  head(plants_bees)
  colnames(plants_bees)
  
  #creating a sum collected coloumn
  plants_bees$sumcol<-plants_bees$Andrenidae+plants_bees$Apidae+
                              plants_bees$Colletidae+ plants_bees$Halictidae
                              +plants_bees$Megachilidae
  
  
  #LOG TRANSFORMATIONS OF THE ABUNDANCE DATA
  #Apidae
  hist(plants_bees$Apidae)
  plants_bees$logApidae<-log(plants_bees$Apidae +1)
  hist(plants_bees$logApidae)
  
  #Andrenidae
  hist(plants_bees$Andrenidae)
  plants_bees$logAndrenidae<-log(plants_bees$Andrenidae +1)
  hist(plants_bees$logAndrenidae)
  
  #Colletidae
  hist(plants_bees$Colletidae)
  plants_bees$logColletidae<-log(plants_bees$Colletidae +1)
  hist(plants_bees$logColletidae)
  
  #Megachilidae
  hist(plants_bees$Megachilidae)
  plants_bees$logMegachilidae<-log(plants_bees$Megachilidae +1)
  hist(plants_bees$logMegachilidae)
  
  #Halictidae
  hist(plants_bees$Halictidae)
  plants_bees$logHalictidae<-log(plants_bees$Halictidae +1)
  hist(plants_bees$logHalictidae)
  
  #Honey bees
  hist(plants_bees$Honey_Bees)
  plants_bees$loghoneybees<-log(plants_bees$Honey_Bees +1)
  hist(plants_bees$loghoneybees)
  
  #Syrpids
  hist(plants_bees$Syrphid)
  plants_bees$logsyrpids<-log(plants_bees$Syrphid +1)
  hist(plants_bees$logsyrpids)
  
  #Collected
  hist(plants_bees$sumcol)
  plants_bees$logsumcol<-log(plants_bees$sumcol +1)
  hist(plants_bees$logsumcol)
  
  #MULTIPLE LINEAR REGRESSION ANALYSIS
  library(car)
  #honeybees
  honeybees<- lm ( loghoneybees ~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(honeybees)
  anova(honeybees, test="F")
  #plot(honey_bees)
  vif(honeybees)
  
  #plots to assess individual relationships
  reg1<-lm(Week_Bloom~loghoneybees, data=plants_bees)
  summary(reg1)
library(ggplot2)  
  
 p1<- ggplot(plants_bees, aes(x=Week_Bloom, y=loghoneybees)) + geom_point()  + stat_smooth(method=lm, level= .99, color="black", se=FALSE) 
 p1 <- p1 + ylab("Honey bees per sample (logx + 1)")
 p1 <-p1 + xlab("Week of Bloom")
 p1 <- p1 + annotate("text",label="r^2 == .12", parse= TRUE, x=39.5, y=1)
  #hb_bb <- hb_bb + ggtitle("a.   2015")
 p1 <- p1 + theme_bw()
 p1 <- p1 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 p1 <- p1 + theme(text= element_text(size = 10))
 p1
  
  
  
  #syrphids
  syrphids<- lm ( logsyrpids ~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(syrphids)
  anova(syrphids, test="F")
  #plot(honey_bees)
  vif(syrphids) #should be the same as for honey bees, same predictors used
  
  #wild bees
  wildbees<- lm ( logsumcol ~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(wildbees)
  anova(wildbees, test="F")
  #plot(honey_bees)
  vif(wildbees)
  
  #apidae
  Apidae<- lm ( logApidae ~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(Apidae)
  anova(Apidae, test="F")
  #plot(honey_bees)
  vif(Apidae)
  
  #andrenidae
  Andrenidae<- lm ( logAndrenidae ~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(Andrenidae)
  anova(Andrenidae, test="F")
  #plot(honey_bees)
  vif(Andrenidae)
  
  #Colletidae
  Colletidae<- lm ( logColletidae ~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(Colletidae)
  anova(Colletidae, test="F")
  #plot(honey_bees)
  vif(Colletidae)
  
  
  #Halictidae
  Halictidae<- lm ( logHalictidae~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(Halictidae)
  anova(Halictidae, test="F")
  #plot(honey_bees)
  vif(Halictidae)
  
  #Megachilidae
  Megachilidae<- lm ( logMegachilidae~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(Megachilidae)
  anova(Megachilidae, test="F")
  #plot(honey_bees)
  vif(Megachilidae)
  
  #Bumblebees only
  #this code will allow us to look at bumblebees specifically
  bee_genus<-dcast(specimen,plant_sp+site ~genus, length)   
  head(bee_genus) 
  colnames(bee_genus)
  
  #merge the plant_t and specimen data files
  plants_bees <- merge(plants_mlr, bee_genus)
  head(plants_bees)
  colnames(plants_bees)
  
  #Bumblebees only
  hist(plants_bees$Bombus)
  plants_bees$logbombus<-log(plants_bees$Bombus +1)
  hist(plants_bees$logbombus)
  
  Bombus<- lm ( logbombus~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(Bombus)
  anova(Bombus, test="F")
  vif(Megachilidae)
  
  #Non-Bombus Apidae
  ###need to create a column that includes Apidae minus Bombus for this analysis
  #family first
  bee_family<-dcast(specimen,plant_sp+site ~family, length)   
  head(bee_family) 
  colnames(bee_family)
  
  #genus second
  bee_genus<-dcast(specimen,plant_sp+site ~genus, length)   
  head(bee_genus) 
  colnames(bee_genus)
  
  #combine family with Genus
  bees <- merge(bee_family, bee_genus)
  head(bees)
  colnames(bees)
  
  #merge this with the plant file
  plants_bees <- merge(plants_mlr, bees)
  head(plants_bees)
  colnames(plants_bees)
  
  #create non-bombus apidae column
  plants_bees$nonbomusapidae<- plants_bees$Apidae-plants_bees$Bombus
  
  
  #non-bombus apidae only
  hist(plants_bees$nonbomusapidae)
  plants_bees$lognonbomus<-log(plants_bees$nonbomusapidae +1)
  hist(plants_bees$lognonbomus)
  
  non_Bombus_apidae<- lm ( lognonbomus~ Week_Bloom+ Floral_Area +Number_Flowers + Flower_Height + Corolla_Width + Hue, data=plants_bees)
  summary(non_Bombus_apidae)
  anova(non_Bombus_apidae, test="F")
  vif(non_Bombus_apidae)
  
}#Multiple linear regressions to see if traits explain visitation



