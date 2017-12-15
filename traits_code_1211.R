####Traits paper analysis

library(reshape2)
library(vegan)
library(MASS) 
library(plyr)
library(ggfortify)
library(pls)
library(Hmisc)
library(car)
library(ggplot2)
library(gridExtra)


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


plants<-read.csv('plant_traits.csv')

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
  plants_cor<- ddply(plants, c("site", "plant_sp"), summarise,
                     Week_Bloom = mean(week_num),
                     Number_Flowers = mean(tot_flw),
                     #avg_flw_area= mean(avg_flor_area),
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
  plants_corr<-plants_cor[3:11]
  
  
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
  model.Ap<- lm(logApidae~ Week_Bloom + Number_Flowers + Floral_Area + Flower_Height + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot,
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

{
  plants<-read.csv("data/plant_traits.csv")
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
                  #Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot),
                  Nectar_Plot=mean(nec_plot))
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
  #plant_data$Pollen_Plot<-NULL
  
  ##PCA #I like this the most
  ord<-rda(plants, distance="bray", autotransform=T, scale=T)
  ord
  plot(ord, display="species", type= "t")
 summary(ord)
  ef <- envfit(ord, plant_data, permu = 1000)
  ef
  plot(ef)
  summary(ef)
  
  fit <- princomp(bee_matrix, cor=TRUE)
  summary(fit) # print variance accounted for 
  loadings(fit) # pc loadings 
}##PCA analysis Family level does include pollen data

{# now lets do some PCA because PCA can utilize trait variables that are correlated.
  plants<-read.csv("data/plant_traits.csv")
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
                  #Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot),
                  Nectar_Plot=mean(nec_plot))
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
  #plants_bees$Anthidiellum<-NULL #0 total
  colnames(plants_bees)
  
  plants_bees$Anthidiellum<-NULL #0 total
  colnames(plants_bees)
  
  
  
  tot_bees<-rowSums(plants_bees[12:41])
  plants_bees_condensed<-plants_bees[tot_bees>5,] #removing interactions zeros
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
}##PCA analysis genus level but does  pollen data

{
  plants<-read.csv('data/plant_traits.csv')
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants$nec_plot<-plants$nectar*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  head(plants)
  
  plants_t<-ddply(plants, c("plant_sp"), summarise,
                  Week_Bloom = mean(week_num),
                  Number_Flowers = mean(tot_flw),
                  #avg_flw_area= mean(avg_flor_area),
                  Floral_Area= mean(tot_flor_area),
                  Flower_Height= mean(avg_tall_flw),
                  #Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot),
                  Nectar_Plot=mean(nec_plot))
  
  plants_t #this the plants trait file
  head(plants_t)
  
  #import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #this creates the specimen matrix with site, year, bloom period attached
  bee_genus<-dcast(specimen,plant_sp ~genus, length)   
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
  colSums(plants_bees[11:41])
  
  #remove genera with low numbers of bees 
  plants_bees$Anthidiellum<-NULL #0 total
  colnames(plants_bees)
  
  
  tot_bees<-rowSums(plants_bees[11:40])
  plants_bees_condensed<-plants_bees[tot_bees>1,] #removing interactions zeros
  colnames(plants_bees_condensed)
  
  
  #creating bee and plant matrix data for PCA analysis
  bee_matrix<- plants_bees_condensed[11:40]
  colnames(bee_matrix)
  head(bee_matrix)
  
  plant_data<- plants_bees_condensed[1:10]
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
  
  plants<-read.csv("data/plant_traits.csv")
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  plants$nec_plot<-plants$nectar*plants$tot_flw
  
  #restructure data
  library(plyr)
  plants_cor<- ddply(plants, c("plant_sp"), summarise,
                     Week_Bloom = mean(week_num),
                     Number_Flowers = mean(tot_flw),
                     #avg_flw_area= mean(avg_flor_area),
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
  
  ###combine plant data with bee data
  #import the bee dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #adding species ricnhess
  bees_sr<-dcast(specimen,plant_sp ~species, length)   
  head(bees_sr) 
  colnames(bees_sr)
  
  bees<-bees_sr[2:136]
  head(bees)
  
  library(vegan)
  SR<-specnumber(bees)
  SR
  
  
  
  
  
  #reshape bee file to merge with plants_cor file
  bees_sum<-dcast(specimen,site+year+week_num+block+plant_sp ~family, length)   
  head(bees_sum) 
  colnames(bees_sum)
  
  
  bee_sum<- ddply(bees_sum, c("plant_sp"), summarise,
                  Andrenidae= mean(Andrenidae),
                  Apidae = mean(Apidae),
                  #avg_flw_area= mean(avg_flor_area),
                  Colletidae= mean(Colletidae),
                  Halictidae= mean(Halictidae),
                  Megachilidae= mean(Megachilidae))
  
  head(bee_sum)
  colnames(bee_sum)
  
  #combine species richness file with bee_sum file
  bee_sum<-cbind(bee_sum,SR)
  colnames(bee_sum)
  #combine these two data files
  plants_bees <- merge(plants_cor, bee_sum)
  head(plants_bees)
  
  
  plants2<- read.csv('data/sare_plants_traits2.csv')
  head(plants2)
  
  #reshape to only have poll data combined with plant data
  library(plyr)
  plants2<- ddply(plants2, c("plant_sp"), summarise,
                  Honey_Bees=mean(tot_hb),
                  Syrphid=mean(tol_syr))
  head(plants2)
  colnames(plants2)
  
  #combine these two data files
  plants_bees <- merge(plants_bees, plants2)
  head(plants_bees)
  
  
  #creating a sum collected coloumn
  plants_bees$sumcolwb<-plants_bees$Andrenidae+plants_bees$Apidae+
    plants_bees$Colletidae+ plants_bees$Halictidae +plants_bees$Megachilidae
  
  head(plants_bees)
  
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
  hist(plants_bees$sumcolwb)
  plants_bees$logsumcolwb<-log(plants_bees$sumcolwb +1)
  hist(plants_bees$logsumcolwb)
  
  
  colnames(plants_bees)
  
  plants_bees<-na.omit(plants_bees)
  plants_bees$avg_per_cov<-NULL
  plants_bees$Nectar_Plot<-NULL
  
  #plants_bees$avg_per_cov<-NULL
  #plants_bees$Chroma<-NULL
  
  #MULTIPLE LINEAR REGRESSION ANALYSIS
  library(rsq)
  
  #this first is to determine variance inflation factors to determine which factors
  #can be used in the model selection processes that follow
  overal_model<-lm(loghoneybees~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                   + Corolla_Width + Chroma + Hue + Pollen_Plot, data=plants_bees)
  
  vif(overal_model)#looks like all factors can be used in the model selection
  
  
  
  {
    
    {  
      
      #honeybees
      #Forward selection method
      #first I need to remove unecessary variables
      colnames(plants_bees)
      plants_hb<-plants_bees[2:9]
      plants_hb2<-plants_bees[24:24]
      loghb<-cbind(plants_hb2,plants_hb)
      colnames((loghb))
      
      {#Dredge first- Model Selection
        honeybees<-lm(loghoneybees~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                      + Corolla_Width + Chroma + Hue + Pollen_Plot, data=loghb)
        
        library(MuMIn)
        options(na.action = "na.fail")
        pos.model = dredge(honeybees)
        best.mods = subset(pos.model, delta < 3)
        best.mods
        
        
        #honeybees
        honeybees<- lm ( loghoneybees ~Flower_Height, data=loghb)
        summary(honeybees)
        anova(honeybees, test="F")
        rsq.partial(honeybees)
        rsq(honeybees)
      } #Dredge
      
      #loghb$Pollen_Plot<-NULL
      #loghb$Nectar_Plot<-NULL
      
      {
        
        honeybees<-lm(loghoneybees~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                      + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=loghb)
        
        ####use
        Fitall<- lm(loghoneybees~., data=loghb)
        Fitstart<-lm(loghoneybees~1, data=loghb)
        step(Fitstart, direction = "forward", scope=formula(Fitall))
        
        #honeybees
        honeybees<- lm ( loghoneybees ~ Floral_Area+Nectar_Plot, data=loghb)
        summary(honeybees)
        anova(honeybees, test="F")
        rsq.partial(honeybees)
        rsq(honeybees)
        #plot(honey_bees)
        AIC(honeybees)
        
        
      }#Forward Selection

      
    }#models  
    {#plots to assess individual relationships
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
      
    }#extra graphics
  }#honeybee model selection and multiple linear regression
  {  
    
    #syrphid flies
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_syr<-plants_bees[2:9]
    plants_syr2<-plants_bees[25:25]
    logsyr<-cbind(plants_syr2,plants_syr)
    colnames((logsyr))
    
    {#Dredge first- Model Selection
      hoverflies<-lm(logsyrpids~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                     + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logsyr)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(hoverflies)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      
      #syrphids
      logsyrpids<- lm ( logsyrpids ~ Floral_Area+Number_Flowers+Week_Bloom, data=logsyr)
      summary(logsyrpids)
      anova(logsyrpids, test="F")
      rsq.partial(logsyrpids)
      rsq(logsyrpids)
    } #Dredge
    
   
  
  }#syrphid flies model
  {
    
    {#wild bees collected
      #Forward selection method
      #first I need to remove unecessary variables
      #first with floral area and not pollen
      colnames(plants_bees)
      plants_wb<-plants_bees[2:10]
      plants_wb2<-plants_bees[26:26]
      logwb<-cbind(plants_wb2,plants_wb)
      colnames((logwb))
      logwb$Pollen_Plot<-NULL
      #loghb$Nectar_Plot<-NULL
      
      
      
      Fitall<- lm(logsumcolwb~., data=logwb)
      Fitstart<-lm(logsumcolwb~1, data=logwb)
      step(Fitstart, direction = "forward", scope=formula(Fitall))
      
      #wild bees
      wb<- lm ( logsumcolwb ~ Floral_Area + avg_per_cov + Corolla_Width + 
                  Hue + Chroma, data = logwb)
      summary(wb)
      anova(wb, test="F")
      rsq.partial(wb)
      rsq(wb)
      AIC(wb)
      
      #switch floral area with pollen
      colnames(plants_bees)
      plants_wb<-plants_bees[3:12]
      plants_wb2<-plants_bees[28:28]
      logwb<-cbind(plants_wb2,plants_wb)
      colnames((logwb))
      #logwb$Pollen_Plot<-NULL
      logwb$Floral_Area<-NULL
      
      
      
      Fitall<- lm(logsumcolwb~., data=logwb)
      Fitstart<-lm(logsumcolwb~1, data=logwb)
      step(Fitstart, direction = "forward", scope=formula(Fitall))
      
      #wild bees
      wb<- lm ( formula = logsumcolwb ~ avg_per_cov + Pollen_Plot + Corolla_Width + 
                  Hue, data = logwb)
      summary(wb)
      anova(wb, test="F")
      rsq.partial(wb)
      rsq(wb)
      AIC(wb)
      
      
    }#with bombus
    {
      #wild bees collected
      #Forward selection method
      #first I need to remove unecessary variables
      #first with floral area and not pollen
      colnames(plants_bees)
      
      plants_bees$wbnobombus<-plants_bees$sumcolwb-plants_bees$Bombus
      
      #Collected
      hist(plants_bees$wbnobombus)
      plants_bees$logwbnobombus<-log(plants_bees$wbnobombus +1)
      hist(plants_bees$logwbnobombus)
      
      colnames(plants_bees)
      plants_wb<-plants_bees[2:9]
      plants_wb2<-plants_bees[30:30]
      logwb<-cbind(plants_wb2,plants_wb)
      colnames((logwb))
      
      
      {#Dredge first- Model Selection
        nonbombuswild<-lm(logwbnobombus~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                          + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
        
        library(MuMIn)
        options(na.action = "na.fail")
        pos.model = dredge(nonbombuswild)
        best.mods = subset(pos.model, delta < 3)
        best.mods
        
        colnames(logwb)
        #wild bees no bombus
        nonbombuswild<- lm(logwbnobombus ~ Floral_Area + Week_Bloom, data=logwb)
        summary(nonbombuswild)
        anova(nonbombuswild, test="F")
        rsq.partial(nonbombuswild)
        rsq(nonbombuswild)
      } #Dredge
      
      #logwb$Pollen_Plot<-NULL
      #loghb$Nectar_Plot<-NULL
      
      
      {
        
        lognonbombuswild<-lm(logwbnobombus~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                             + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=logwb)
        
        Fitall<- lm(logwbnobombus~., data=logwb)
        Fitstart<-lm(logwbnobombus~1, data=logwb)
        step(Fitstart, direction = "forward", scope=formula(Fitall))
        
        #wild bees
        wb<- lm (logwbnobombus ~ Floral_Area, data = logwb)
        summary(wb)
        anova(wb, test="F")
        rsq.partial(wb)
        rsq(wb)
        AIC(wb)
      }#forward selection
      
    
    }#without bombus
    

    
  }#wild bees total model
  { 
    
    #Apidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[3:11]
    plants_wb2<-plants_bees[21:21]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    
    {#Dredge first- Model Selection
      logApidae<-lm(logApidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                    + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(logApidae)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      logAndrenidae<- lm( logApidae ~ Corolla_Width+ Floral_Area+Flower_Height+Hue, data=logwb)
      summary(logAndrenidae)
      anova(logAndrenidae, test="F")
      rsq.partial(logAndrenidae)
      rsq(logAndrenidae)
    } #Dredge
    
    
  
    
    
    {  
    Fitall<- lm(logApidae~., data=logwb)
    Fitstart<-lm(logApidae~1, data=logwb)
    step(Fitstart, direction = "forward", scope=formula(Fitall))
    
    #Apidae only
    wb<- lm ( logApidae ~ Floral_Area + Flower_Height + Hue + 
                avg_per_cov + Number_Flowers + Corolla_Width, data = logwb)
    summary(wb)
    anova(wb, test="F")
    rsq.partial(wb)
    rsq(wb)
    AIC(wb)
    
    }#forward selection
  
    
    ###BETTER MODEL IS MODEL WITH FLORAL AREA 
    
    
    
  }#Apidae model
  {
    #Andrenidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[20:20]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    
    {#Dredge first- Model Selection
      logAndrenidae<-lm(logAndrenidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                        + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(logAndrenidae)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      logAndrenidae<- lm( logAndrenidae ~Corolla_Width+Pollen_Plot, data=logwb)
      summary(logAndrenidae)
      anova(logAndrenidae, test="F")
      rsq.partial(logAndrenidae)
      rsq(logAndrenidae)
    } #Dredge
    
    
    #logwb$Pollen_Plot<-NULL
    #loghb$Nectar_Plot<-NULL
    
    {
      logAndrenidae<-lm(logAndrenidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                        + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      Fitall<- lm(logAndrenidae~., data=logwb)
      Fitstart<-lm(logAndrenidae~1, data=logwb)
      step(Fitstart, direction = "forward", scope=formula(Fitall))
      
      #Andrenidae only
      wb<- lm (formula = logAndrenidae ~ Pollen_Plot+ Floral_Area, data = logwb)
      summary(wb)
      anova(wb, test="F")
      rsq.partial(wb)
      rsq(wb)
      AIC(wb)
      
    }#forward selection
  
    
  }#Andrenidae model
  { 
    
    #Colletidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[21:21]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    
    
    {#Dredge first- Model Selection
      logColletidae<-lm(logColletidae~Pollen_Plot + Number_Flowers + Floral_Area + Flower_Height
                        + Corolla_Width + Chroma + Hue + Week_Bloom, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(logColletidae)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      logColletidae<- lm( logColletidae ~ Corolla_Width, data=logwb)
      summary(logColletidae)
      anova(logColletidae, test="F")
      rsq.partial(logColletidae)
      rsq(logColletidae)
    } #Dredge
    
    
    #logwb$Pollen_Plot<-NULL
    #loghb$Nectar_Plot<-NULL
    
    {
      logColletidae<-lm(logColletidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                        + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      
      Fitall<- lm(logColletidae~., data=logwb)
      Fitstart<-lm(logColletidae~1, data=logwb)
      step(Fitstart, direction = "forward", scope=formula(Fitall))
      
      #Colletidae only
      wb<- lm (logColletidae ~ Corolla_Width + Pollen_Plot + Floral_Area, data = logwb)
      summary(wb)
      anova(wb, test="F")
      rsq.partial(wb)
      rsq(wb)
      AIC(wb)
      
    }#forward selection
  
    
    #Models with Floral area and Pollen performed the same
    
  }#Colletidae model
  {
    #Halictidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[23:23]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    {#Dredge first- Model Selection
      logHalictidae<-lm(logHalictidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                        + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(logHalictidae)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      logHalictidae<- lm( logHalictidae ~ Floral_Area+Week_Bloom, data=logwb)
      summary(logHalictidae)
      anova(logHalictidae, test="F")
      rsq.partial(logHalictidae)
      rsq(logHalictidae)
    } #Dredge
    
    
    #logwb$Pollen_Plot<-NULL
    #loghb$Nectar_Plot<-NULL
    
    { 
    logHalictidae<-lm(logHalictidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                      + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=logwb)
    
    Fitall<- lm(logHalictidae~., data=logwb)
    Fitstart<-lm(logHalictidae~1, data=logwb)
    step(Fitstart, direction = "forward", scope=formula(Fitall))
    
    #Halictidae only
    wb<- lm (logHalictidae ~ Floral_Area, data = logwb)
    summary(wb)
    anova(wb, test="F")
    rsq.partial(wb)
    rsq(wb)
    AIC(wb)
    
    }#forward selection
   
    
    
    ##Model with Floral area is better model
    
    
  }#Halictidae model
  {
    #Megachilidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[22:22]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    
    
    {#Dredge first- Model Selection
      logMegachilidae<-lm(logMegachilidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                          + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(logMegachilidae)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      logMegachilidae<- lm( logMegachilidae ~ Week_Bloom, data=logwb)
      summary(logMegachilidae)
      anova(logMegachilidae, test="F")
      rsq.partial(logMegachilidae)
      rsq(logMegachilidae)
    } #Dredge
    
    {
    logMegachilidae<-lm(logMegachilidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                        + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=logwb)
    
    
    
    Fitall<- lm(logMegachilidae~., data=logwb)
    Fitstart<-lm(logMegachilidae~1, data=logwb)
    step(Fitstart, direction = "forward", scope=formula(Fitall))
    
    #Megachilidae only
    wb<- lm ( logMegachilidae ~ Week_Bloom + Nectar_Plot +Floral_Area + Pollen_Plot, 
              data = logwb)
    summary(wb)
    anova(wb, test="F")
    rsq.partial(wb)
    rsq(wb)
    AIC(wb)
    }#forward selection

  }#megachilidae model
  { 
    #Bumblebees only
    #this code will allow us to look at bumblebees specifically
    
    bee_genus<-dcast(specimen,site+year+week_num+block+plant_sp~genus, length)   
    head(bee_genus) 
    colnames(bee_genus)
    
    bee_genus<- ddply(bee_genus, c("plant_sp"), summarise,
                      Bombus= mean(Bombus))
    
    head(bee_genus)
    colnames(bee_genus)
    
    #merge the plant_t and specimen data files
    plants_bees <- merge(plants_bees, bee_genus)
    head(plants_bees)
    colnames(plants_bees)
    
    #Bumblebees only
    hist(plants_bees$Bombus)
    plants_bees$logbombus<-log(plants_bees$Bombus +1)
    hist(plants_bees$logbombus)
    
    #bumbles collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[28:28]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    {#Dredge first- Model Selection
      bombus<-lm(logbombus~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                 + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(bombus)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #honeybees
      bombus<- lm ( logbombus ~Flower_Height+Hue, data=logwb)
      summary(bombus)
      anova(bombus, test="F")
      rsq.partial(bombus)
      rsq(bombus)
    } #Dredge
    
    #logwb$Pollen_Plot<-NULL
    #loghb$Nectar_Plot<-NULL
    
    {
      logbombus<-lm(logbombus~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                    + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      #bumblebees only
      Fitall<- lm(logbombus~., data=logwb)
      Fitstart<-lm(logbombus~1, data=logwb)
      step(Fitstart, direction = "forward", scope=formula(Fitall))
      
      #BUmbles only
      wb<- lm (logbombus ~ Flower_Height + Hue + Floral_Area 
               , data = logwb)
      summary(wb)
      anova(wb, test="F")
      rsq.partial(wb)
      rsq(wb)
      AIC(wb)
      
    }#forward selection
    

  }#bumblebees only model
  { 
    colnames(plants_bees)
    #Non-Bombus Apidae
    ###need to create a column that includes Apidae minus Bombus for this analysis
    #create non-bombus apidae column
    plants_bees$nonbomusapidae<- plants_bees$Apidae-plants_bees$Bombus
    colnames(plants_bees)
    
    
    #non-bombus apidae only
    hist(plants_bees$nonbomusapidae)
    plants_bees$lognonbomusapidae<-log(plants_bees$nonbomusapidae +1)
    hist(plants_bees$lognonbomusapidae)
    
    #non-bombus apidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[32:32]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    
    {#Dredge first- Model Selection
      nonbombusapidae<-lm(lognonbomusapidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                          + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(nonbombusapidae)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      nonbombusapidae<- lm( lognonbomusapidae ~ Corolla_Width+ Floral_Area+Week_Bloom, data=logwb)
      summary(nonbombusapidae)
      anova(nonbombusapidae, test="F")
      rsq.partial(nonbombusapidae)
      rsq(nonbombusapidae)
    } #Dredge
    
    #logwb$Pollen_Plot<-NULL
    #loghb$Nectar_Plot<-NULL
    
    {
      
      lognonbombusapidae<-lm(lognonbomusapidae~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                             + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      #non-bombus apidae only
      Fitall<- lm(lognonbomusapidae~., data=logwb)
      Fitstart<-lm(lognonbomusapidae~1, data=logwb)
      step(Fitstart, direction = "forward", scope=formula(Fitall))
      
      #non-bombus apidae only
      wb<- lm (formula = lognonbomusapidae ~ Floral_Area + Corolla_Width + 
                 Week_Bloom, data = logwb)
      summary(wb)
      anova(wb, test="F")
      rsq.partial(wb)
      rsq(wb)
      AIC(wb)
      
    }#forward selection
   
    
  }#non-bombus apidae model
  {
    colnames(plants_bees)
    
    #species richness
    hist(plants_bees$SR)
    plants_bees$logsr<-log(plants_bees$SR +1)
    hist(plants_bees$logsr)
    
    #non-bombus apidae collected
    #Forward selection method
    #first I need to remove unecessary variables
    #first with floral area and not pollen
    colnames(plants_bees)
    
    plants_wb<-plants_bees[2:9]
    plants_wb2<-plants_bees[15:15]
    logwb<-cbind(plants_wb2,plants_wb)
    colnames((logwb))
    
    
    {#Dredge first- Model Selection
      logsr<-lm(SR~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
                + Corolla_Width + Chroma + Hue + Pollen_Plot, data=logwb)
      
      library(MuMIn)
      options(na.action = "na.fail")
      pos.model = dredge(logsr)
      best.mods = subset(pos.model, delta < 3)
      best.mods
      
      colnames(logwb)
      #wild bees no bombus
      logsr<- lm( SR ~ Floral_Area, data=logwb)
      summary(logsr)
      anova(logsr, test="F")
      rsq.partial(logsr)
      rsq(logsr)
    } #Dredge
    
    
    {
    
    logsr<-lm(SR~Week_Bloom + Number_Flowers + Floral_Area + Flower_Height
              + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=logwb)
    #species richness only
    Fitall<- lm(SR~., data=logwb)
    Fitstart<-lm(SR~1, data=logwb)
    step(Fitstart, direction = "forward", scope=formula(Fitall))
    
    #non-bombus apidae only
    wb<- lm (logsr ~ Floral_Area, data = logwb)
    summary(wb)
    anova(wb, test="F")
    rsq.partial(wb)
    rsq(wb)
    AIC(wb)
    
    }#forward selection

  }#species richness model
  
  {
    
    #honeybees
    
    honeybees<-lm(loghoneybees~ Floral_Area, data=plants_bees)
    summary(honeybees)
    plot(honeybees)
    
    bumblebees<-lm(logbombus~ Floral_Area, data=plants_bees)
    summary(bumblebees)
    plot(bumblebees)
    
    nonbombus<-lm(logwbnobombus~ Floral_Area, data=plants_bees)
    summary(nonbombus)
    plot(nonbombus)
    
    
    syrphids<-lm(logsyrpids~Floral_Area, data=plants_bees)
    summary(syrphids)
    
    
  }#simple linear regressions of broad pol groups
  
}#Multiple linear regressions to see trait if traits explain vistation #using taxa means/plantsp only

{
colnames(plants_bees)
bees_mat<-plants_bees[c(25,59,63,61,21,22,24,23,16)]
colnames(bees_mat)
plant_mat<-plants_bees[c(5,7,6,8,3,10)]
colnames(plant_mat)
bee_plant<-cbind(bees_mat,plant_mat)
colnames(bee_plant)

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
res2<-rcorr(as.matrix(bee_plant), type = c("pearson"))
flattenCorrMatrix(res2$r, res2$P)
#performance metrics
library("PerformanceAnalytics")
chart.Correlation(bee_plant, histogram=TRUE, pch=19)
}#matrix of bee and plant traits that were found to be significant

{#import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #SWMREC FIRST
  #this creates the specimen matrix with site, year, bloom period attached
  specimen<-dcast(specimen,family+genus+species ~plant_sp, length)   
  head(specimen) 
  colnames(specimen)
  #write.csv(swmrec1,"swmrec1.csv")
  
  
  ###NOT SURE YET
  tot_bees_by_plant<-rowSums(specimen[4:55])
  plants_bees_condensed<-specimen[tot_bees_by_plant>2,] #removing interactions zeros
  colnames(plants_bees_condensed)
  rowSums(plants_bees_condensed[4:57])
  
  write.csv(plants_bees_condensed,"swmrec2.csv")
  
  
  
  
  
  #creating bee and plant matrix data for PCA analysis
  all_plants<- specimen[4:57]
  colnames(swmrec_plants)
  head(swmrec_plants)
  
  all_bees<- specimen[1:3]
  colnames(swmrec_bees)
  
  #plot it
  ord<-metaMDS(all_plants, autotransform=FALSE)
  ord
  most_abund<-colSums(all_bees)>100
  plot(ord, disp='sites', type="n", cex=0.5)
  #title(main="All Sites", cex.main=1.5, adj=0)
  #display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
  points(ord, display="sites", select=which(all_bees$family=="Apidae"), pch=19, col="#a6cee3", cex=0.75)
  points(ord, display="sites", select=which(all_bees$family=="Andrenidae"), pch=19, col="#1f78b4", cex=0.75)
  points(ord, display="sites", select=which(all_bees$family=="Halictidae"), pch=19, col="#b2df8a", cex=0.75)
  points(ord, display="sites", select=which(all_bees$family=="Megachilidae"), pch=19, col="#33a02c", cex=0.75)
  points(ord, display="sites", select=which(all_bees$family=="Colletidae"), pch=19, col="#e3bba6", cex=0.75)
  levels(all_bees$family)=c("Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae")
  ordiellipse(ord, all_bees$family, draw="polygon", col="#a6cee3", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Apidae")
  ordiellipse(ord, all_bees$family, draw="polygon", col="#1f78b4", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Andrenidae")
  ordiellipse(ord, all_bees$family, draw="polygon", col="#b2df8a", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Halictidae")
  ordiellipse(ord, all_bees$family, draw="polygon", col="#33a02c", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Megachilidae")
  ordiellipse(ord, all_bees$family, draw="polygon", col="#e3bba6", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Colletidae")
  legend(.8,1, title=NULL, pch=c(19,19,19,19,19), ncol=1, text.width=0.8,
         col=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#e3bba6"), 
         cex=0.8, legend=c( "Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae" ))
  
  adonis(all_plants ~all_bees$family, method="bray",permutations=100)
}#ALL SITES NMDS
{#import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  #subset by site
  swmrec<-subset(specimen,site == "SWMREC")
  crc<-subset(specimen,site == "CRC")
  nwmhrc<-subset(specimen,site == "NWMHRC")
  
  
  
  #SWMREC FIRST
  #this creates the specimen matrix with site, year, bloom period attached
  swmrec1<-dcast(swmrec,family+genus+species ~plant_sp, length)   
  head(swmrec1) 
  colnames(swmrec1)
  #write.csv(swmrec1,"swmrec1.csv")
  
  
  ###NOT SURE YET
  tot_bees_by_plant<-rowSums(swmrec1[4:55])
  plants_bees_condensed<-swmrec1[tot_bees_by_plant>2,] #removing interactions zeros
  colnames(plants_bees_condensed)
  rowSums(plants_bees_condensed[4:55])
  
  write.csv(plants_bees_condensed,"swmrec2.csv")
  
  
  
  
  
  #creating bee and plant matrix data for PCA analysis
  swmrec_plants<- swmrec1[4:55]
  colnames(swmrec_plants)
  head(swmrec_plants)
  
  swmrec_bees<- swmrec1[1:3]
  colnames(swmrec_bees)
  
  #plot it
  ord<-metaMDS(swmrec_plants, autotransform=FALSE)
  ord
  most_abund<-colSums(swmrec_bees)>100
  plot(ord, disp='sites', type="n", cex=0.5)
  title(main="SWMREC", cex.main=1.5, adj=0)
  #display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
  points(ord, display="sites", select=which(swmrec_bees$family=="Apidae"), pch=19, col="#a6cee3", cex=0.75)
  points(ord, display="sites", select=which(swmrec_bees$family=="Andrenidae"), pch=19, col="#1f78b4", cex=0.75)
  points(ord, display="sites", select=which(swmrec_bees$family=="Halictidae"), pch=19, col="#b2df8a", cex=0.75)
  points(ord, display="sites", select=which(swmrec_bees$family=="Megachilidae"), pch=19, col="#33a02c", cex=0.75)
  points(ord, display="sites", select=which(swmrec_bees$family=="Colletidae"), pch=19, col="#e3bba6", cex=0.75)
  levels(swmrec_bees$family)=c("Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae")
  ordiellipse(ord, swmrec_bees$family, draw="polygon", col="#a6cee3", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Apidae")
  ordiellipse(ord, swmrec_bees$family, draw="polygon", col="#1f78b4", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Andrenidae")
  ordiellipse(ord, swmrec_bees$family, draw="polygon", col="#b2df8a", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Halictidae")
  ordiellipse(ord, swmrec_bees$family, draw="polygon", col="#33a02c", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Megachilidae")
  ordiellipse(ord, swmrec_bees$family, draw="polygon", col="#e3bba6", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Colletidae")
  legend(.7,1.2, title=NULL, pch=c(19,19,19,19,19), ncol=1, text.width=0.8,
         col=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#e3bba6"), 
         cex=0.8, legend=c( "Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae" ))
  
  adonis(swmrec_plants ~swmrec_bees$family, method="bray",permutations=100)
  
  ##NOW CRC
  #CRC Second
  #this creates the specimen matrix with site, year, bloom period attached
  crc1<-dcast(crc,family+genus+species ~plant_sp, length)   
  head(crc1) 
  colnames(crc1)
  #write.csv(crc1,"crc1.csv")
  
  
  ###NOT SURE YET
  tot_bees_by_plant<-rowSums(swmrec1[4:55])
  plants_bees_condensed<-swmrec1[tot_bees_by_plant>2,] #removing interactions zeros
  colnames(plants_bees_condensed)
  rowSums(plants_bees_condensed[4:55])
  
  write.csv(plants_bees_condensed,"swmrec2.csv")
  
  
  
  
  
  #creating bee and plant matrix data for PCA analysis
  crc_plants<- crc1[4:52]
  colnames(crc_plants)
  head(crc_plants)
  
  crc_bees<- crc1[1:3]
  colnames(crc_bees)
  
  #plot it
  ord<-metaMDS(crc_plants, autotransform=FALSE)
  ord
  most_abund<-colSums(crc_bees)>100
  plot(ord, disp='sites', type="n", cex=0.5)
  title(main="CRC", cex.main=1.5, adj=0)
  #display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
  points(ord, display="sites", select=which(crc_bees$family=="Apidae"), pch=19, col="#a6cee3", cex=0.75)
  points(ord, display="sites", select=which(crc_bees$family=="Andrenidae"), pch=19, col="#1f78b4", cex=0.75)
  points(ord, display="sites", select=which(crc_bees$family=="Halictidae"), pch=19, col="#b2df8a", cex=0.75)
  points(ord, display="sites", select=which(crc_bees$family=="Megachilidae"), pch=19, col="#33a02c", cex=0.75)
  points(ord, display="sites", select=which(crc_bees$family=="Colletidae"), pch=19, col="#e3bba6", cex=0.75)
  levels(crc_bees$family)=c("Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae")
  ordiellipse(ord, crc_bees$family, draw="polygon", col="#a6cee3", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Apidae")
  ordiellipse(ord, crc_bees$family, draw="polygon", col="#1f78b4", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Andrenidae")
  ordiellipse(ord, crc_bees$family, draw="polygon", col="#b2df8a", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Halictidae")
  ordiellipse(ord, crc_bees$family, draw="polygon", col="#33a02c", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Megachilidae")
  ordiellipse(ord, crc_bees$family, draw="polygon", col="#e3bba6", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Colletidae")
  legend(1,1, title=NULL, pch=c(19,19,19,19,19), ncol=1, text.width=0.8,
         col=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#e3bba6"), 
         cex=0.8, legend=c( "Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae" ))
  
  adonis(crc_plants ~ crc_bees$family, method="bray",permutations=1000)
  ##NOW NWMHRC
  #NWMHRC Third
  #this creates the specimen matrix with site, year, bloom period attached
  nwmhrc1<-dcast(nwmhrc,family+genus+species ~plant_sp, length)   
  head(nwmhrc1) 
  colnames(nwmhrc1)
  #write.csv(nwmhrc1,"nwmhrc1.csv")
  
  
  ###NOT SURE YET
  #tot_bees_by_plant<-rowSums(swmrec1[4:55])
  #plants_bees_condensed<-swmrec1[tot_bees_by_plant>2,] #removing interactions zeros
  #colnames(plants_bees_condensed)
  #rowSums(plants_bees_condensed[4:55])
  
  #write.csv(plants_bees_condensed,"swmrec2.csv")
  
  
  
  
  
  #creating bee and plant matrix data for PCA analysis
  nwmhrc_plants<- nwmhrc1[4:51]
  colnames(nwmhrc_plants)
  head(nwmhrc_plants)
  
  nwmhrc_bees<- nwmhrc1[1:3]
  colnames(nwmhrc_bees)
  
  #plot it
  ord<-metaMDS(nwmhrc_plants, autotransform=FALSE)
  ord
  most_abund<-colSums(nwmhrc_bees)>100
  plot(ord, disp='sites', type="n", cex=0.5)
  title(main="NWMHRC", cex.main=1.5, adj=0)
  #display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
  points(ord, display="sites", select=which(nwmhrc_bees$family=="Apidae"), pch=19, col="#a6cee3", cex=0.75)
  points(ord, display="sites", select=which(nwmhrc_bees$family=="Andrenidae"), pch=19, col="#1f78b4", cex=0.75)
  points(ord, display="sites", select=which(nwmhrc_bees$family=="Halictidae"), pch=19, col="#b2df8a", cex=0.75)
  points(ord, display="sites", select=which(nwmhrc_bees$family=="Megachilidae"), pch=19, col="#33a02c", cex=0.75)
  points(ord, display="sites", select=which(nwmhrc_bees$family=="Colletidae"), pch=19, col="#e3bba6", cex=0.75)
  levels(nwmhrc_bees$family)=c("Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae")
  ordiellipse(ord, nwmhrc_bees$family, draw="polygon", col="#a6cee3", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Apidae")
  ordiellipse(ord, nwmhrc_bees$family, draw="polygon", col="#1f78b4", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Andrenidae")
  ordiellipse(ord, nwmhrc_bees$family, draw="polygon", col="#b2df8a", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Halictidae")
  ordiellipse(ord, nwmhrc_bees$family, draw="polygon", col="#33a02c", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Megachilidae")
  ordiellipse(ord, nwmhrc_bees$family, draw="polygon", col="#e3bba6", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="Colletidae")
  legend(1,1.7, title=NULL, pch=c(19,19,19,19,19), ncol=1, text.width=0.8,
         col=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#e3bba6"), 
         cex=0.8, legend=c( "Apidae", "Andrenidae", "Halictidae", "Megachilidae","Colletidae" ))
  
  adonis(nwmhrc_plants ~ nwmhrc_bees$family, method="bray",permutations=1000)
  
  
}#SITE SPECIFIC NMDS PLOTS
{
  plants<-read.csv("data/plant_traits.csv")
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants$nec_plot<-plants$nectar*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  head(plants)
  
  plants_t<-ddply(plants, c("site","plant_sp"), summarise,
                  Week_Bloom = mean(week_num),
                  Number_Flowers = mean(tot_flw),
                  #avg_flw_area= mean(avg_flor_area),
                  Floral_Area= mean(tot_flor_area),
                  Flower_Height= mean(avg_tall_flw),
                  #Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot),
                  Nectar_Plot=mean(nec_plot))
  plants_t #this the plants trait file
  head(plants_t)
  
  #import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  
  colnames(specimen) #make sure it looks correct
  
  #this creates the specimen matrix with site, year, bloom period attached
  bee_genus<-dcast(specimen,site+plant_sp ~genus, length)   
  head(bee_genus) 
  colnames(bee_genus)
  
  #merge the plant_t and specimen data files
  plants_bees <- merge(plants_t, bee_genus)
  head(plants_bees)
  
  #remove NAs because NMDS can't hang
  plants_bees<- na.omit(plants_bees)
  head(plants_bees)
  colnames(plants_bees)
  
  #look at observations with fewer bees- we need to remove ones that didn't capture much
  
  plants_bees<-plants_bees[which(rowSums(plants_bees[12:42])>5),]
  
  #look at genus sums- will remove genera with very low numbers
  colSums(plants_bees[12:42])
  
  
  #remove genera with low numbers of bees 
  #plants_bees$Anthidiellum<-NULL #0 total
  colnames(plants_bees)
  
  
  
  #tot_bees<-colSums(plants_bees[12:146])
  #plants_bees_condensed<-plants_bees[tot_bees<2,] #removing interactions zeros
  #colSums(plants_bees_condensed[12:146])
  
  
  
  #creating bee and plant matrix data for PCA analysis
  bee_matrix<- plants_bees[12:42]
  colnames(bee_matrix)
  head(bee_matrix)
  #1 per genus
  bee_matrix$Anthidiellum<-NULL
  bee_matrix$Dufourea<-NULL
  bee_matrix$Calliopsis<-NULL
  bee_matrix$Dianthidium<-NULL
  bee_matrix$Peponapis<-NULL
  #two per genus
  bee_matrix$Augochlora<-NULL
  bee_matrix$Coelioxys<-NULL
  bee_matrix$Protandrena<-NULL
  bee_matrix$Triepeolus<-NULL
  # three to 25
  bee_matrix$Anthophora<-NULL
  bee_matrix$Colletes<-NULL
  bee_matrix$Eucera<-NULL
  bee_matrix$Nomada<-NULL
  bee_matrix$`Psuedo panurgus`<-NULL
  bee_matrix$Heriades<-NULL
  bee_matrix$Hoplitis<-NULL
  bee_matrix$Osmia<-NULL
  bee_matrix$Perdita<-NULL
  bee_matrix$`Sphecodes `<-NULL
  bee_matrix$Xylocopa<-NULL
  bee_matrix$Anthidium<-NULL
  
  plant_data<- plants_bees[1:11]
  colnames(plant_data)
  
  #remove unneccesary information
  plant_data$plant_sp<-NULL
  plant_data$site<-NULL
  plant_data$year<-NULL
  plant_data$block<-NULL
  plant_data$Nectar_Plot<-NULL

  #plant_data$Pollen_Plot<-NULL
  
  ##PCA #I like this the most
  ord<-rda(bee_matrix, distance="bray", autotransform=T, scale=T)
  ordiplot(ord, display="species", type= "t")
  summary(ord)
  
  ef <- envfit(ord, plant_data, permu = 1000, repel = TRUE)
  ef
 labels(ef)
  plot(ef)
  scores(ef, "sites")
  summary(ef)
  com<-adonis(bee_matrix ~Week_Bloom+Number_Flowers+Floral_Area+Flower_Height+Corolla_Width+Chroma+Hue+Pollen_Plot, data=plant_data,  method="bray",permutations=1000)
  com
  
  summary(com)
  
  aov(com, type=F)
  
  
}##PCA all sites combined
{
  #import specimen dataset
  specimen <- read.csv('dat/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  
  #create bee matrix
  bee_matrix<-dcast(specimen,site+plant_sp ~genus, length)   
  head(bee_matrix) 
  colnames(bee_matrix)
  
  
  #import plant data
  plants<-read.csv("data/plant_traits.csv")
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants$nec_plot<-plants$nectar*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  head(plants)
  
  plants_t<-ddply(plants, c("site","plant_sp"), summarise,
                  Week_Bloom = mean(week_num),
                  Number_Flowers = mean(tot_flw),
                  #avg_flw_area= mean(avg_flor_area),
                  Floral_Area= mean(tot_flor_area),
                  Flower_Height= mean(avg_tall_flw),
                  #Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot),
                  Nectar_Plot=mean(nec_plot))
  plants_t #this the plants trait file
  head(plants_t)
  
  #combine the plants and the bee files
  plants_bees<-merge(plants_t,bee_matrix)
  head(plants_bees)
  colnames(plants_bees)
  
  
  
  
  #subset by site
  swmrec<-subset(plants_bees,site == "SWMREC")
  crc<-subset(plants_bees,site == "CRC")
  nwmhrc<-subset(plants_bees,site == "NWMHRC")
  
  
  {
    head(swmrec)
    colnames(swmrec)
    
    
    
    #creating bee and plant matrix data for PCA/NMDS analysis
    swmrec_traits<- swmrec[3:10]
    colnames(swmrec_traits)
    head(swmrec_traits)
    
    swmrec_bees<- swmrec[12:42]
    colnames(swmrec_bees)
    library(vegan)
    
    tot_bees_by_plant<-colSums(swmrec_bees[1:31])
    swmrec_condensed<-swmrec_bees[,tot_bees_by_plant>2] #removing interactions zeros
    colnames(swmrec_condensed)

    
    ##PCA #I like this the most
    ord<-rda(swmrec_condensed, distance="bray", autotransform=T, scale=T)
    ord
    plot(ord, display="sites", type= "p")
    summary(ord)

    ef <- envfit(ord, swmrec_traits, permu = 1000)
    ef
    plot(ef)
    scores(ef, "vectors")
    print(scores)
    
    
    table<-summary(ord)
    table
    write.csv(table,"table.csv")
    
    #community analysis
    swmrec_com<-adonis(swmrec_condensed ~ Week_Bloom + Number_Flowers+ Floral_Area + Flower_Height
                       + Corolla_Width + Chroma + Hue + Pollen_Plot, data=swmrec_traits,  method="bray",permutations=1000)
    swmrec_com
    
  }#SWMREC PCA/ADONIS
  {
    head(crc)
    colnames(crc)
    
    
    
    #creating bee and plant matrix data for PCA/NMDS analysis
    crc_traits<- crc[3:11]
    colnames(crc_traits)
    head(crc_traits)
    
    crc_bees<- crc[12:146]
    colnames(crc_bees)
    library(vegan)
    
    tot_bees_by_plant<-colSums(crc_bees[1:135])
    crc_condensed<-crc_bees[,tot_bees_by_plant>2] #removing interactions zeros
    colnames(crc_condensed)
    
    
    ##PCA #I like this the most
    ord<-rda(crc_condensed, distance="bray", autotransform=T, scale=T)
    ord
    plot(ord, display="sites", type= "p")
    summary(ord)
    
    ef <- envfit(ord, crc_traits, permu = 1000)
    ef
    plot(ef)
    scores(ef, "vectors")
    print(scores)
    
    
    table<-summary(ord)
    table
    write.csv(table,"table.csv")
    
    #community analysis
    crc_com<-adonis(crc_condensed ~ Week_Bloom + Number_Flowers+ Floral_Area + Flower_Height
                       + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=crc_traits,  method="bray",permutations=1000)
    crc_com
  }#CRC PCA/ADONIS
  {
    head(nwmhrc)
    colnames(nwmhrc)
    
    
    
    #creating bee and plant matrix data for PCA/NMDS analysis
    nwmhrc_traits<- nwmhrc[3:11]
    colnames(nwmhrc_traits)
    head(nwmhrc_traits)
    
    nwmhrc_bees<- nwmhrc[12:146]
    colnames(nwmhrc_bees)
    library(vegan)
    
    tot_bees_by_plant<-colSums(nwmhrc_bees[1:135])
    nwmhrc_condensed<-nwmcrc_bees[,tot_bees_by_plant>2] #removing interactions zeros
    colnames(nwmhrc_condensed)
    
    
    ##PCA #I like this the most
    ord<-rda(nwmhrc_condensed, distance="bray", autotransform=T, scale=T)
    ord
    plot(ord, display="sites", type= "p")
    summary(ord)
    
    ef <- envfit(ord, nwmhrc_traits, permu = 1000)
    ef
    plot(ef)
    scores(ef, "vectors")
    print(scores)
    
    
    table<-summary(ord)
    table
    write.csv(table,"table.csv")
    
    #community analysis
    nwmhrc_com<-adonis(nwmhrc_condensed ~ Week_Bloom + Number_Flowers+ Floral_Area + Flower_Height
                    + Corolla_Width + Chroma + Hue + Pollen_Plot + Nectar_Plot, data=nwmhrc_traits,  method="bray",permutations=1000)
    nwmhrc_com
  }#NWMHRC PCA/ADONIS
}#SITE SPECIFIC PCA / ADONIS
{
  ##import specimen dataset
  specimen <- read.csv('data/all_bees3.csv') 
  head(specimen)
  colnames(specimen) #make sure it looks correct
  
  
  #create bee matrix
  bee_matrix<-dcast(specimen,site+plant_sp ~family, length)   
  head(bee_matrix) 
  colnames(bee_matrix)
  
  
  #import plant data
  plants<-read.csv("data/plant_traits.csv")
  head(plants)
  plants$pol_plot<-plants$pol.unit*plants$tot_flw
  plants$nec_plot<-plants$nectar*plants$tot_flw
  plants<-na.omit(plants)
  plants$per_tar_cov<-as.numeric(plants$per_tar_cov)
  head(plants)
  
  plants_t<-ddply(plants, c("site","plant_sp"), summarise,
                  Week_Bloom = mean(week_num),
                  Number_Flowers = mean(tot_flw),
                  #avg_flw_area= mean(avg_flor_area),
                  Floral_Area= mean(tot_flor_area),
                  Flower_Height= mean(avg_tall_flw),
                  #Percent_Tar_Cov= mean(per_tar_cov),
                  Corolla_Width=mean(avg_fw),
                  Chroma=mean(chroma),
                  Hue=mean(hue),
                  #avg_pol_flw= mean(pol.unit),
                  Pollen_Plot=mean(pol_plot),
                  Nectar_Plot=mean(nec_plot))
  plants_t #this the plants trait file
  head(plants_t)
  colnames(plants_t)
  plants<-plants_t[c(3:11)]
  #combine the plants and the bee files
  plants_bees<-merge(plants_t,bee_matrix)
  head(plants_bees)
  colnames(plants_bees)
  
  #subset to get only the floral area data
  Floral_area<-plants_bees[c(2,5,10,12:16)]
  write.csv(Floral_area, "floralarea.csv")
  
  library("reshape")
  #reshape for figure
  reshape<-melt(Floral_area, id=c("plant_sp","Floral_Area","Pollen_Plot"))
  reshape$specimen<-reshape$variable
  reshape$variable<-NULL
  write.csv(reshape,"reshaped.csv")
  
  #makefigure
  library(ggplot2)
  xyplot <- ggplot(reshape, aes(Pollen_Plot , value))+
    geom_point(aes())+
    xlab("Average Floral Area")+
    ylab("Total number collected")+
    theme_bw()+ #removes grey background
    
    #stat_smooth(method=lm, size=.1)
    
    geom_smooth(aes(colour= factor(specimen)),method="lm", se=FALSE) #adds treadline to each factor
  
  xyplot
  
}#GGPLOT of floral area and visitation by bee families

