setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")
setwd("~/Dropbox/SESYNC/SESYNC_RACs/For R Package")
#all data
pdata<-read.csv('pplots_example.csv')

#no replicate
pdata2<-subset(pdata, replicate==1)
pdata2<-pdata2[,-3]

#no treatment
pdata3<-subset(pdata, treatment=="N1P0")
pdata3<-pdata3[-2]

#no time
pdata4<-subset(pdata, time==2002)
pdata4 <- pdata4[-1]

#adding fake block
#rep_trt <- unique(subset(df, select = c(replicate.var, treatment.var)))
#write.csv(rep_trt,"forfakeblocks.csv")
fb<-read.csv("forfakeblocks2.csv")
pdata5<-merge(fb, pdata, by=c("treatment","replicate"))

##fake blocks no time
pdata6<-subset(pdata5, time==2002)
pdata6<-pdata6[-4]

replicate.var <- 'replicate'
treatment.var <- 'treatment'
species.var <- 'species'
time.var <- 'time'
abundance.var <- 'abundance'
block.var <- 'block'



#community_diversity
#works with time and rep.
test1.shan<-community_diversity(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance")
test1.simp<-community_diversity(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance", diversity = "Simpson")
#works with no rep
test2.simp<-community_diversity(pdata2, time.var="time", abundance.var = "abundance", diversity = "Simpson")
test2.shan<-community_diversity(pdata2, time.var="time", abundance.var = "abundance")
#works with no time
test3.shan<-community_diversity(pdata, replicate.var = "replicate", abundance.var = "abundance")
test3.simp<-community_diversity(pdata, replicate.var = "replicate", abundance.var = "abundance", diversity = "Simpson")

#community_structure
test1.eq<-community_structure(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance")
test1.esimp<-community_structure(pdata, time.var="time", replicate.var = "replicate", abundance.var = "abundance", evenness = "SimpEven")
#works with no rep
test2.esimp<-community_structure(pdata2, time.var="time", abundance.var = "abundance", evenness = "SimpEven")
test2.eq<-community_structure(pdata2, time.var="time", abundance.var = "abundance")
#works with no time
test3.esimp<-community_structure(pdata4, abundance.var = "abundance", replicate.var = "replicate", evenness = "SimpEven")
test3.eq<-community_structure(pdata4, abundance.var = "abundance", replicate.var = "replicate")


#RAC_changes
#works with time and rep.
test1<-RAC_change(pdata, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance")
#works with no rep
test2<-RAC_change(pdata2, time.var="time",species.var = "species", abundance.var = "abundance")

#multivariate change
#with treatment
test2<-multivariate_change(pdata, time.var="time", replicate.var = "replicate", treatment.var = "treatment", species.var = "species", abundance.var = "abundance")
#without treatment
test1<-multivariate_change(pdata3, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance")

#curve change
test1<-curve_change(pdata, time.var="time", abundance.var = "abundance", species.var = "species", replicate.var="replicate")
test2<-curve_change(pdata2, time.var="time", abundance.var = "abundance", species.var = "species")

#abundance change
test1<-abundance_change(pdata, time.var="time", abundance.var = "abundance", species.var = "species", replicate.var="replicate")
test2<-abundance_change(pdata2, time.var="time", abundance.var = "abundance", species.var = "species")

####
#RAC difference
####

#with blocks no time
test1<-RAC_difference(df = pdata6, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', block.var = "block", replicate.var = "replicate")
#with blocks and time
test2<-RAC_difference(df = pdata5, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', block.var = "block", replicate.var = "replicate", time.var = "time")

#pooling by treatment no time
test3<-RAC_difference(df = pdata4, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', pool="YES", replicate.var = "replicate")
#pooling by treatment with time
test4<-RAC_difference(df = pdata, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', pool="YES", replicate.var = "replicate", time.var = "time")

#all replicates no trt no time
test5<-RAC_difference(df=pdata4, species.var = "species", abundance.var = "abundance", replicate.var = "replicate")
#all replicates time
test6<-RAC_difference(df=pdata, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", time.var="time")

#all replicates with trt but no time
test7<-RAC_difference(df=pdata4, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", treatment.var = 'treatment')
#all replicates time
test8<-RAC_difference(df=pdata, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", time.var="time", treatment.var= 'treatment')

####
#Abundance difference
####

#with blocks no time
test1<-abundance_difference(df = pdata6, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', block.var = "block", replicate.var = "replicate")
#with blocks and time
test2<-abundance_difference(df = pdata5, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', block.var = "block", replicate.var = "replicate", time.var = "time")

#pooling by treatment no time
test3<-abundance_difference(df = pdata4, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', pool="YES", replicate.var = "replicate")
#pooling by treatment with time
test4<-abundance_difference(df = pdata, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', pool="YES", replicate.var = "replicate", time.var = "time")

#all replicates no trt no time
test5<-abundance_difference(df=pdata4, species.var = "species", abundance.var = "abundance", replicate.var = "replicate")
#all replicates time
test6<-abundance_difference(df=pdata, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", time.var="time")

#all replicates with trt but no time
test7<-abundance_difference(df=pdata4, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", treatment.var = 'treatment')
#all replicates time
test8<-abundance_difference(df=pdata, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", time.var="time", treatment.var= 'treatment')


####
#curve difference
####

#with blocks no time
test1<-curve_difference(df = pdata6, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', block.var = "block", replicate.var = "replicate")
#with blocks and time
test2<-curve_difference(df = pdata5, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', block.var = "block", replicate.var = "replicate", time.var = "time")

#pooling by treatment no time
test3<-curve_difference(df = pdata4, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', pool="YES", replicate.var = "replicate")
#pooling by treatment with time
test4<-curve_difference(df = pdata, species.var = "species", abundance.var = "abundance", treatment.var = 'treatment', pool="YES", replicate.var = "replicate", time.var = "time")

#all replicates no trt no time
test5<-curve_difference(df=pdata4, species.var = "species", abundance.var = "abundance", replicate.var = "replicate")
#all replicates time
test6<-curve_difference(df=pdata, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", time.var="time")

#all replicates with trt but no time
test7<-curve_difference(df=pdata4, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", treatment.var = 'treatment')
#all replicates time
test8<-curve_difference(df=pdata, species.var = "species", abundance.var = "abundance", replicate.var = "replicate", time.var="time", treatment.var= 'treatment')

#multivariate difference
#without time
test2<-multivariate_difference(pdata4, replicate.var = "replicate", treatment.var = "treatment", species.var = "species", abundance.var = "abundance")
#with time
test1<-multivariate_difference(pdata, time.var="time", replicate.var = "replicate", species.var = "species", abundance.var = "abundance", treatment.var = "treatment")
