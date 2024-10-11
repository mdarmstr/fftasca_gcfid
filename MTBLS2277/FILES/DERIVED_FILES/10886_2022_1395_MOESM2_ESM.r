### Analysis of Tribolium castaneum chemical profiles - submitted to Journal of Chemical Ecology, manuscript title:
#     "Immune stimulation alters chemical profiles of adult Tribolium castaneum"
# Lai Ka LO, Reshma R, Lisa J. TEWES, Barbara MILUTINOVIC, Caroline M?LLER, Joachim KURTZ

#############################################################################
#############################################################################
#############################################################################
########## PERMANOVA on gland secretion profiles of Tribolium castaneum 
########## of different immune status (naive, wounded, primed)
#############################################################################

#loading required package
library(vegan) 		

# writing function for pairwise comparison according to 
# Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4. https://github.com/pmartinezarbizu/pairwiseAdonis

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
library(vegan)

co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()


for(elem in 1:ncol(co)){
if(sim.function == 'daisy'){
library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
} else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model =c(F.Model,ad$aov.tab[1,4]);
R2 = c(R2,ad$aov.tab[1,5]);
p.value = c(p.value,ad$aov.tab[1,6])
}
p.adjusted = p.adjust(p.value,method=p.adjust.m)
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
print("Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1")
return(pairw.res)

} 


#############################################################################
# previously subsetted datasets and separate tables with treatment group 
# levels were loaded, subset is indicated by the name of the datafiles
#############################################################################

#### COMPLETE dataset 24 h after treatment
data <-read.csv("secretions_data24h_all_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("secretions_treatments24h_all_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Sex, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA (testing additive effects after interaction was not significant)
adonis(data_ws~Sex*Treatment, permutations = 10000, method="kulczynski")
adonis(data_ws~Sex+Treatment, permutations = 10000, method="kulczynski")

#############################################################################
#### FEMALES ONLY dataset 24 h after treatment
data <-read.csv("secretions_data24h_females_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("secretions_treatments24h_females_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA
adonis(data_ws~Treatment, permutations = 10000, method="kulczynski")

# pairwise comparison after significant treatment effect
pairwise.adonis(data_ws,Treatment, sim.method="kulczynski", p.adjust.m="BH")


#############################################################################
#### MALES ONLY dataset 24 h after treatment
data <-read.csv("secretions_data24h_males_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("secretions_treatments24h_males_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA (no significant treatment effect)
adonis(data_ws~Treatment, permutations = 10000, method="kulczynski")


#############################################################################
#### COMPLETE dataset 72 h after treatment
data <-read.csv("secretions_data72h_all_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("secretions_treatments72h_all_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Sex, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA (testing additive effects after interaction was not significant)
adonis(data_ws~Sex*Treatment, permutations = 10000, method="kulczynski")
adonis(data_ws~Sex+Treatment, permutations = 10000, method="kulczynski")

#############################################################################
#### FEMALES ONLY dataset 72 h after treatment
data <-read.csv("secretions_data72h_females_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("secretions_treatments72h_females_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA (no significant treatment effect)
adonis(data_ws~Treatment, permutations = 10000, method="kulczynski")

#############################################################################
#### MALES ONLY dataset 72 h after treatment
data <-read.csv("secretions_data72h_males_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("secretions_treatments72h_males_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA (no significant treatment effect)
adonis(data_ws~Treatment, permutations = 10000, method="kulczynski")


#############################################################################
#############################################################################
#############################################################################
########## PERMANOVA on cuticular hydrocarbon profiles of Tribonium castaneum 
########## of different immune status (naive, wounded, primed)
#############################################################################

#### COMPLETE dataset
data <-read.csv("CHC_data_all_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("CHC_treatments_all_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Sex, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA (testing additive effects after interaction was not significant)
adonis(data_ws~Sex*Treatment, permutations = 10000, method="kulczynski")
adonis(data_ws~Sex+Treatment, permutations = 10000, method="kulczynski")

#############################################################################
#### Comparison of sex differnces in NAIVE beetles only
data <-read.csv("CHC_data_naive_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("CHC_treatments_naive_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Sex, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA
adonis(data_ws~Sex, permutations = 10000, method="kulczynski")

#############################################################################
#### FEMALES ONLY 
data <-read.csv("CHC_data_females_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("CHC_treatments_females_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA
adonis(data_ws~Treatment, permutations = 10000, method="kulczynski")

# pairwise comparison after significant treatment effect
pairwise.adonis(data_ws,Treatment, sim.method="kulczynski", p.adjust.m="BH")

#############################################################################
#### MALES ONLY dataset 24 h after treatment
data <-read.csv("CHC_data_males_relative.csv", header=T, row.names=1, check.names=F)
attach(treatments <-read.csv("CHC_treatments_males_relative.csv", header=T))

# wisconsin-square-root transformation
data_ws <- wisconsin(sqrt(data))

# test of homogeneity of the multivariate variance spread among the treatment groups
mod <- betadisper(vegdist(data_ws,method="kulczynski"), Treatment, type="median")
permutest(mod, pairwise = TRUE, permutations = 10000)

# PERMANOVA
adonis(data_ws~treat, permutations = 10000, method="kulczynski")

# pairwise comparison after significant treatment effect
pairwise.adonis(data_ws,Treatment, sim.method="kulczynski", p.adjust.m="BH")

#############################################################################

#############################################################################
#############################################################################
#############################################################################
### Analysis of T.castaneum secretion profiles
#############################################################################
#############################################################################

# Since the initial NMDS plots and analysis (Fig. 1a, b) showed large sex differences and a
# non-significant model for for the data acquired 72 h after the treatment, we continued 
# our analysis of compounds on the data acquired 24 h after the treatment and separately
# for each sex

# To identify candidate compounds that contribute to treatment differences we use random forest
# analysis on relative compound abundance.


rm(list=ls())# removes all info in workspace 
getwd()#confirms that everything is reset
options(scipen=999)# stops R from writing long numbers in sci. format. 
setwd("C:/Users/Bmilutin/Dropbox/LaikaPhD/Rfiles/publication_script")

#libraries
library('randomForest')
library("multcomp")
library("car")

#Load the 24 h dataset (relative compound amount)
data=read.csv("secretions_relative.csv", header=T)
data24=droplevels(subset(data,Timepoint == "24h"))
View(data24)
summary(data24)

data24[sapply(data24, is.character)] <- lapply(data24[sapply(data24, is.character)], 
                                           as.factor)

### Random forest:
#24 h Females
data24F=droplevels(subset(data24,Sex == "female"))
treatment=data24F$Treatment
data24F=data24F[,5:30]#remove character variables

set.seed(245)
rf=randomForest(treatment ~., data= data24F, importance=T, ntrees=500, keep.forest=FALSE) 
rf 
#Important compound suggestion plot:
varImpPlot(rf, sort=TRUE, n.var=min(26, nrow(rf$importance)), type=1, scale=F, main="female relative compounds")
  #feature RI1015 seems contributing to treatment differences in females 

#24 h Males
data24M=droplevels(subset(data24,Sex == "male"))
treatment=data24M$Treatment
data24M=data24M[,5:30]#remove character variables

set.seed(245)
rf=randomForest(treatment ~., data= data24M, importance=T, ntrees=500, keep.forest=FALSE) 
rf 
#Important compound suggestion plot:
varImpPlot(rf, sort=TRUE, n.var=min(26, nrow(rf$importance)), type=1, scale=F, main="male relative compounds")
  #feature RI1473 seems contributing to treatment differences in males


### Individual models on two rf-suggested (RI1015, RI1472) and three candidate secretion compounds (MBQ, EBQ, C15)
#Reshape data for narrow format
library("reshape")
data=melt(data24, id = 1:4)
names(data)[5]<-paste("compound")#rename column 5
names(data)[6]<-paste("abundance")#rename column 6
str(data)
head(data, 10)

#relevel so that wounded and primed are compared to naive 
data$treatment=relevel(data$Treatment,ref="naive")

#MBQ
dataMBQ=droplevels(subset(data,compound == "Methyl.1.4.benzoquinone"))
dataMBQ$Treatment=relevel(dataMBQ$Treatment,ref="naive")
boxplot(abundance~Treatment*Sex,data=dataMBQ)
#MODEL
full=lm(abundance~Treatment*Sex,data=dataMBQ)
null=lm(abundance~1,data=dataMBQ)
Anova(full)
anova(null, full)
summary(full)

#model diagnostics
diagnostics.plot<-function(mod.res){
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 1, 0.5))
  hist(residuals(mod.res), probability=T, xlab="", ylab="", main="")
  mtext(text="histogram of residuals", side=3, line=0)
  x=seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  qqnorm(residuals(mod.res), main="", pch=19)
  qqline(residuals(mod.res))
  mtext(text="qq-plot of residuals", side=3, line=0)
  plot(fitted(mod.res), residuals(mod.res), pch=19)
  abline(h=0, lty=2)
  mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))#potential problem
max(cooks.distance(full))
max(abs(dfbeta(full)))#potential problem

#Transform to see if this gets better
library("bestNormalize")
bn = bestNormalize(dataMBQ$abundance)# center_scale(x) best transformation
hist(bn$x.t)
dataMBQ=cbind(dataMBQ, bn$x.t)
ncol(dataMBQ) #[1] 7
names(dataMBQ)[7]<-paste("transf")
head(dataMBQ, 10)
#new model with box cox transformed data
full=lm(transf~Treatment*Sex,data=dataMBQ)
null=lm(transf~1,data=dataMBQ)
#diagnostics
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))# better
max(cooks.distance(full))
max(abs(dfbeta(full)))# better
Anova(full)
anova(null, full)
# Transformed data model gives similar values as the model with the raw data, so we can stay with the first.
#The model for relative MBQ compound is not significant


#EBQ
dataEBQ=droplevels(subset(data,compound == "Ethyl.1.4.benzoquinone"))
dataEBQ$Treatment=relevel(dataEBQ$Treatment,ref="naive")
boxplot(abundance~Treatment*Sex,data=dataEBQ)
#MODEL
full=lm(abundance~Treatment*Sex,data=dataEBQ)
null=lm(abundance~1,data=dataEBQ)
Anova(full)
anova(null, full)
summary(full)
#diagnostics
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))
max(cooks.distance(full))
max(abs(dfbeta(full)))
#The model for relative EBQ compound is not significant

#C-15
dataC15=droplevels(subset(data,compound == "X1.C15.ene"))
dataC15$Treatment=relevel(dataC15$Treatment,ref="naive")
boxplot(abundance~Treatment*Sex,data=dataC15)
#MODEL
full=lm(abundance~Treatment*Sex,data=dataC15)
null=lm(abundance~1,data=dataC15)
Anova(full)
anova(null, full)
summary(full)
#diagnostics
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))#potential problem
max(cooks.distance(full))
max(abs(dfbeta(full)))
#transform to see if coefficients look better
bn = bestNormalize(dataC15$abundance)# Standardized asinh(x) best transformation
hist(bn$x.t)
dataC15=cbind(dataC15, bn$x.t)
ncol(dataC15) #[1] 7
names(dataC15)[7]<-paste("transf")
head(dataC15, 10)

#new model with transformed data
full=lm(transf~Treatment*Sex,data=dataC15)
null=lm(transf~1,data=dataC15)
anova(null, full)
summary(full)
#diagnostics
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))#better
max(cooks.distance(full))
max(abs(dfbeta(full)))
#The model for relative C-15 compound is not significant


#RI1015
dataRI1015=droplevels(subset(data,compound == "RI1015"))
dataRI1015$Treatment=relevel(dataRI1015$Treatment,ref="naive")
boxplot(abundance~Treatment*Sex,data=dataRI1015)
#MODEL
full=lm(abundance~Treatment*Sex,data=dataRI1015)
null=lm(abundance~1,data=dataRI1015)
Anova(full)
anova(null, full)
summary(full)
#diagnostics
diagnostics.plot(full)#not good
cor.test(fitted(full),abs(residuals(full)))#not good

#transform to see if diagnostics look better
bn = bestNormalize(dataRI1015$abundance)# Standardized sqrt(x + a) best transformation
hist(bn$x.t)#looks good
dataRI1015=cbind(dataRI1015, bn$x.t)
ncol(dataRI1015) #[1] 7
names(dataRI1015)[7]<-paste("transf")
head(dataRI1015, 10)

#new model with transformed data
full=lm(transf~treatment*sex,data=dataRI1015)
null=lm(transf~1,data=dataRI1015)
anova(null, full)
summary(full)
#diagnostics
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#not good
car::ncvTest(full)#not good
#Since model diagnostics is not good, we perform separate Kruskal-Wallis and Pairwise Wilcoxon Rank Sum tests
# for each sex

dataRI1015f=droplevels(subset(dataRI1015,Sex == "female"))
dataRI1015m=droplevels(subset(dataRI1015,Sex == "male"))
#Females
kruskal.test(dataRI1015f$abundance~dataRI1015f$Treatment)
#Kruskal test for females is significant, we continue with Benjamini-Hochberg posthoc tests
pairwise.wilcox.test(dataRI1015f$abundance,dataRI1015f$Treatment,p.adjust.method="BH")#both wounded and primed females differ from males, but not between each other
#Males
kruskal.test(dataRI1015m$abundance~dataRI1015m$Treatment)#model for males is not significant


#RI1473
dataRI1473=droplevels(subset(data,compound == "RI1473"))
dataRI1473$Treatment=relevel(dataRI1473$Treatment,ref="naive")
boxplot(abundance~Treatment*Sex,data=dataRI1473)
#Since naive and primed males have zero or very low of this compound, we perform separate Kruskal-Wallis and 
# Pairwise Wilcoxon Rank Sum tests for each sex immediately

dataRI1473f=droplevels(subset(dataRI1473,Sex == "female"))
dataRI1473m=droplevels(subset(dataRI1473,Sex == "male"))
#Females
kruskal.test(dataRI1473f$abundance~dataRI1473f$Treatment)#model for females is not significant
#Males
kruskal.test(dataRI1473m$abundance~dataRI1473m$Treatment)
pairwise.wilcox.test(dataRI1473m$abundance,dataRI1473m$Treatment,p.adjust.method="BH")#Relative abundance of compound
# RI1473 differs between naive and wounded males and only trending between primed and wounded

#overall model p-values were corrected using Benjamini-Hochberg procedure:
p.adjust(c(0.2081, 0.3906,0.659,0.007305,0.4228,0.1798,0.009694), method="BH")
#[1] 0.3641750 0.4932667 0.6590000 0.0339290 0.4932667 0.3641750 0.0339290

# Absolute abundance of three candidate secretion compounds (MBQ, EBQ, C15)
rm(list=ls())# removes all info in workspace 
getwd()#confirms that everything is reset
options(scipen=999)# stops R from writing long numbers in sci. format. 
setwd("C:/Users/Bmilutin/Dropbox/LaikaPhD/Rfiles/publication_script")

data=read.csv("secretions_absolute.csv", header=T)
str(data)
data24=droplevels(subset(data,Timepoint == "24h")) # for 24 h time point
data72=droplevels(subset(data,Timepoint == "72h")) # for 72 h time point

#re-shape the data for long format
library("reshape")
data=melt(data24, id = 1:4) # for 24 h time point
data=melt(data72, id = 1:4) # for 72 h time point
names(data)[5]<-paste("compound")#rename column 5
names(data)[6]<-paste("abundance")#rename column 6
str(data)
head(data, 10)

data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], 
                                           as.factor)
levels(data$compound)
#MBQ - absolute
dataMBQ=droplevels(subset(data,compound == "Methyl.1.4.benzoquinone"))
dataMBQ$Treatment=relevel(dataMBQ$Treatment,ref="naive")
#Model
full=lm(abundance~Treatment*Sex,data=dataMBQ)
null=lm(abundance~1,data=dataMBQ)
Anova(full)
anova(null, full)
summary(full)
#assumptions
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#not good

#transform
library(bestNormalize)
bn = bestNormalize(dataMBQ$abundance)
hist(bn$x.t)
dataMBQ=cbind(dataMBQ, bn$x.t)
ncol(dataMBQ)
#[1] 7
names(dataMBQ)[7]<-paste("transf")
head(dataMBQ)
#model with transformed data
full=lm(transf~Treatment*Sex,data=dataMBQ)
null=lm(transf~1,data=dataMBQ)
#assumptions
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#better
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))
max(cooks.distance(full))
max(abs(dfbeta(full)))
#use the model on transformed data
anova(null, full)#overall model significant
Anova(full)#interaction not significant
summary(full)

#Separate models on females and males to obtain posthocs
dataMBQf=droplevels(subset(dataMBQ,Sex == "female"))
dataMBQm=droplevels(subset(dataMBQ,Sex == "male"))
#F
full=lm(transf~Treatment,data=dataMBQf)
null=lm(transf~1,data=dataMBQf)
anova(full,null)#model for females significant
summary(full)
boxplot(abundance~Treatment,data=dataMBQf)
summary(glht(full, linfct=mcp(Treatment="Tukey")), test= adjusted("BH"))
#wounded is different from naives, and trending versus primed (corrected p-value 0.052)

#M
full=lm(transf~Treatment,data=dataMBQm)
null=lm(transf~1,data=dataMBQm)
anova(full,null)#model for males is not significant
summary(full)

#EBQ - absolute
dataEBQ=droplevels(subset(data,compound == "Ethyl.1.4.benzoquinone"))
dataEBQ$Treatment=relevel(dataEBQ$Treatment,ref="naive")
#Model
full=lm(abundance~Treatment*Sex,data=dataEBQ)
null=lm(abundance~1,data=dataEBQ)
Anova(full)
anova(null, full)
summary(full)
#assumptions
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#not good

#transform
library(bestNormalize)
bn = bestNormalize(dataEBQ$abundance)
hist(bn$x.t)
dataEBQ=cbind(dataEBQ, bn$x.t)
ncol(dataEBQ)
#[1] 7
names(dataEBQ)[7]<-paste("transf")
head(dataEBQ)
#model with transformed data
full=lm(transf~Treatment*Sex,data=dataEBQ)
null=lm(transf~1,data=dataEBQ)
#assumptions
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#better
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))
max(cooks.distance(full))
max(abs(dfbeta(full)))
#use the model on transformed data
anova(null, full)#overall model significant
Anova(full)#interaction not significant
summary(full)

#Separate models on females and males to obtain posthocs
dataEBQf=droplevels(subset(dataEBQ,Sex == "female"))
dataEBQm=droplevels(subset(dataEBQ,Sex == "male"))
#F
full=lm(transf~Treatment,data=dataEBQf)
null=lm(transf~1,data=dataEBQf)
anova(full,null)#model for females significant
summary(full)
boxplot(abundance~Treatment,data=dataEBQf)
summary(glht(full, linfct=mcp(Treatment="Tukey")), test= adjusted("BH"))
#wounded is different from naives (similar to MBQ), and trending versus primed (corrected p-value 0.058)

#M
full=lm(transf~Treatment,data=dataEBQm)
null=lm(transf~1,data=dataEBQm)
anova(full,null)#model for males is not significant
summary(full)


#C15 - absolute
dataC15=droplevels(subset(data,compound == "X1.C15.ene"))
dataC15$Treatment=relevel(dataC15$Treatment,ref="naive")
#Model
full=lm(abundance~Treatment*Sex,data=dataC15)
null=lm(abundance~1,data=dataC15)
Anova(full)
anova(null, full)
summary(full)
#assumptions
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#again not good

#transform
library(bestNormalize)
bn = bestNormalize(dataC15$abundance)
hist(bn$x.t)
dataC15=cbind(dataC15, bn$x.t)
ncol(dataC15)
#[1] 7
names(dataC15)[7]<-paste("transf")
head(dataC15)
#model with transformed data
full=lm(transf~Treatment*Sex,data=dataC15)
null=lm(transf~1,data=dataC15)
#assumptions
diagnostics.plot(full)
cor.test(fitted(full),abs(residuals(full)))#better
car::ncvTest(full)
max(abs(dffits(full)))
coefficients(full)+t(apply(dfbeta(full),2,range))
max(cooks.distance(full))
max(abs(dfbeta(full)))
#use the model on transformed data
anova(null, full)#overall model significant
Anova(full)#interaction not significant
summary(full)

#Separate models on females and males to obtain posthocs
dataC15f=droplevels(subset(dataC15,Sex == "female"))
dataC15m=droplevels(subset(dataC15,Sex == "male"))
#F
full=lm(transf~Treatment,data=dataC15f)
null=lm(transf~1,data=dataC15f)
anova(full,null)#model for females significant
summary(full)
boxplot(abundance~Treatment,data=dataC15f)
summary(glht(full, linfct=mcp(Treatment="Tukey")), test= adjusted("BH"))
#wounded is different from naives (similar to MBQ), and trending versus primed (corrected p-value 0.058)

#M
full=lm(transf~Treatment,data=dataC15m)
null=lm(transf~1,data=dataC15m)
anova(full,null)#model for males is not significant
summary(full)

#Overall models were corrected Benjamini-Hochberg
p.adjust(c(0.0117, 0.02185,0.001716), method="BH")
