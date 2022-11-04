library(LogisticDx)
#install.packages("matrixStats")
library(matrixStats)
#install.packages('BAGofT')
#library(BAGofT)
#install.packages('blorr') 
library(blorr)   #BIC
#install.packages('DescTools')#
library(DescTools)#all types of R2 page 427   #lecessies test page 253
#install.packages('modEvA')
library(modEvA)  #threshMeasures(model = fm, thresh=0.5) p39
#install.packages('rms')
library(rms) #lecessies test page 160

#install.packages('MKmisc')
library(MKmisc) #for cHCH different link

## Install package BiocManager
#install.packages("BiocManager")
## Use BiocManager to install limma
#BiocManager::install("limma")
#install.packages('limma')
#library(limma)


#################running settings A1 t0 A5################
############## Ommission of quadratic form ###############
######### model contains one linear predictor#############

##############################################################
#                          Matrices of power and null
##########################                          ####################
############################   run it Only at the first time#########
gofs_T=matrix(nr=1,nc=9)       ######################################
gofs_F=matrix(nr=1,nc=9)       ######################################
#                              #run it Only at the first time#######
auc_mean_T=matrix(nr=1,nc=3)   #######################################
auc_mean_F=matrix(nr=1,nc=3)   #######################################
auc_std_T=matrix(nr=1,nc=3)    #######################################
auc_std_F=matrix(nr=1,nc=3)    #######################################
#                              #run it Only at the first time#######
r2_mean_T=matrix(nr=1,nc=10)   ########################################
r2_mean_F=matrix(nr=1,nc=10)   ########################################
r2_std_T=matrix(nr=1,nc=10)    ########################################
r2_std_F=matrix(nr=1,nc=10)    ########################################
#                              #run it Only at the first time#######
bic_mean_T=matrix(nr=1,nc=1)   ########################################
bic_mean_F=matrix(nr=1,nc=1)   ########################################
bic_std_T=matrix(nr=1,nc=1)    ########################################
bic_std_F=matrix(nr=1,nc=1)    ########################################
#                              #run it Only at the first time#######
HL_T=matrix(nr=1,nc=3)         ########################################
HL_F=matrix(nr=1,nc=3)         ########################################
##                             ######################################## 
HL_T_2=matrix(nr=1,nc=3)       ######################################## 
HL_F_2=matrix(nr=1,nc=3)       ######################################## 
#                              #run it Only at the first time#######                              #run it Only at the first time#######
cvh2_T_res=matrix(nr=1,nc=1)   ########################################
cvh2_F_res=matrix(nr=1,nc=1)   ########################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
#########################################################################

###############################
############################### Matriceses of simulation



sim=1000 #number of simulations
pval_T=matrix(0,nr=sim,nc=9) #pvalue of true model    #LogisticDx package
pval_F=matrix(0,nr=sim,nc=9)#pvalue of false model

auc_T=matrix(0,nr=sim,nc=3)                           #LogisticDx package
auc_F=matrix(0,nr=sim,nc=3)


r2_T=matrix(0,nr=sim,nc=10)                          # desk tool package
r2_F=matrix(0,nr=sim,nc=10)

bic_T=matrix(0,nr=sim,nc=1)                          # blorr package
bic_F=matrix(0,nr=sim,nc=1)


HL_gp_T=matrix(0,nr=sim,nc=3)                       # desk tool package
HL_gp_F=matrix(0,nr=sim,nc=3)

HL_gp_2_T=matrix(0,nr=sim,nc=3)                    #MKmisc package for chch with glm
HL_gp_2_F=matrix(0,nr=sim,nc=3)

cvh2_T= matrix(0,nr=sim,nc=1)                        # rms package
cvh2_F= matrix(0,nr=sim,nc=1)                  #lecessie van holwegian test 

#for test collection
#tests_T=matrix(0,nr=sim,nc=13)
#tests_F=matrix(0,nr=sim,nc=13)
#NOT_tests_T
#Not_tests_f


##################################settings of quadratic term omission

n=c(100,200,300,400,500,600,700,800,900,1000) 
#b0=c(-1.1,-2,-2.3,-2.7,-3.2)
#b1=c(1.3,1,.9,.7,.6)
#b2=c(0,0.2,.3,.4,.5)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

################################## simulation


#setseed_check=matrix(nr=1,nc=1)


#                                     ######   
#x1=runif(n[2],-3,3)                  ######Change these settings of n and betas
#x2=x1^2                               ###### where n has 10 settings
#XB1 = -1.963 + 0.982*x1 + 0.218*x2     ###### and betas have 5 settings
#                                   ######                                             


#set.seed(4)


for (j in 1:10){
  set.seed(4)
  x1=runif(n[j],-6,6)               ######Change these settings of n and betas
  x2=rnorm(n[j],0,2.25)                               ###### where n has 10 settings
  x3=rchisq(n[j],4)
  XB =  0.267*x1 + 0.267*x2+0.217*x3
  pr = 1/(1+exp(-XB))
   
for (k in 1:sim){
  #set.seed(4)
  
  
  y= rbinom(length(x1),1,pr)
  td=data.frame(y,x1,x2,x3) #td=True Dataset
  tm=glm(y ~ x1+x2+x3,data=td, family=binomial(link='logit')) #tm= True model
  fd=data.frame(y,x1,x2) #fd= False data
  fm=glm(y~x1+x2,data=fd,family=binomial(link = 'logit')) #fm = False model
  
  
  gt=gof(tm)
  gf=gof(fm)                           #logistic dx
  pval_T[k,]=as.matrix(gt$gof[,5])
  pval_F[k,]=as.matrix(gf$gof[,5])
  
  
  auc_T[k,]=t(as.matrix(gt$auc))
  auc_F[k,]=t(as.matrix(gf$auc))
  
  r2_T[k,]=t(as.matrix(PseudoR2(tm,c('McFadden','McFaddenAdj','CoxSnell',
                                     'Nagelkerke','AldrichNelson','VeallZimmermann',
                                     'McKelveyZavoina','Efron','Tjur','AIC'))))
  r2_F[k,]=t(as.matrix(PseudoR2(fm,c('McFadden','McFaddenAdj','CoxSnell',
                                     'Nagelkerke','AldrichNelson','VeallZimmermann',
                                     'McKelveyZavoina','Efron','Tjur','AIC'))))
    
  bic_T[k]=t(as.matrix(blr_model_fit_stats(tm)[11]))
  bic_F[k]=t(as.matrix(blr_model_fit_stats(fm)[11]))
  
  
  HLc_F=as.numeric(HosmerLemeshowTest(fit = fitted(fm), obs = y, X = cbind(x1))$C[3])
  HLh_F=as.numeric(HosmerLemeshowTest(fit = fitted(fm), obs = y, X = cbind(x1))$H[3])
  cvh1_F=as.numeric(HosmerLemeshowTest(fit = fitted(fm), obs = y, X = cbind(x1))$gof[2])
  HL_gp_F[k,]=t(as.matrix(cbind( HLc_F,  HLh_F,cvh1_F)))
  #desctool page 253) 
   
   
  HLc_T=as.numeric(HosmerLemeshowTest(fit = fitted(tm), obs = y, X = cbind(x1,x2))$C[3])
  HLh_T=as.numeric(HosmerLemeshowTest(fit = fitted(tm), obs = y, X = cbind(x1,x2))$H[3])
  cvh1_T=as.numeric(HosmerLemeshowTest(fit = fitted(tm), obs = y, X = cbind(x1,x2))$gof[2])
  HL_gp_T[k,]=t(as.matrix(cbind(HLc_T,HLh_T,cvh1_T)))
 
  
  HLchch_F=as.numeric(HLgof.test(fit = fitted(fm), obs = y, X = model.matrix(y ~ x1))$gof[2])
  HLc2_F=as.numeric(HLgof.test(fit = fitted(fm), obs = y, X = model.matrix(y ~ x1))$C[3])
  HLh2_F=as.numeric(HLgof.test(fit = fitted(fm), obs = y, X = model.matrix(y ~ x1))$H[3])
  HL_gp_2_F[k,]=t(as.matrix(cbind(HLc2_F,HLh2_F,HLchch_F)))
  #MKmisc package for chch with glm
  
  HLchch_T=as.numeric(HLgof.test(fit = fitted(tm), obs = y, X = model.matrix(y ~ x1+x2))$gof[2])
  HLc2_T=as.numeric(HLgof.test(fit = fitted(tm), obs = y, X = model.matrix(y ~ x1+x2))$C[3])
  HLh2_T=as.numeric(HLgof.test(fit = fitted(tm), obs = y, X = model.matrix(y ~ x1+x2))$H[3])
  HL_gp_2_T[k,]=t(as.matrix(cbind(HLc2_T,HLh2_T,HLchch_T)))
 
  
  
  tm2 <- lrm(y ~ x1+x2+x3 , x=TRUE, y=TRUE)
  cvh2_T[k]=as.matrix(resid(tm2, "gof")[5])
  fm2 <- lrm(y ~ x1+x2 , x=TRUE, y=TRUE)
  cvh2_F[k]=as.matrix(resid(fm2, "gof")[5])
   
  
  
  #pred_T <- predict(tm, type = "response")
  #pred_F <- predict(fm, type = "response")
  
  #perfMeasures(pred_T, truth = td$y, namePos = 1)[,2]   ##MKmisk package
  #perfScores(pred_T, truth = td$y, namePos = 1) [,2]
  
  #perfMeasures(pred_F, truth = fd$y, namePos = 1)[,2]
  #perfScores(pred_F, truth = fd$y, namePos = 1) [,2]
  
  
  
  #st=sig(tm)
  #sf=sig(fm)  ## logisticdx package sif function
  #st$Wald
  #st$LR
  #st$score
  
}


################################################

#                      getting the power and null from tests
#                                      and
#                      mean and std from non tests

f=matrix(0,nr=sim,nc=9) #f= logical matrix of p value >0.05 in false model
t=matrix(0,nr=sim,nc=9) #t= logical matrix of p value >0.05 in true model
for (u in 1:9){ (f[,u]=pval_F[,u]>0.05)
                (t[,u]=pval_T[,u]>0.05)}                      
power=as.matrix(1-colMeans(f))
null_hyp=as.matrix(colMeans(t))
gofs_F=rbind(gofs_F,t(power))
gofs_T=rbind(gofs_T,t(null_hyp))

auc_mean_F=rbind(auc_mean_F,colMeans(auc_F))
auc_std_F=rbind(auc_std_F,colSds(auc_F))
auc_mean_T=rbind(auc_mean_T,colMeans(auc_T))
auc_std_T=rbind(auc_std_T,colSds(auc_T))

r2_mean_F=rbind(r2_mean_F,colMeans(r2_F))
r2_std_F=rbind(r2_std_F,colSds(r2_F))
r2_mean_T=rbind(r2_mean_T,colMeans(r2_T))
r2_std_T=rbind(r2_std_T,colSds(r2_T))

bic_mean_F=rbind(bic_mean_F,colMeans (as.matrix((as.numeric(bic_F)))))
bic_std_F=rbind(bic_std_F,colSds (as.matrix((as.numeric(bic_F)))))
bic_mean_T=rbind(bic_mean_T,colMeans (as.matrix((as.numeric(bic_T)))))
bic_std_T=rbind(bic_std_T,colSds (as.matrix((as.numeric(bic_T)))))

fff=matrix(0,nr=sim,nc=3) #f= logical matrix of p value >0.05 in false model
ttt=matrix(0,nr=sim,nc=3) #t= logical matrix of p value >0.05 in true model
for (u in 1:3){ (fff[,u]=HL_gp_F[,u]>0.05)
                (ttt[,u]=HL_gp_T[,u]>0.05)}                      
HL_power=as.matrix(1-colMeans(fff))
HL_null_hyp=as.matrix(colMeans(ttt))
HL_F=rbind(HL_F,t(HL_power))
HL_T=rbind(HL_T,t(HL_null_hyp))



ffff=matrix(0,nr=sim,nc=3) #f= logical matrix of p value >0.05 in false model
tttt=matrix(0,nr=sim,nc=3) #t= logical matrix of p value >0.05 in true model
for (u in 1:3){ (ffff[,u]=HL_gp_2_F[,u]>0.05)
  (tttt[,u]=HL_gp_2_T[,u]>0.05)}                      
HL_power_2=as.matrix(1-colMeans(ffff))
HL_null_hyp_2=as.matrix(colMeans(tttt))
HL_F_2=rbind(HL_F_2,t(HL_power_2))
HL_T_2=rbind(HL_T_2,t(HL_null_hyp_2))


ff=matrix(0,nr=sim,nc=1) #f= logical matrix of p value >0.05 in false model
tt=matrix(0,nr=sim,nc=1) #t= logical matrix of p value >0.05 in true model
for (u in 1:sim){ (ff[u]=cvh2_F[u]>0.05)
                  (tt[u]=cvh2_T[u]>0.05) }
cvh2_power=as.matrix(1-colMeans(ff))
cvh2_null_hyp=as.matrix(colMeans(tt))
cvh2_F_res=rbind(cvh2_F_res,t(cvh2_power))
cvh2_T_res=rbind(cvh2_T_res,t(cvh2_null_hyp))
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
######################  transforming matrices to dataframes 
############################       and 
###################           writing csv files



#Tests=as.data.frame(cbind(rbind(gofs_F,gofs_T),rbind(HL_F,HL_T),
                         #rbind (cvh2_F_res,cvh2_T_res)))

#Not_Tests=as.data.frame(cbind(rbind(auc_mean_F,auc_std_F,auc_mean_T,auc_std_T),
                              #rbind(r2_mean_F,r2_std_F,r2_mean_T,r2_std_T),
                              #rbind(bic_mean_F,bic_std_F,bic_mean_T,bic_std_T)))
#write.csv(Tests,'Tests.csv')
#write.csv(Not_Tests,'Not_Tests.csv')


Tests_Mariam=as.data.frame(cbind(gofs_F,gofs_T,HL_F,HL_T,HL_F_2,HL_T_2,cvh2_F_res,cvh2_T_res))
Not_Tests_Mariam=as.data.frame(cbind(auc_mean_F,auc_std_F,auc_mean_T,auc_std_T,r2_mean_F,r2_std_F,
                              r2_mean_T,r2_std_T,bic_mean_F,bic_std_F,bic_mean_T,bic_std_T))
write.csv(Tests_Mariam,'Tests_Mariam_quad.csv')
write.csv(Not_Tests_Mariam,'Not_Tests_Mariam_quad.csv')

###############################################################################

#####################

