library(semPlot)
library(lavaan) 
library(MASS)
library(leaps)
library(readr)

MEMEX <- read_csv("/Users/simonwdavis/Dropbox/Duke/Memex/Writing/Network_Neuroscience_paper/SEM_modeling/MEMEX_471_FA_62Ss.csv")

# ROIs: 1-2: ATL; 3-4: CST; 5-6: Cingulum; 7-8: CingHipp; 9: Genu; 10: Splenium; 11-12: IFOF; 13-14: ILF; 15-16: SLF; 17-18: UF; 19-20: SLFtemp; 21: fornix

MEMEX$DTI1 = (scale(MEMEX$DTI1) + scale(MEMEX$DTI2))/2
MEMEX$DTI3 = (scale(MEMEX$DTI3) + scale(MEMEX$DTI4))/2
MEMEX$DTI5 = (scale(MEMEX$DTI5) + scale(MEMEX$DTI6))/2
MEMEX$DTI7 = (scale(MEMEX$DTI7) + scale(MEMEX$DTI8))/2
MEMEX$DTI9 = scale(MEMEX$DTI9)
MEMEX$DTI10 = scale(MEMEX$DTI10)
MEMEX$DTI11 = (scale(MEMEX$DTI11) + scale(MEMEX$DTI12))/2
MEMEX$DTI13 = (scale(MEMEX$DTI13) + scale(MEMEX$DTI14))/2
MEMEX$DTI15 = (scale(MEMEX$DTI15) + scale(MEMEX$DTI16))/2
MEMEX$DTI17 = (scale(MEMEX$DTI17) + scale(MEMEX$DTI18))/2
MEMEX$DTI19 = (scale(MEMEX$DTI19) + scale(MEMEX$DTI20))/2
MEMEX$DTI21 = scale(MEMEX$DTI21)

MEMEX$ITM1 = (scale(MEMEX$ITM1) + scale(MEMEX$ITM2))/2
MEMEX$ITM3 = (scale(MEMEX$ITM3) + scale(MEMEX$ITM4))/2
MEMEX$ITM5 = (scale(MEMEX$ITM5) + scale(MEMEX$ITM6))/2
MEMEX$ITM7 = (scale(MEMEX$ITM7) + scale(MEMEX$ITM8))/2
MEMEX$ITM9 = scale(MEMEX$ITM9)
MEMEX$ITM10 = scale(MEMEX$ITM10)
MEMEX$ITM11 = (scale(MEMEX$ITM11) + scale(MEMEX$ITM12))/2
MEMEX$ITM13 = (scale(MEMEX$ITM13) + scale(MEMEX$ITM14))/2
MEMEX$ITM15 = (scale(MEMEX$ITM15) + scale(MEMEX$ITM16))/2
MEMEX$ITM17 = (scale(MEMEX$ITM17) + scale(MEMEX$ITM18))/2
MEMEX$ITM19 = (scale(MEMEX$ITM19) + scale(MEMEX$ITM20))/2
MEMEX$ITM21 = scale(MEMEX$ITM21)

MEMEX$SRC1 = (scale(MEMEX$SRC1) + scale(MEMEX$SRC2))/2
MEMEX$SRC3 = (scale(MEMEX$SRC3) + scale(MEMEX$SRC4))/2
MEMEX$SRC5 = (scale(MEMEX$SRC5) + scale(MEMEX$SRC6))/2
MEMEX$SRC7 = (scale(MEMEX$SRC7) + scale(MEMEX$SRC8))/2
MEMEX$SRC9 = scale(MEMEX$SRC9)
MEMEX$SRC10 = scale(MEMEX$SRC10)
MEMEX$SRC11 = (scale(MEMEX$SRC11) + scale(MEMEX$SRC12))/2
MEMEX$SRC13 = (scale(MEMEX$SRC13) + scale(MEMEX$SRC14))/2
MEMEX$SRC15 = (scale(MEMEX$SRC15) + scale(MEMEX$SRC16))/2
MEMEX$SRC17 = (scale(MEMEX$SRC17) + scale(MEMEX$SRC18))/2
MEMEX$SRC19 = (scale(MEMEX$SRC19) + scale(MEMEX$SRC20))/2
MEMEX$SRC21 = scale(MEMEX$SRC21)

MEMEX$Src = scale(MEMEX$Src)
MEMEX$Item = scale(MEMEX$Item)
MEMEX$agecat <- ifelse(MEMEX$Age > 40, c("older"), c("younger")) 


############################################################
########## Full models #####################################
############################################################

SEM_SRC.model <- '
WM =~ DTI17 + DTI7 +  DTI21 + DTI9
SRC =~ SRC17 + SRC7 + SRC21 + SRC9
WM ~~ SRC
Src ~ WM + SRC
'
fit4 <- sem(SEM_SRC.model, data=subset(MEMEX, Age > 40),  estimator="ML") 
# fit44 <- cv_regsem(fit5, pars_pen=NULL, n.lambda=15, type="lasso", jump=0.3) , modindices = TRUE
summary(fit4, fit.measures = TRUE, standardized=TRUE, rsquare=TRUE)
fitmeasures(fit4)
semPaths(fit4, "std", rotation = 2)
anova(fit4, fit4uc)

###################
SEM_ITEM.model <- '
WM =~ DTI5 + DTI11 + DTI13 + DTI17
ITM =~ ITM5 + ITM11 + ITM13 + ITM17
WM ~~ ITM
Item ~ WM + ITM
'
fit5 <- sem(SEM_ITEM.model, data=subset(MEMEX, Age > 40), estimator="ML")
# fit55 <- cv_regsem(fit5, max.try=3, n.lambda=15, lambda=0.1, type="lasso", gradFun="ram")
summary(fit5, fit.measures = TRUE, standardized=TRUE,rsquare=TRUE)
fitmeasures(fit5)
semPaths(fit5, "std", rotation = 2)


