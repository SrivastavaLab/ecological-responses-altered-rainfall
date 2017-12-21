#import fulldata from 01_accessing_data.R

library(lme4)
library(lmerTest)
library(MASS)
library(vcd)
library(fitdistrplus)
library(car)
library(knitr)
library(visreg)
library(vegan)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(Hmisc)
library(quantreg)
library(MuMIn)
library(readr)
library(tibble)

agreement <- Agreement::agreement
#need to use package agreement function agreement to do Lin's concordance and partition acuracy and precision
#Lawrence Lin, A. S Hedayat, Bikas Sinha, Min Yang. Journal of the American Statistical Associa-
#tion. March 1, 2002, 97(457): 257-270



#=====don't skip!making new variables========
#cardoso has an "engulfing" monopelopia (Diptera 175) - switch to piercing - this has now been done
#macae had a grazer, this has now been corrected


#======Dsquared function====
Dsquared <- function(model = NULL, 
                     obs = NULL, 
                     pred = NULL, 
                     family = NULL, # needed only when 'model' not provided
                     adjust = FALSE, 
                     npar = NULL) { # needed only when 'model' not provided
  # version 1.4 (31 Aug 2015)
  
  model.provided <- ifelse(is.null(model), FALSE, TRUE)
  
  if (model.provided) {
    if (!("glm" %in% class(model))) stop ("'model' must be of class 'glm'.")
    if (!is.null(pred)) message("Argument 'pred' ignored in favour of 'model'.")
    if (!is.null(obs)) message("Argument 'obs' ignored in favour of 'model'.")
    obs <- model$y
    pred <- model$fitted.values
    
  } else { # if model not provided
    if (is.null(obs) | is.null(pred)) stop ("You must provide either 'obs' and 'pred', or a 'model' object of class 'glm'.")
    if (length(obs) != length(pred)) stop ("'obs' and 'pred' must be of the same length (and in the same order).")
    if (is.null(family)) stop ("With 'obs' and 'pred' arguments (rather than a model object), you must also specify one of two model family options: 'binomial' or 'poisson' (in quotes).")
    else if (!is.character(family)) stop ("Argument 'family' must be provided as character (i.e. in quotes: 'binomial' or 'poisson').")
    else if (length(family) != 1 | !(family %in% c("binomial", "poisson"))) stop ("'family' must be either 'binomial' or 'poisson' (in quotes).")
    
    if (family == "binomial") {
      if (any(!(obs %in% c(0, 1)) | pred < 0 | pred > 1)) stop ("'binomial' family implies that 'obs' data should be binary (with values 0 or 1) and 'pred' data should be bounded between 0 and 1.")
      link <- log(pred / (1 - pred))  # logit
    }  # end if binomial
    
    else if (family == "poisson") {
      if (any(obs %%1 != 0)) stop ("'poisson' family implies that 'obs' data should consist of whole numbers.")
      link <- log(pred)
    }  # end if poisson
    
    model <- glm(obs ~ link, family = family)
  }  # end if model not provided
  
  D2 <- (model$null.deviance - model$deviance) / model$null.deviance
  
  if (adjust) {
    if (model.provided) {
      n <- length(model$y)
      #p <- length(model$coefficients)
      p <- attributes(logLik(model))$df
    } else {
      if (is.null(npar)) stop ("Adjusted D-squared from 'obs' and 'pred' values (rather than a model object) requires specifying the number of parameters in the underlying model ('npar').")
      n <- length(na.omit(obs))
      p <- npar
    }  # end if model.provided else
    
    D2 <- 1 - ((n - 1) / (n - p)) * (1 - D2)
  }  # end if adjust
  
  return (D2)
}

#==========multisite model comparisons

aic.lmx<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol)+site, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+log(mu.scalar), family=family, data = dataset)
  m2<-glm(y~log(maxvol)+site+log(k.scalar), family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), family=family, data = dataset)
  m4<-glm(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), family=family, data = dataset)
  m5<-glm(y~log(maxvol)+site*log(mu.scalar), family=family, data = dataset)
  m6<-glm(y~log(maxvol)+site*log(k.scalar), family=family, data = dataset)
  m7<-glm(y~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m8<-glm(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m9<-glm(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m10<-glm(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m11<-glm(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m12<-glm(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m13<-glm(y~log(maxvol)+site*log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m14<-glm(y~log(maxvol)+site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m15<-glm(y~log(maxvol)+site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m16<-glm(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  aic.mod<-c(m0$aic, m1$aic, m2$aic, m3$aic, m4$aic, m5$aic, m6$aic, m7$aic, m8$aic, m9$aic, m10$aic,m11$aic, m12$aic, m13$aic, m14$aic, m15$aic, m16$aic)
   print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16))
}

aic.lmxnb<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site, data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+log(mu.scalar), data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+log(k.scalar), data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), data = dataset)
  m5<-glm.nb(y~log(maxvol)+site*log(mu.scalar),  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site*log(k.scalar), data = dataset)
  m7<-glm.nb(y~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m8<-glm.nb(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m9<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar),  data = dataset)
  m10<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m11<-glm.nb(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m12<-glm.nb(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m13<-glm.nb(y~log(maxvol)+site*log(mu.scalar)*log(k.scalar),  data = dataset)
  m14<-glm.nb(y~log(maxvol)+site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m15<-glm.nb(y~log(maxvol)+site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m16<-glm.nb(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16))
}

aic.lmxnb.init<-function(y,init.theta, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site, init.theta=init.theta, data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+log(mu.scalar),init.theta=init.theta, data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+log(k.scalar),init.theta=init.theta, data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), init.theta=init.theta,data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), init.theta=init.theta,data = dataset)
  m5<-glm.nb(y~log(maxvol)+site*log(mu.scalar), init.theta=init.theta, data = dataset)
  m6<-glm.nb(y~log(maxvol)+site*log(k.scalar), init.theta=init.theta,data = dataset)
  m7<-glm.nb(y~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), init.theta=init.theta,data = dataset)
  m8<-glm.nb(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), init.theta=init.theta,data = dataset)
  m9<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar), init.theta=init.theta, data = dataset)
  m10<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), init.theta=init.theta,data = dataset)
  m11<-glm.nb(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), init.theta=init.theta,data = dataset)
  m12<-glm.nb(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), init.theta=init.theta,data = dataset)
  m13<-glm.nb(y~log(maxvol)+site*log(mu.scalar)*log(k.scalar),init.theta=init.theta,  data = dataset)
  m14<-glm.nb(y~log(maxvol)+site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), init.theta=init.theta,data = dataset)
  m15<-glm.nb(y~log(maxvol)+site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), init.theta=init.theta,data = dataset)
  m16<-glm.nb(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), init.theta=init.theta,data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16))
}

aic.lmxnb.add<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site, data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+log(mu.scalar), data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+log(k.scalar), data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), data = dataset)
  m9<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar),  data = dataset)
  m10<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m11<-glm.nb(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m12<-glm.nb(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m9, m10, m11, m12))
}

aic.lmx.add<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol)+site, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+log(mu.scalar), family=family,data = dataset)
  m2<-glm(y~log(maxvol)+site+log(k.scalar),family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), family=family,data = dataset)
  m4<-glm(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), family=family,data = dataset)
  m9<-glm(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar),family=family,data = dataset)
  m10<-glm(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family,data = dataset)
  m11<-glm(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family,data = dataset)
  m12<-glm(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family,data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m9, m10, m11, m12))
}

aic.lmx.nb.best<-function(a, b, scalar)
{
  nbset<-aic.lmxnb(round(b[,a]*scalar), b)  #set of nb models
  newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
  bestrain<-aic.lmx(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
  print(bestrain[1])
}


aic.hydro<-function(y, formula, family, dataset)
{
  m0<-glm(formula, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+cv.depth, family=family, data = dataset)
  m2<-glm(y~log(maxvol)+site+prop.overflow.days, family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+prop.driedout.days, family=family, data = dataset)
  m4<-glm(y~log(maxvol)+site+mean.depth, family=family, data = dataset)
  m5<-glm(y~log(maxvol)+site+long_dry, family=family, data = dataset)
  m6<-glm(y~log(maxvol)+site+last_wet, family=family, data = dataset)
  m7<-glm(y~log(maxvol)+site+change_cv_temp, family=family, data = dataset)
  m8<-glm(y~log(maxvol)+site+change_mean_temp, family=family, data = dataset)
  m17<-glm(y~log(maxvol)+site*cv.depth, family=family, data = dataset)
  m18<-glm(y~log(maxvol)+site*prop.overflow.days, family=family, data = dataset)
  m19<-glm(y~log(maxvol)+site*prop.driedout.days, family=family, data = dataset)
  m20<-glm(y~log(maxvol)+site*mean.depth, family=family, data = dataset)
  m21<-glm(y~log(maxvol)+site*long_dry, family=family, data = dataset)
  m22<-glm(y~log(maxvol)+site*last_wet, family=family, data = dataset)
  m23<-glm(y~log(maxvol)+site*change_cv_temp, family=family, data = dataset)
  m24<-glm(y~log(maxvol)+site*change_mean_temp, family=family, data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m17,m18,m19,m20,m21,m22,m23,m24))
}

aic.hydro.notemp<-function(y, formula, family, dataset)
{
  m0<-glm(formula, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+cv.depth, family=family, data = dataset)
  m2<-glm(y~log(maxvol)+site+prop.overflow.days, family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+prop.driedout.days, family=family, data = dataset)
  m4<-glm(y~log(maxvol)+site+mean.depth, family=family, data = dataset)
  m5<-glm(y~log(maxvol)+site+long_dry, family=family, data = dataset)
  m6<-glm(y~log(maxvol)+site+last_wet, family=family, data = dataset)
  m17<-glm(y~log(maxvol)+site*cv.depth, family=family, data = dataset)
  m18<-glm(y~log(maxvol)+site*prop.overflow.days, family=family, data = dataset)
  m19<-glm(y~log(maxvol)+site*prop.driedout.days, family=family, data = dataset)
  m20<-glm(y~log(maxvol)+site*mean.depth, family=family, data = dataset)
  m21<-glm(y~log(maxvol)+site*long_dry, family=family, data = dataset)
  m22<-glm(y~log(maxvol)+site*last_wet, family=family, data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m17,m18,m19,m20,m21,m22))
}

aic.hydro.pure<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol)+site, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+cv.depth, family=family, data = dataset)
  m2<-glm(y~log(maxvol)+site+prop.overflow.days, family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+prop.driedout.days, family=family, data = dataset)
  m4<-glm(y~log(maxvol)+site+mean.depth, family=family, data = dataset)
  m5<-glm(y~log(maxvol)+site+long_dry, family=family, data = dataset)
  m6<-glm(y~log(maxvol)+site+last_wet, family=family, data = dataset)
  m17<-glm(y~log(maxvol)+site*cv.depth, family=family, data = dataset)
  m18<-glm(y~log(maxvol)+site*prop.overflow.days, family=family, data = dataset)
  m19<-glm(y~log(maxvol)+site*prop.driedout.days, family=family, data = dataset)
  m20<-glm(y~log(maxvol)+site*mean.depth, family=family, data = dataset)
  m21<-glm(y~log(maxvol)+site*long_dry, family=family, data = dataset)
  m22<-glm(y~log(maxvol)+site*last_wet, family=family, data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m17,m18,m19,m20,m21,m22))
}

aic.hydro.add<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol)+site, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+cv.depth, family=family, data = dataset)
  m2<-glm(y~log(maxvol)+site+prop.overflow.days, family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+prop.driedout.days, family=family, data = dataset)
  m4<-glm(y~log(maxvol)+site+mean.depth, family=family, data = dataset)
  m5<-glm(y~log(maxvol)+site+long_dry, family=family, data = dataset)
  m6<-glm(y~log(maxvol)+site+last_wet, family=family, data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6))
}

aic.hydro.nb<-function(y, formula, dataset)
{
  m0<-glm.nb(formula,data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+cv.depth,  data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+prop.overflow.days, data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+prop.driedout.days,  data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+mean.depth,  data = dataset)
  m5<-glm.nb(y~log(maxvol)+site+long_dry,  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site+last_wet,  data = dataset)
  m7<-glm.nb(y~log(maxvol)+site+change_cv_temp,  data = dataset)
  m8<-glm.nb(y~log(maxvol)+site+change_mean_temp,  data = dataset)
  m17<-glm.nb(y~log(maxvol)+site*cv.depth,  data = dataset)
  m18<-glm.nb(y~log(maxvol)+site*prop.overflow.days, data = dataset)
  m19<-glm.nb(y~log(maxvol)+site*prop.driedout.days,  data = dataset)
  m20<-glm.nb(y~log(maxvol)+site*mean.depth,  data = dataset)
  m21<-glm.nb(y~log(maxvol)+site*long_dry,  data = dataset)
  m22<-glm.nb(y~log(maxvol)+site*last_wet,  data = dataset)
  m23<-glm.nb(y~log(maxvol)+site*change_cv_temp,  data = dataset)
  m24<-glm.nb(y~log(maxvol)+site*change_mean_temp,  data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m17,m18,m19,m20,m21,m22,m23,m24))
}

aic.hydro.nb.notemp<-function(y, formula, dataset)
{
  m0<-glm.nb(formula,data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+cv.depth,  data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+prop.overflow.days, data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+prop.driedout.days,  data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+mean.depth,  data = dataset)
  m5<-glm.nb(y~log(maxvol)+site+long_dry,  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site+last_wet,  data = dataset)
  m17<-glm.nb(y~log(maxvol)+site*cv.depth,  data = dataset)
  m18<-glm.nb(y~log(maxvol)+site*prop.overflow.days, data = dataset)
  m19<-glm.nb(y~log(maxvol)+site*prop.driedout.days,  data = dataset)
  m20<-glm.nb(y~log(maxvol)+site*mean.depth,  data = dataset)
  m21<-glm.nb(y~log(maxvol)+site*long_dry,  data = dataset)
  m22<-glm.nb(y~log(maxvol)+site*last_wet,  data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6,  m17,m18,m19,m20,m21,m22))
}


aic.hydro.nb.pure<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site,data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+cv.depth,  data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+prop.overflow.days, data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+prop.driedout.days,  data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+mean.depth,  data = dataset)
  m5<-glm.nb(y~log(maxvol)+site+long_dry,  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site+last_wet,  data = dataset)
  m17<-glm.nb(y~log(maxvol)+site*cv.depth,  data = dataset)
  m18<-glm.nb(y~log(maxvol)+site*prop.overflow.days, data = dataset)
  m19<-glm.nb(y~log(maxvol)+site*prop.driedout.days,  data = dataset)
  m20<-glm.nb(y~log(maxvol)+site*mean.depth,  data = dataset)
  m21<-glm.nb(y~log(maxvol)+site*long_dry,  data = dataset)
  m22<-glm.nb(y~log(maxvol)+site*last_wet,  data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m17,m18,m19,m20,m21,m22))
}

aic.hydro.nb.add<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site,  data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+cv.depth,  data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+prop.overflow.days, data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+prop.driedout.days,  data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+mean.depth,  data = dataset)
  m5<-glm.nb(y~log(maxvol)+site+long_dry,  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site+last_wet,  data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6))
}

aic.hydro.nb.best<-function(a, b, scalar)
{
nbset<-aic.hydro.nb.pure(round(b[,a]*scalar), b)  #set of nb models
newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
besthydro<-aic.hydro.pure(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
print(besthydro[1])
}


aic.site<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol), family=family, data = dataset)
  m1<-glm(y~log(maxvol)+log(mu.scalar), family=family, data = dataset)
  m2<-glm(y~log(maxvol)+log(k.scalar), family=family, data = dataset)
  m3<-glm(y~log(maxvol)+log(mu.scalar)+I(log(mu.scalar)^2), family=family, data = dataset)
  m4<-glm(y~log(maxvol)+log(k.scalar)+I(log(k.scalar)^2), family=family, data = dataset)
  m9<-glm(y~log(maxvol)+log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m10<-glm(y~log(maxvol)+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m11<-glm(y~log(maxvol)+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m12<-glm(y~log(maxvol)+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m9, m10, m11, m12))
}

aic.sitenb<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol), data = dataset)
  m1<-glm.nb(y~log(maxvol)+log(mu.scalar), data = dataset)
  m2<-glm.nb(y~log(maxvol)+log(k.scalar), data = dataset)
  m3<-glm.nb(y~log(maxvol)+log(mu.scalar)+I(log(mu.scalar)^2), data = dataset)
  m4<-glm.nb(y~log(maxvol)+log(k.scalar)+I(log(k.scalar)^2), data = dataset)
  m9<-glm.nb(y~log(maxvol)+log(mu.scalar)*log(k.scalar), data = dataset)
  m10<-glm.nb(y~log(maxvol)+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m11<-glm.nb(y~log(maxvol)+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m12<-glm.nb(y~log(maxvol)+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m9, m10, m11, m12))
}


aic.siterain.nb.best<-function(a, b, scalar)
{
  nbset<-aic.sitenb(round(b[,a]*scalar), b)  #set of nb models
  newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
  best.site.rain<-aic.site(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
  print(best.site.rain[1])
}


aic.site.hydro.nb<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol),  data = dataset)
  m1<-glm.nb(y~log(maxvol)+cv.depth,  data = dataset)
  m2<-glm.nb(y~log(maxvol)+prop.overflow.days, data = dataset)
  m3<-glm.nb(y~log(maxvol)+prop.driedout.days,  data = dataset)
  m4<-glm.nb(y~log(maxvol)+mean.depth,  data = dataset)
  m5<-glm.nb(y~log(maxvol)+long_dry,  data = dataset)
  m6<-glm.nb(y~log(maxvol)+last_wet,  data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6))
}

aic.site.hydro<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol),  family=family, data = dataset)
  m1<-glm(y~log(maxvol)+cv.depth,  family=family, data = dataset)
  m2<-glm(y~log(maxvol)+prop.overflow.days, family=family,data = dataset)
  m3<-glm(y~log(maxvol)+prop.driedout.days,  family=family,data = dataset)
  m4<-glm(y~log(maxvol)+mean.depth, family=family, data = dataset)
  m5<-glm(y~log(maxvol)+long_dry,  family=family,data = dataset)
  m6<-glm(y~log(maxvol)+last_wet,  family=family,data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6))
}

aic.sitehydro.nb.best<-function(a, b, scalar)
{
  nbset<-aic.site.hydro.nb(round(b[,a]*scalar), b)  #set of nb models
  newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
  best.site.hydro<-aic.site.hydro(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
  print(best.site.hydro[1])
}


abs.check<-function(conting, abs){
  m.conting<-conting
  m.abs<-abs
  print(aicset<-model.sel(m.conting,m.abs))
}



#rule: type 3 if sig interactions, type 2 if just sig main effects. Start with type 2 (most power), if sig int switch to type 3

#useful check for NAs

datacheck<-function(a)
  {
  sum(length(which((a)>0)))
}
tapply(fulldata$Odonata_bio, fulldata$site, datacheck)

aic.lmxnb.x<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site, data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+log(mu.scalar), data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+log(k.scalar), data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), data = dataset)
  m5<-glm.nb(y~log(maxvol)+site*log(mu.scalar),  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site*log(k.scalar), data = dataset)
  m7<-glm.nb(y~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m8<-glm.nb(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m9<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar),  data = dataset)
  m10<-glm.nb(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m11<-glm.nb(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m12<-glm.nb(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m13<-glm.nb(y~log(maxvol)+site*log(mu.scalar)*log(k.scalar),  data = dataset)
  m14<-glm.nb(y~log(maxvol)+site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m15<-glm.nb(y~log(maxvol)+site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m16<-glm.nb(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m17<-glm.nb(y~site, data = dataset)
  m18<-glm.nb(y~site+log(mu.scalar), data = dataset)
  m19<-glm.nb(y~site+log(k.scalar), data = dataset)
  m20<-glm.nb(y~site+log(mu.scalar)+I(log(mu.scalar)^2), data = dataset)
  m21<-glm.nb(y~site+log(k.scalar)+I(log(k.scalar)^2), data = dataset)
  m22<-glm.nb(y~site*log(mu.scalar),  data = dataset)
  m23<-glm.nb(y~site*log(k.scalar), data = dataset)
  m24<-glm.nb(y~site*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m25<-glm.nb(y~site*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m26<-glm.nb(y~site+log(mu.scalar)*log(k.scalar),  data = dataset)
  m27<-glm.nb(y~site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m28<-glm.nb(y~site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m29<-glm.nb(y~site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m30<-glm.nb(y~site*log(mu.scalar)*log(k.scalar),  data = dataset)
  m31<-glm.nb(y~site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), data = dataset)
  m32<-glm.nb(y~site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
  m33<-glm.nb(y~site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), data = dataset)
   print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, 
                           m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29, m30, m31, m32, m33))
}

aic.lmx.x<-function(y, family, dataset)
{
  m0<-glm(y~log(maxvol)+site, family=family, data = dataset)
  m1<-glm(y~log(maxvol)+site+log(mu.scalar), family=family, data = dataset)
  m2<-glm(y~log(maxvol)+site+log(k.scalar), family=family, data = dataset)
  m3<-glm(y~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2), family=family, data = dataset)
  m4<-glm(y~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), family=family, data = dataset)
  m5<-glm(y~log(maxvol)+site*log(mu.scalar), family=family, data = dataset)
  m6<-glm(y~log(maxvol)+site*log(k.scalar), family=family, data = dataset)
  m7<-glm(y~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m8<-glm(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m9<-glm(y~log(maxvol)+site+log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m10<-glm(y~log(maxvol)+site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m11<-glm(y~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m12<-glm(y~log(maxvol)+site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m13<-glm(y~log(maxvol)+site*log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m14<-glm(y~log(maxvol)+site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m15<-glm(y~log(maxvol)+site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m16<-glm(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m17<-glm(y~site, family=family, data = dataset)
  m18<-glm(y~site+log(mu.scalar), family=family, data = dataset)
  m19<-glm(y~site+log(k.scalar), family=family, data = dataset)
  m20<-glm(y~site+log(mu.scalar)+I(log(mu.scalar)^2), family=family, data = dataset)
  m21<-glm(y~site+log(k.scalar)+I(log(k.scalar)^2), family=family, data = dataset)
  m22<-glm(y~site*log(mu.scalar), family=family, data = dataset)
  m23<-glm(y~site*log(k.scalar), family=family, data = dataset)
  m24<-glm(y~site*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m25<-glm(y~site*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m26<-glm(y~site+log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m27<-glm(y~site+log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m28<-glm(y~site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m29<-glm(y~site+(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m30<-glm(y~site*log(mu.scalar)*log(k.scalar), family=family, data = dataset)
  m31<-glm(y~log(maxvol)+site*log(mu.scalar)*(log(k.scalar)+I(log(k.scalar)^2)), family=family, data = dataset)
  m32<-glm(y~log(maxvol)+site*log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  m33<-glm(y~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)), family=family, data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16,m17,m18,m19,m20,
                          m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32, m33))
}

aic.hydro.nb.purex<-function(y, dataset)
{
  m0<-glm.nb(y~log(maxvol)+site,data = dataset)
  m1<-glm.nb(y~log(maxvol)+site+cv.depth,  data = dataset)
  m2<-glm.nb(y~log(maxvol)+site+prop.overflow.days, data = dataset)
  m3<-glm.nb(y~log(maxvol)+site+prop.driedout.days,  data = dataset)
  m4<-glm.nb(y~log(maxvol)+site+mean.depth,  data = dataset)
  m5<-glm.nb(y~log(maxvol)+site+long_dry,  data = dataset)
  m6<-glm.nb(y~log(maxvol)+site+last_wet,  data = dataset)
  m17<-glm.nb(y~log(maxvol)+site*cv.depth,  data = dataset)
  m18<-glm.nb(y~log(maxvol)+site*prop.overflow.days, data = dataset)
  m19<-glm.nb(y~log(maxvol)+site*prop.driedout.days,  data = dataset)
  m20<-glm.nb(y~log(maxvol)+site*mean.depth,  data = dataset)
  m21<-glm.nb(y~log(maxvol)+site*long_dry,  data = dataset)
  m22<-glm.nb(y~log(maxvol)+site*last_wet,  data = dataset)
  m23<-glm.nb(y~site,data = dataset)
  m24<-glm.nb(y~site+cv.depth,  data = dataset)
  m25<-glm.nb(y~site+prop.overflow.days, data = dataset)
  m26<-glm.nb(y~site+prop.driedout.days,  data = dataset)
  m27<-glm.nb(y~site+mean.depth,  data = dataset)
  m28<-glm.nb(y~site+long_dry,  data = dataset)
  m29<-glm.nb(y~site+last_wet,  data = dataset)
  m30<-glm.nb(y~site*cv.depth,  data = dataset)
  m31<-glm.nb(y~site*prop.overflow.days, data = dataset)
  m32<-glm.nb(y~site*prop.driedout.days,  data = dataset)
  m33<-glm.nb(y~site*mean.depth,  data = dataset)
  m34<-glm.nb(y~site*long_dry,  data = dataset)
  m35<-glm.nb(y~site*last_wet,  data = dataset)
  print(aicset<-model.sel(m0, m1, m2, m3, m4, m5, m6, m17,m18,m19,m20,m21,m22,
                          m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35))
}


mean.na<-function(y){mean(y,na.rm=TRUE)}
