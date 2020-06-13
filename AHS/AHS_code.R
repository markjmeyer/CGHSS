###### lag analysis ######
# Code to accompany Boyce et al. (2020):
# Financial assistance for health security: effects of international financial
# 	assistance on capacities for preventing, detecting, and responding to
# 	public health emergencies

### source BPPR code ###
source('bppr.R')

### load data files ###
comLag    <- readRDS('comLag.RData')
legLag    <- readRDS('legLag.RData')
zooLag    <- readRDS('zooLag.RData')
labLag    <- readRDS('labLag.RData')
survLag   <- readRDS('survLag.RData')
hrLag     <- readRDS('hrLag.RData')
prepLag   <- readRDS('prepLag.RData')
riskLag   <- readRDS('riskLag.RData')
eresLag   <- readRDS('eresLag.RData')
coorLag   <- readRDS('coorLag.RData')
amrLag    <- readRDS('amrLag.RData')
foodLag   <- readRDS('foodLag.RData')
bioLag    <- readRDS('bioLag.RData')


##### overall capacity score #####
modComb <- bppr(fixed = comOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = comLag, up = 1000, nu = 0)
sum(abs(modComb$geweke) > 3)

round(cbind(t(apply(modComb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))), 
            apply(modComb$alphas > 0, 2, mean)), 4)



##### legislation capacities #####
modLegb <- bppr(fixed = legOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = legLag, up = 1000)
sum(abs(modLegb$geweke) > 3)

round(cbind(t(apply(modLegb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modLegb$alphas > 0, 2, mean)), 4)



##### zoonotic capacities #####
modZoob <- bppr(fixed = zooOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = zooLag, up = 1000)
sum(abs(modZoob$geweke) > 3)

round(cbind(t(apply(modZoob$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modZoob$alphas > 0, 2, mean)), 4)



##### laboratory capacities #####
modLabb <- bppr(fixed = labOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = labLag, up = 1000)
sum(abs(modLabb$geweke) > 3)

round(cbind(t(apply(modLabb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modLabb$alphas > 0, 2, mean)), 4)



##### surveillance capacities #####
modSurb <- bppr(fixed = survOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = survLag, up = 1000)
sum(abs(modSurb$geweke) > 3)

round(cbind(t(apply(modSurb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))), 
            apply(modSurb$alphas > 0, 2, mean)), 4)



##### workforce capacities #####
modHrb <- bppr(fixed = hrOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
               data = hrLag, up = 1000)
sum(abs(modHrb$geweke) > 3)

round(cbind(t(apply(modHrb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modHrb$alphas > 0, 2, mean)), 4)



##### preparedness capacities #####
modPreb <- bppr(fixed = preOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = prepLag, up = 1000)
sum(abs(modPreb$geweke) > 3)

round(cbind(t(apply(modPreb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))), 
            apply(modPreb$alphas > 0, 2, mean)), 4)



##### risk communication capacities #####
modRisb <- bppr(fixed = commOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = riskLag, up = 1000)
sum(abs(modRisb$geweke) > 3)

round(cbind(t(apply(modRisb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modRisb$alphas > 0, 2, mean)), 4)



##### emergency response capacities #####
modEreb <- bppr(fixed = respOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = eresLag, up = 1000)
sum(abs(modEreb$geweke) > 3)

round(cbind(t(apply(modEreb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modEreb$alphas > 0, 2, mean)), 4)



##### coordination capacities #####
modCoob <- bppr(fixed = coorOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = coorLag, up = 1000)
sum(abs(modCoob$geweke) > 3)

round(cbind(t(apply(modCoob$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modCoob$alphas > 0, 2, mean)), 4)



##### AMR capacities #####
modAmrb <- bppr(fixed = respOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = amrLag, up = 1000)
sum(abs(modAmrb$geweke) > 3)

round(cbind(t(apply(modAmrb$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modAmrb$alphas > 0, 2, mean)), 4)



##### food safety capacities #####
# THREE LEVEL, NOT FOUR #
modFoob <- bppr(fixed = foodOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = foodLag, up = 1000, cut.start = c(0, 0.5))
sum(abs(modFoob$geweke) > 3)

round(cbind(t(apply(modFoob$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modFoob$alphas > 0, 2, mean)), 4)



##### biosafety capacities #####
modBiob <- bppr(fixed = labOrd ~ pop0 + gdp0 + factor(year) + factor(region), random = ~ lag0 + lag1 + lag2,
                data = bioLag, up = 1000)
sum(abs(modBiob$geweke) > 3)

round(cbind(t(apply(modBiob$alphas, 2, quantile, probs = c(0.5, 0.025, 0.975))),
            apply(modBiob$alphas > 0, 2, mean)), 4)


###### Marginal Effects ######

##### legistation capacity #####
# lag 1 #
summary(legLag$lag1)
nx    <- 100
mmes_leg  <- hmes_leg   <- lmes_leg   <- matrix(0, nrow = nx, ncol = 3)
for(x in 1:nx){
  mex   <- me(eff = x-1, model = modLegb, data = legLag)
  
  mmes_leg[x,]  <- mex$mme
  lmes_leg[x,]  <- mex$lme
  hmes_leg[x,]  <- mex$hme
}

# summary estimates #
legSumM  <- me(eff = median(legLag$lag1), model = modLegb, data = legLag)
c(legSumM$mme[2], legSumM$lme[2], legSumM$hme[2])

legSum1  <- me(eff = quantile(legLag$lag1, probs = 0.25), model = modLegb, data = legLag)
c(legSum1$mme[2], legSum1$lme[2], legSum1$hme[2])

legSum3  <- me(eff = quantile(legLag$lag1, probs = 0.75), model = modLegb, data = legLag)
c(legSum3$mme[2], legSum3$lme[2], legSum3$hme[2])



##### laboratory capacity #####
# lag 1 #
summary(labLag$lag1)
nx    <- 50
mmes_lab  <- hmes_lab   <- lmes_lab   <- matrix(0, nrow = nx, ncol = 3)
for(x in 1:nx){
  mex   <- me(eff = x-1, model = modLabb, data = labLag)
  
  mmes_lab[x,]  <- mex$mme
  lmes_lab[x,]  <- mex$lme
  hmes_lab[x,]  <- mex$hme
}

# summary estimates #
labSumM  <- me(eff = median(labLag$lag1), model = modLabb, data = labLag)
c(labSumM$mme[2], labSumM$lme[2], labSumM$hme[2])

labSum1  <- me(eff = quantile(labLag$lag1, probs = 0.25), model = modLabb, data = labLag)
c(labSum1$mme[2], labSum1$lme[2], labSum1$hme[2])

labSum3  <- me(eff = quantile(labLag$lag1, probs = 0.75), model = modLabb, data = labLag)
c(labSum3$mme[2], labSum3$lme[2], labSum3$hme[2])



##### workforce capacity #####
# lag 1 #
summary(hrLag$lag1)
nx    <- 50
mmes_hr  <- hmes_hr   <- lmes_hr   <- matrix(0, nrow = nx, ncol = 3)
for(x in 1:nx){
  mex   <- me(eff = x-1, model = modHrb, data = hrLag)
  
  mmes_hr[x,]  <- mex$mme
  lmes_hr[x,]  <- mex$lme
  hmes_hr[x,]  <- mex$hme
}

# summary estimates #
hrSumM  <- me(eff = median(hrLag$lag1), model = modHrb, data = hrLag)
c(hrSumM$mme[2], hrSumM$lme[2], hrSumM$hme[2])

hrSum1  <- me(eff = quantile(hrLag$lag1, probs = 0.25), model = modHrb, data = hrLag)
c(hrSum1$mme[2], hrSum1$lme[2], hrSum1$hme[2])

hrSum3  <- me(eff = quantile(hrLag$lag1, probs = 0.75), model = modHrb, data = hrLag)
c(hrSum3$mme[2], hrSum3$lme[2], hrSum3$hme[2])



##### risk communication capacity #####
# lag 1 #
summary(riskLag$lag1)
nx    <- 20
mmes_ris  <- hmes_ris   <- lmes_ris   <- matrix(0, nrow = nx, ncol = 3)
for(x in 1:nx){
  mex   <- me(eff = x-1, model = modRisb, data = riskLag)
  
  mmes_ris[x,]  <- mex$mme
  lmes_ris[x,]  <- mex$lme
  hmes_ris[x,]  <- mex$hme
}

# summary estimates #
risSumM  <- me(eff = median(riskLag$lag1), model = modRisb, data = riskLag)
c(risSumM$mme[2], risSumM$lme[2], risSumM$hme[2])

risSum1  <- me(eff = quantile(riskLag$lag1, probs = 0.25), model = modRisb, data = riskLag)
c(risSum1$mme[2], risSum1$lme[2], risSum1$hme[2])

risSum3  <- me(eff = quantile(riskLag$lag1, probs = 0.75), model = modRisb, data = riskLag)
c(risSum3$mme[2], risSum3$lme[2], risSum3$hme[2])



##### Graphics #####

# pdf('sig_capacity_lags.pdf')
par(mfrow = c(2,2))

xp  <- 41
matplot(0:(xp-1), cbind(mmes_leg[1:xp,2], lmes_leg[1:xp,2], hmes_leg[1:xp,2]), type = 'n', ylab = '',
        xlab = 'Disbursed Funds (in $1m)', main = 'Legislation Capacity\nLag of One Year')
title(ylab = 'Adjusted Probability\nof Level Change', line = 2)
abline(v = seq(0, xp, by = 10), h = axTicks(2), col = rgb(0.75, 0.75, 0.75, alpha = 0.5), lty = 3)
polygon(c(0:(xp-1), (xp-1):0), c(lmes_leg[1:xp,2], rev(hmes_leg[1:xp,2])), 
        col = rgb(0/255, 61/255, 165/255, alpha = 0.5), border = NA)
lines(0:(xp-1), mmes_leg[1:xp,2], col = rgb(0/255, 61/255, 165/255), lwd = 2)

xp  <- 16
matplot(0:(xp-1), cbind(mmes_lab[1:xp,2], lmes_lab[1:xp,2], hmes_lab[1:xp,2]), type = 'n', ylab = '',
        xlab = 'Disbursed Funds (in $1m)', main = 'Laboratory Capacity\nLag of One Year')
title(ylab = 'Adjusted Probability\nof Level Change', line = 2)
abline(v = seq(0, xp, by = 5), h = axTicks(2), col = rgb(0.75, 0.75, 0.75, alpha = 0.5), lty = 3)
polygon(c(0:(xp-1), (xp-1):0), c(lmes_lab[1:xp,2], rev(hmes_lab[1:xp,2])), col = rgb(0/255, 61/255, 165/255, alpha = 0.5), 
        border = NA)
lines(0:(xp-1), mmes_lab[1:xp,2], col = rgb(0/255, 61/255, 165/255), lwd = 2)

xp  <- 41
matplot(0:(xp-1), cbind(mmes_hr[1:xp,2], lmes_hr[1:xp,2], hmes_hr[1:xp,2]), type = 'n', ylab = '',
        xlab = 'Disbursed Funds (in $1m)', main = 'Workforce Capacity\nLag of One Year')
title(ylab = 'Adjusted Probability\nof Level Change', line = 2)
abline(v = seq(0, xp, by = 10), h = axTicks(2), col = rgb(0.75, 0.75, 0.75, alpha = 0.5), lty = 3)
polygon(c(0:(xp-1), (xp-1):0), c(lmes_hr[1:xp,2], rev(hmes_hr[1:xp,2])), col = rgb(0/255, 61/255, 165/255, alpha = 0.5), 
        border = NA)
lines(0:(xp-1), mmes_hr[1:xp,2], col = rgb(0/255, 61/255, 165/255), lwd = 2)

xp  <- 16
matplot(0:(xp-1), cbind(mmes_ris[1:xp,2], lmes_ris[1:xp,2], hmes_ris[1:xp,2]), type = 'n', ylab = '',
        xlab = 'Disbursed Funds (in $1m)', main = 'Risk Communication Capacity\nLag of One Year')
title(ylab = 'Adjusted Probability\nof Level Change', line = 2)
abline(v = seq(0, xp, by = 5), h = axTicks(2), col = rgb(0.75, 0.75, 0.75, alpha = 0.5), lty = 3)
polygon(c(0:(xp-1), (xp-1):0), c(lmes_ris[1:xp,2], rev(hmes_ris[1:xp,2])), col = rgb(0/255, 61/255, 165/255, alpha = 0.5), 
        border = NA)
lines(0:(xp-1), mmes_ris[1:xp,2], col = rgb(0/255, 61/255, 165/255), lwd = 2)
# dev.off()




