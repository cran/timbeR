## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, warning = FALSE, message = FALSE, fig.align = 'center',
  out.width = '75%', fig.dim = c(8,5)
)

## -----------------------------------------------------------------------------
# install.packages("devtools")
# options(download.file.method = "libcurl")
# devtools::install_github('sergiocostafh/timbeR')

library(dplyr)
library(timbeR)

glimpse(tree_scaling)

## -----------------------------------------------------------------------------
library(ggplot2)

tree_scaling <- tree_scaling %>% 
  mutate(did = di/dbh,
         hih = hi/h)

ggplot(tree_scaling, aes(x = hih, y = did, group = tree_id))+
  geom_point()+
  labs(x = 'hi / h',
       y = 'di / dbh')

## -----------------------------------------------------------------------------
poli5 <- lm(did~hih+I(hih^2)+I(hih^3)+I(hih^4)+I(hih^5),tree_scaling)
summary(poli5)

tree_scaling <- tree_scaling %>% 
  mutate(di_poli = predict(poli5)*dbh)

poli_rmse <- tree_scaling %>% 
  summarise(RMSE = sqrt(sum((di_poli-di)^2)/mean(di_poli))) %>% 
  pull(RMSE) %>% 
  round(2)

ggplot(tree_scaling,aes(x=hih))+
  geom_point(aes(y = (di_poli-di)/di_poli*100))+
  geom_hline(aes(yintercept = 0))+
  scale_y_continuous(limits=c(-60,60), breaks = seq(-100,100,20))+
  scale_x_continuous(limits=c(0,1))+
  labs(x = 'hi / h', y = 'Residuals (%)',
       title = '5th degree polynomial taper function (Schöepfer, 1966)',
       subtitle = 'Dispersion of residuals along the stem',
       caption = paste0('Root Mean Squared Error = ', poli_rmse,'%'))+
  theme(plot.title.position = 'plot')

## -----------------------------------------------------------------------------
library(minpack.lm)

bi <-  nlsLM(di ~ dbh*(log(sin((pi/2)*(hih)))/(log(sin((pi/2)*(1.3/h)))))**
             (b0+b1*sin((pi/2)*(hih))+b2*cos((3*pi/2)*(hih))+b3*(sin((pi/2)*(hih))/(hih))+
                b4*dbh+b5*(hih)*dbh**0.5+b6*(hih)*h**0.5),
           data=tree_scaling,
           start=list(b0=1.8,b1=-0.2,b2=-0.04,b3=-0.9,b4=-0.0006,b5=0.07,b6=-.14))
summary(bi)

tree_scaling <- tree_scaling %>% 
  mutate(di_bi = predict(bi))

bi_rmse <- tree_scaling %>% 
  summarise(RMSE = sqrt(sum((di_bi-di)^2)/mean(di_bi))) %>% 
  pull(RMSE) %>% 
  round(2)

ggplot(tree_scaling,aes(x=hih))+
  geom_point(aes(y = (di_bi-di)/di_bi*100))+
  geom_hline(aes(yintercept = 0))+
  scale_y_continuous(limits=c(-60,60), breaks = seq(-100,100,20))+
  scale_x_continuous(limits=c(0,1))+
  labs(x = 'hi / h', y = 'Residuals (%)',
       title = 'Bi (2000) trigonometric variable-form taper function',
       subtitle = 'Dispersion of residuals along the stem',
       caption = paste0('Root Mean Squared Error = ', bi_rmse,'%'))+
  theme(plot.title.position = 'plot')


## -----------------------------------------------------------------------------
kozak <- nlsLM(di ~
                 b0*(dbh**b1)*(h**b2)*((1-hih**(1/4))/(1-(p^(1/3))))**(b3*hih**4+b4*(1/exp(dbh/h))+b5*((1-hih**(1/4))/(1-(p^(1/3))))**0.1+b6*
                                                                         (1/dbh)+b7*(h**(1-
                                                                                           hih**(1/3)))+b8*((1-hih**(1/4))/(1-(p^(1/3))))),
               start=list(b0=1.00,b1=.97,b2=.03,b3=.49,b4=-
                            0.87,b5=0.50,b6=3.88,b7=0.03,b8=-0.19, p = .1),
               data = tree_scaling,
               control = nls.lm.control(maxiter = 1000, maxfev = 2000)
)
summary(kozak)

tree_scaling <- tree_scaling %>% 
  mutate(di_kozak = predict(kozak))

kozak_rmse <- tree_scaling %>% 
  summarise(RMSE = sqrt(sum((di_kozak-di)^2)/mean(di_kozak))) %>% 
  pull(RMSE) %>% 
  round(2)

ggplot(tree_scaling, aes(x=hih))+
  geom_point(aes(y = (di_kozak-di)/di_kozak*100))+
  geom_hline(aes(yintercept = 0))+
  scale_y_continuous(limits=c(-100,100), breaks = seq(-100,100,20))+
  scale_x_continuous(limits=c(0,1))+
  labs(x = 'hi / h', y = 'Residuals (%)',
       title = 'Kozak (2004) variable-form taper function',
       subtitle = 'Dispersion of residuals along the stem',
       caption = paste0('Root Mean Squared Error = ', kozak_rmse,'%'))+
  theme(plot.title.position = 'plot')


## -----------------------------------------------------------------------------
dbh <- 25
h <- 20

## -----------------------------------------------------------------------------
coef_poli <- coef(poli5)
coef_bi <- coef(bi)
coef_kozak <- coef(kozak)[-10]
p_kozak <- coef(kozak)[10]

## -----------------------------------------------------------------------------
hi <- 15

poly5_di(dbh, h, hi, coef_poli)
bi_di(dbh, h, hi, coef_bi)
kozak_di(dbh, h, hi, coef_kozak, p = p_kozak)


## -----------------------------------------------------------------------------
hi <- seq(0.1,h,.1)

ggplot(mapping=aes(x=hi))+
  geom_line(aes(y=poly5_di(dbh, h, hi, coef_poli), linetype = '5th degree polynomial'))+
  geom_line(aes(y=bi_di(dbh, h, hi, coef_bi), linetype = 'Bi (2000)'))+
  geom_line(aes(y=kozak_di(dbh, h, hi, coef_kozak, p_kozak), linetype = 'Kozak (2004)'))+
  scale_linetype_manual(name = 'Fitted models', values = c('solid','dashed','dotted'))+
  labs(x = 'hi (m)',
       y = 'Predicted di (cm)')


## -----------------------------------------------------------------------------
di <- 10

poly5_hi(dbh, h, di, coef_poli)
bi_hi(dbh, h, di, coef_bi)
kozak_hi(dbh, h, di, coef_kozak, p_kozak)

## -----------------------------------------------------------------------------
poly5_vol(dbh, h, coef_poli)
bi_vol(dbh, h, coef_bi)
kozak_vol(dbh, h, coef_kozak, p_kozak)

## -----------------------------------------------------------------------------
hi = 15
h0 = .2

poly5_vol(dbh, h, coef_poli, hi, h0)
bi_vol(dbh, h, coef_bi, hi, h0)
kozak_vol(dbh, h, coef_kozak, p_kozak, hi, h0)

## -----------------------------------------------------------------------------
assortments <- data.frame(
  NAME = c('15-25','4-15'),
  SED = c(15,4),
  MINLENGTH = c(2.65,2),
  MAXLENGTH = c(2.65,4.2),
  LOSS = c(5,5)
)


## -----------------------------------------------------------------------------
poly5_logs(dbh, h, coef_poli, assortments)
bi_logs(dbh, h, coef_bi, assortments)
kozak_logs(dbh, h, coef_kozak, p_kozak, assortments)

## -----------------------------------------------------------------------------
poly5_logs_plot(dbh, h, coef_poli, assortments)
bi_logs_plot(dbh, h, coef_bi, assortments)
kozak_logs_plot(dbh, h, coef_kozak, p_kozak, assortments)

## -----------------------------------------------------------------------------
# Using mapply

tree_data <- data.frame(dbh = c(18.3, 23.7, 27.2, 24.5, 20, 19.7),
                   h = c(18, 24, 28, 24, 18.5, 19.2))

assortment_vol <- mapply(
  poly5_logs,
  dbh = tree_data$dbh,
  h = tree_data$h,
  SIMPLIFY = T,
  MoreArgs = list(
    coef = coef_poli,
    assortments = assortments,
    stump_height = 0.2,
    total_volume = T,
    only_vol = T
  )
) %>%
  t()


assortment_vol


# Bindind tree_data and volumes output

library(tidyr)

cbind(tree_data, assortment_vol) %>% 
  unnest()

## -----------------------------------------------------------------------------
# Using pmap

library(purrr)

tree_data %>% 
  mutate(coef = list(coef_poli),
         assortments = list(assortments),
         stump_height = 0.2,
         total_volume = T,
         only_vol = T) %>% 
  mutate(assortment_vol = pmap(.,poly5_logs)) %>% 
  select(dbh, h, assortment_vol) %>% 
  unnest(assortment_vol)

