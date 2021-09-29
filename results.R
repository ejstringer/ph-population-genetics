#.------------------------------------------------------------------------------
# analysis for chapter 1: temporal population genetics


# libraries --------------------------------------------------------------------
library(tidyverse)
library(dartR)
library(emR)
library(gt) # tables
library(lme4) # models
library(lubridate)

# load data --------------------------------------------------------------------
ph <- emR::em.data.ph_filtered
booms <- emR::em.data.peaks %>%
  rename(sinceBoom = sincePeak.adj, boomNo = peakNo.adj) %>% 
  select(trip, rain, npp, capEstimates, sinceBoom, boomNo) %>% 
  mutate(boomNo = plyr::revalue(boomNo, c("peak2" = "boom2007",
                                          "peak3" = "boom2010",
                                          "peak4" = "boom2016"))) %>% 
  filter(lubridate::ymd(trip) >= lubridate::ymd("2006-01-01") & 
         lubridate::ymd(trip) <= lubridate::ymd("2017-05-01"))
fst <- emR::em.data.ph_fst_grids[[4]] %>% 
  filter(peakNo.adj != "peak1") %>% 
  select(Population1, Population2, pairs, Fst, fst.log, 
         weights, n1, n2, trip, rain, npp, capEstimates,
         sincePeak.adj, peakNo.adj, metres) %>% 
  rename(sinceBoom = sincePeak.adj, boomNo = peakNo.adj) %>% 
  mutate(sinceBoom.log = log(sinceBoom), 
         boomNo = plyr::revalue(boomNo, c("peak2" = "boom2007",
                                          "peak3" = "boom2010",
                                          "peak4" = "boom2016")))
  
he <- emR::em.data.ph_he_grids[[4]] %>% 
  filter(peakNo.adj != "peak1") %>% 
  select(He, varHe, weights, trip, gridId, n, rain, npp, capEstimates,
         sincePeak.adj, peakNo.adj) %>% 
  rename(sinceBoom = sincePeak.adj, boomNo = peakNo.adj) %>% 
  mutate(sinceBoom.log = log(sinceBoom),
         boomNo = plyr::revalue(boomNo, c("peak2" = "boom2007",
                                          "peak3" = "boom2010",
                                          "peak4" = "boom2016")))  




# rainfall and captures --------------------------------------------------------

panDates <- ymd(booms$trip[booms$sinceBoom == 0])

ggplot(booms, aes(ymd(trip), rain)) +
  geom_area(aes(y = rain), fill = "blue", alpha = 0.5)+
  geom_area(aes(y = capEstimates * 10), lwd = 1.1, alpha = 0.3)+
  geom_vline(xintercept = panDates, size = 2, colour = "red",
             linetype = 2, alpha=0.9)+
  theme_classic()+
  xlab("time") +
  ylab("rain (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y") 


# simulation -------------------------------------------------------------------
res <- em.data.simBoombust_dispersal[[6]]

 ggplot(res, aes(x=gen, y=fst))+
  geom_line()+geom_point()+
  theme_classic()+
  xlab("generations") +
  ylab("Fst")+
  ylim(c(-0.01,0.20))
  


 ggplot(res, aes(x=gen, y=he))+
  geom_line()+geom_point()+
  theme_classic()+
  xlab("generations") +
  ylab("He")+
  ylim(c(0.08,0.15))

# map --------------------------------------------------------------------------
 ph@other$latlong <- ph@other$latlon
 gl.map.interactive(ph)
 
# sample -----------------------------------------------------------------------
 
 sam <- rbind(fst[, c("Population1", "n1", "trip", "boomNo")],
              fst[, c("Population2", "n2", "trip", "boomNo")] %>% 
                rename(Population1 = Population2, n1 = n2)) %>% 
   unique() %>% select(trip, Population1, n1, boomNo)
 names(sam) <- c("trip", "gridId", "sample size", "boomNo")
 
 sample.table <-function(sam, peak, header, subheader){
   tb <- sam %>% filter(boomNo == peak) %>% 
     select(-boomNo) %>%  
     group_by(trip) %>% 
     gt(rowname_col = "gridId") %>% 
     tab_header(
       title = md(paste0("**", header, "**")),
       subtitle = html(paste0("<em>",subheader, "</em>"))
     )
   return(tb) 
 }
 
 sample.table(sam, "boom2007", "2007 Boom", "September 2007 to April 2010")
 sample.table(sam, "boom2010", "2010 Boom", "November 2010 to April 2015")
 sample.table(sam, "boom2016", "2016 Boom", "April 2016 to May 2017")
 


# fst plots and models ---------------------------------------------------------

ggplot(fst, aes(sinceBoom.log, fst.log)) +
  geom_jitter(colour = "#f8766d") + 
  facet_grid(~boomNo)+
  theme_bw() + 
  xlab("log months since boom") +
  ylab("log Fst")+
  geom_smooth(method = "lm") 

boxplot(fst$fst.log)
hist(fst$fst.log)

# fit separate models for each boom with site as a random effect and accounting for distance = log(dist)
m1 <- lmer(fst.log ~ sinceBoom.log + log(metres) + (1|pairs), data = filter(fst, boomNo == "boom2007"))
m2 <- lmer(fst.log ~ sinceBoom.log + log(metres) + (1|pairs), data = filter(fst, boomNo == "boom2010"))
m3 <- lmer(fst.log ~ sinceBoom.log + log(metres) + (1|pairs), data = filter(fst, boomNo == "boom2016"))

sjPlot::tab_model(m1,m2,m3,
                  dv.labels = c("2007 Boom", "2010 Boom", "2016 Boom"), collapse.ci = F)

              # fit as one model !!! needs work
              fst$ran <- paste(fst$pairs, fst$boomNo)
              m <- lmer(fst.log ~ sinceBoom.log + log(metres) + boomNo + (1|ran), data = fst)
              
              sjPlot::tab_model(m, dv.labels = "full model", collapse.ci = F)

# he plots and models ----------------------------------------------------------
ggplot(he, aes(sinceBoom.log, He))+
  geom_jitter(colour = "#f8766d") + 
  facet_grid(~boomNo)+
  theme_bw() + 
  xlab("log months since boom") +
  ylab("exp heterozygosity")+
  geom_smooth(method = "lm")

boxplot(he$He)
hist(scale(he$He))

mm1 <- lm(He ~ sinceBoom.log, data = filter(he, boomNo == "boom2007"), weights = weights)
mm2 <- lm(He ~ sinceBoom.log, data = filter(he, boomNo == "boom2010"), weights = weights)
mm3 <- lm(He ~ sinceBoom.log, data = filter(he, boomNo == "boom2016"), weights = weights)
sjPlot::tab_model(mm1,mm2,mm3,
                  dv.labels = c("2007 Boom", "2010 Boom", "2016 Boom"), collapse.ci = F)
mm <-  lm(He ~ sinceBoom.log + boomNo,data = he , weights = weights)
sjPlot::tab_model(mm, dv.labels = "full model", collapse.ci = F)


# predicted values -------------------------------------------------------------



