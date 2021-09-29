#.------------------------------------------------------------------------------
# analysis for chapter 1: temporal population genetics


# libraries --------------------------------------------------------------------
library(tidyverse)
library(dartR)
library(emR)
library(gt) # tables
library(lme4) # models

# load data --------------------------------------------------------------------
ph <- emR::em.data.ph_filtered
peaks <- emR::em.data.peaks
fst <- emR::em.data.ph_fst_grids[[4]] %>% filter(peakNo != "peak1")# min sample size 4
he <- emR::em.data.ph_he_grids[[4]] %>% filter(peakNo != "peak1")
glimpse(peaks)
glimpse(fst)
glimpse(he)
boxplot(fst$fst.log)
boxplot(he$He)
# map --------------------------------------------------------------------------
ph@other$latlong <- ph@other$latlon
gl.map.interactive(ph)

# sample -----------------------------------------------------------------------
  
sam <- rbind(fst[, c("Population1", "n1", "trip", "peakNo")],
           fst[, c("Population2", "n2", "trip", "peakNo")] %>% 
        rename(Population1 = Population2, n1 = n2)) %>% 
  unique() %>% select(trip, Population1, n1, peakNo)
names(sam) <- c("trip", "gridId", "sample size", "peakNo")

sample.table <-function(sam, peak, header, subheader){
  tb <- sam %>% filter(peakNo == peak) %>% 
    select(-peakNo) %>%  
    group_by(trip) %>% 
    gt(rowname_col = "gridId") %>% 
    tab_header(
      title = md(paste0("**", header, "**")),
      subtitle = html(paste0("<em>",subheader, "</em>"))
  )
 return(tb) 
}

sample.table(sam, "peak2", "2007 Boom", "May 2007 to September 2009")
sample.table(sam, "peak3", "2010 Boom", "April 2010 to April 2015")
sample.table(sam, "peak4", "2016 Boom", "April 2016 to May 2017")

# rainfall and captures --------------------------------------------------------
mypeaks <- filter(peaks, 
                  lubridate::ymd(trip) >= lubridate::ymd("2006-01-01") & 
                  lubridate::ymd(trip) <= lubridate::ymd("2017-05-01"))
mypeaks$peakNo <- mypeaks$peakNo.adj
mypeaks$sincePeak <- mypeaks$sincePeak.adj


panmictic <- mypeaks$sincePeak == 0 & complete.cases(mypeaks$peakNo)
panDates <- lubridate::ymd(mypeaks$trip[panmictic])
df <- data.frame(panDates, y1 = 0, y2 = 50)
rp <- ggplot(mypeaks, aes(lubridate::ymd(trip), rain)) + 
  geom_area(aes(y = rain), fill = "blue", alpha = 0.5) + 
  geom_area(data = filter(mypeaks, complete.cases(capEstimates)), 
            aes(lubridate::ymd(trip), capEstimates * 10), lwd = 1.1, 
            alpha = 0.3) + theme_classic() + xlab("time") + 
  ylab("rain (mm)") +
  geom_vline(xintercept = panDates, size = 2, 
             colour = "red", linetype = 2, alpha=0.9)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y") 

rp

# simulation -------------------------------------------------------------------
res <- em.data.simBoombust_dispersal[[6]]
glimpse(res)

 ggplot(res, aes(x=gen, y=fst))+
  geom_line()+geom_point()+
  theme_classic()+
  xlab("generations") +
  ylab("Fst")
  # ylim(c(-0.01,0.20))+
  # xlim(c(-4, 20))


 ggplot(res, aes(x=gen, y=he))+
  geom_line()+geom_point()+
  theme_classic()+
  xlab("generations") +
  ylab("He")

# exploratory plots ------------------------------------------------------------


 ggplot(fst, aes(log(fst$sincePeak.adj), fst.log)) +
   geom_jitter(colour = "#f8766d") + 
   facet_grid(~peakNo.adj)+
   theme_bw() + 
   xlab("log months since peak") +
   ylab("log Fst")+
   geom_smooth(method = "lm") 

  ggplot(he, aes(log(sincePeak.adj), He))+
    geom_jitter(colour = "#f8766d") + 
    facet_grid(~peakNo.adj)+
    theme_bw() + 
    xlab("log months since peak") +
    ylab("exp heterozygosity")+
    geom_smooth(method = "lm")
   
boxplot(he$he)
# models -----------------------------------------------------------------------
fst$log.sincePeak <- log(fst$sincePeak.adj)
fst$ran <- paste(fst$pairs, fst$peakNo.adj)

hist(fst$fst.log)
# fit separate models for each peak with site as a random effect and accounting for distance = log(dist)
m1 <- lmer(fst.log ~ log.sincePeak + log(metres) + (1|pairs), data = filter(fst, peakNo == "peak2"))
m2 <- lmer(fst.log ~ log.sincePeak + log(metres) + (1|pairs), data = filter(fst, peakNo == "peak3"))
m3 <- lmer(fst.log ~ log.sincePeak + log(metres) + (1|pairs), data = filter(fst, peakNo == "peak4" & fst.log > -8))

m <- lmer(fst.log ~ log.sincePeak + peakNo.adj + capEstimates + (1|ran), data = fst)
isSingular(m)

summary(m1)
summary(m2)
summary(m3)

sjPlot::tab_model(m1,m2,m3,
                  dv.labels = c("boom1", "boom2", "boom3"), collapse.ci = T)
sjPlot::tab_model(m, dv.labels = "full model", collapse.ci = T)

he$log.sincePeak <- log(he$sincePeak.adj)
he$he <- he$He
hist(scale(he$He))

m1 <- lm(He ~ log.sincePeak, data = filter(he, peakNo == "peak2"), weights = weights)
m2 <- lm(He ~ log.sincePeak, data = filter(he, peakNo == "peak3"), weights = weights)
m3 <- lm(He ~ log.sincePeak, data = filter(he, peakNo == "peak4"), weights = weights)
sjPlot::tab_model(m1,m2,m3,
                  dv.labels = c("boom1", "boom2", "boom3"), collapse.ci = T)



# predicted values -------------------------------------------------------------



