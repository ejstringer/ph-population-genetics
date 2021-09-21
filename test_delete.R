

fis <- emR::em.data.ph_fis_grids

library(shiny)
ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      actionButton("switch_tab", "Go to the second tab")
    ),
    
    
    mainPanel(
    navlistPanel(id="inTabset", "right",
    tabPanel("Tab 1",actionButton("switch_tab2", "Go to the third tab")
    ),
    tabPanel("Tab 2", tableOutput("table")),
    tabPanel("Tab 3", "there!"))
    )
)
)

server <- function(input, output, session) {
  
  observeEvent(input$switch_tab, {
    updateTabsetPanel(session, "inTabset",selected = "Tab 2")
  })
  observeEvent(input$switch_tab2, {
    updateTabsetPanel(session, "inTabset",selected = "Tab 3")
  })
  
  output$table <- renderTable({
    head(fis[[1]])
  })
  
}
shinyApp(ui = ui, server = server)


fis <- emR::em.data.ph_fis_grids[[4]]
fst <- emR::em.data.ph_fst_grids[[4]]
library(tidyverse)
ggplot(filter(fis, peakNo == "peak3"),
       aes((He))) + geom_histogram()
hist(log(fis[fis$He > 0.12,]$He))
ggplot(filter(fst, peakNo != "peak1"), aes(capEstimates, Fst)) +
  geom_point() +
  facet_grid(~peakNo) +
  geom_smooth(method = "lm")
tapply(fis$x, fis$peakNo, length)
pherm <- emR::em.data.ph_filtered
lm(He ~ capEstimates, filter(fis, peakNo == "peak2"), weights = 1/varHe)%>% summary
lm(He ~ capEstimates, filter(fis, peakNo == "peak3"), weights = 1/varHe)%>% summary
lm(He ~ capEstimates, filter(fis, peakNo == "peak4"), weights = 1/varHe)%>% summary




lm(He ~ capEstimates, fis, weights = 1/varHe) %>% summary

tripdate <- pherm@other$ind.metrics$trip
table(tripdate) %>% table
index <- table(tripdate) < 6
levels(tripdate) == names(index)
tooSmall <- levels(tripdate)[index]
library(dartR)
pherm2 <- gl.drop.pop(pherm, pop.list = tooSmall, as.pop = "trip")
class(pherm@other$ind.metrics$trip)
pop(pherm) <- pherm@other$ind.metrics$trip  

pherm@pop %>% table %>% length

fstTrip <- gl.fst.pop(pherm)

ggplot(fst, aes("fst", log(fstScaled))) +
  geom_boxplot( outlier.colour = "red")


ggplot(filter(fst, log(fstScaled) > -4),
       aes(xlog, log(fstScaled))) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_grid(~peakNo)
min(fst$ylog)
lm(log(fstScaled) ~ xlog + peakNo, filter(fst, log(fstScaled) > -4)) %>% 
  summary

lm(log(Fst) ~ capEstimates, fst) %>% summary
lm(ylog ~ capEstimates, fst) %>% summary

fst$fstScaled <- (scale(fst$Fst)[,1] + 1.19740)

fst[,c("capEstimates","trip")]

filter(fst, y == 0) %>% select(Population1, Population2,
                               trip, capEstimates, Fst)


seh <- emR::em.data.session 

fst$session <- paste0(fst$trip, "_", fst$Population1)
fst$sesh2 <- paste0(fst$trip, "_", fst$Population1)

fst2 <- fst %>% left_join(seh[,c("session", "phAbundance")],
                  by = c("sesh" = "session")) %>% 
  left_join(seh[,c("session", "phAbundance")],
            by = c("sesh2" = "session")) %>% 
  mutate(abundance = (phAbundance.x + phAbundance.y)/2)


fst2[, c("trip", "phAbundance.x",  "phAbundance.y")]


ggplot(fst2, aes(abundance, ylog)) +
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(~peakNo)


ggplot(fst, aes(x,y)) +
  geom_point() + geom_line()
