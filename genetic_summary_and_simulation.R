library(shiny)
library(tidyverse)
library(plotly)
library(shinyjs)
library(emR)


fst <- em.data.ph_fst_grids
fis <- em.data.ph_fis_grids
fit <- em.data.ph_fit_grids

simD <- em.data.simBoombust_dispersal
simC <- em.data.simBoombust_constant

# UI -------------
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(width = 3,
      useShinyjs(), 
      titlePanel("Input Options:"),
      selectInput(inputId = "stat", label = "Genetic summary statistic", 
                  choices = c("Fst", "He", "Ho", "fis", "fit")),
      
      
      sliderInput(inputId = "nSize", label = "min sample size",
                  min = 1, max = 10,value = 4),
      checkboxInput(inputId = "even", label = "even sampling", value = FALSE),
      checkboxInput(inputId = "xlog", label = "log sincePeak", value = TRUE),
      checkboxInput(inputId =  "facet", label = "seperate by peak", value = F),
      checkboxInput(inputId = "combine", label = "combine peaks", F),
      
      
      (radioButtons(inputId = "dispersal", label = "dispersal rate",inline = T,
                          choices = c("0.02" = 1, "0.05" = 2, "0.10" =3,
                                      "0.15" = 4,"0.20" = 5,"0.25" = 6,"0.40" = 7),
                    selected = 6)),
      (checkboxInput(inputId = "constant", 
                           label = "constant dispersal", value = F)),
   
     
    ),
    # main panel --------------------
    mainPanel(
      
      tabsetPanel(id="inTabset", type = "tabs",
                  tabPanel("Plot", plotlyOutput("statplot")),
                  tabPanel("Drop", plotlyOutput("dropplot")),
                  tabPanel("Sim", plotlyOutput("simplot"))),
      verbatimTextOutput("model")
      
      
    )
  )
)

# server ------------------------------
server <- function(input, output, session) {
  
  # observe -------------
  
  #output$text <- renderText(input$dispersal)
  observeEvent(input$inTabset,{
    if (input$inTabset == "Plot") {
      show("stat")
      show("nSize")
      show("even")
      show("xlog")
      show("facet")
      hide("combine")
      hide("dispersal")
      hide("constant")
    }
    if (input$inTabset == "Drop") {
      show("stat")
      show("nSize")
      show("even")
      hide("xlog")
      hide("facet")
      show("combine")
      hide("dispersal")
      hide("constant")
    }
    if (input$inTabset == "Sim") {
      show("stat")
      hide("nSize")
      hide("even")
      hide("xlog")
      hide("facet")
      hide("combine")
      show("dispersal")
      show("constant")
    }
    
    
  })
  
  
  geneticData <- reactive({
    
    if (input$stat == "Fst") gData <- fst
    if (input$stat == "fit") gData <- fit
    if (input$stat %in% c("He", "Ho", "fis")) gData <- fis
    
    sampleSize <- input$nSize
    
    if (input$even) sampleSize <- sampleSize + 10 
    
    
    df <- gData[[sampleSize]]
    df$yvalue <- df[, input$stat]
    
    if (input$xlog) df$xvalue <- df$xlog else df$xvalue <- df$x 
    
    df
    
    
  })
  
  simData <- reactive({
    
    if(input$constant) simdata <- simC else simdata <- simD
    simdata <- simdata[[as.numeric(input$dispersal)]]
    
    simdata
    
  })
  
  # plot f stats ----------------------------------
  output$statplot <- renderPlotly({
    
    pdata <- geneticData()  
    xlabel <- ifelse(input$xlog, "months since peak [log]", "months since peak")
    
    ggplot(pdata, aes(xvalue, yvalue)) +
      geom_point() +
      geom_smooth(method = "lm")+
      xlab(xlabel)+
      ylab(input$stat)+
      theme_bw()-> p
    
    if (input$facet) p <- p + facet_grid(~peakNo)
    ggplotly(p) 
    
  })
  
  # plot sim -------------------------------------
  output$simplot <- renderPlotly({
    
    sdata <- simData()
    sdata$y <- sdata[,tolower(input$stat)]
    ggplot(sdata, aes(gen, y)) +
      geom_point(aes(colour = dispersal, size = popIncrease)) +
      geom_line()+
      xlab("generation")+
      ylab(input$stat)+
      theme_classic()-> p
    
    ggplotly(p) 
    #plot(as.numeric(sdata$gen), sdata$fis)
    
  })
  
  # drop plot -----------------------------------
  dropData <- reactive({
    gdata <- geneticData()
    
    gdatadiff2 <- gdata %>% filter(peakNo == "peak3" | peakNo == "peak4",
                                   sincePeak == 3 | sincePeak == 60) %>% 
      select(peakNo, sincePeak, yvalue) %>% 
      mutate(drop = "drop2", mixing = ifelse(sincePeak < 5, TRUE, FALSE))
    
    gdatadiff1 <- gdata %>% filter( peakNo == "peak2" & sincePeak > 15 |
                                      peakNo == "peak3" & sincePeak == 2) %>% 
      select(peakNo, sincePeak, yvalue) %>% 
      mutate(drop = "drop1", mixing = ifelse(sincePeak < 5, TRUE, FALSE))
    
    gdatadiff <- rbind(gdatadiff1, gdatadiff2)
    gdatadiff$mixEvent <- paste0(gdatadiff$drop, "_", gdatadiff$mixing)
    
    if (input$combine){
      gdatadiff$mixEvent <- gdatadiff$mixing
      gdatadiff$drop <- gdatadiff$mixing
    }
    
    gdatadiff
    
  })
  
  
  output$dropplot <- renderPlotly({
   
    
    ddata <- dropData()
    
    ggplot(ddata, aes(mixEvent, yvalue, colour = drop))+
      geom_boxplot()+
      theme_classic()+
      theme(legend.position = "none")+
      ylab(input$stat)-> pp
    ggplotly(pp)
    
  })
  
  # model -------------------
  output$model <- renderText({
    gdata <- geneticData()
    gdata <- filter(gdata, peakNo != "peak1") %>% 
             mutate(peakNo = factor(peakNo))
    sdata <- simData()
    ddata <- dropData()
    if (input$inTabset == "Plot") {
      
          if(input$facet){
            mlinear <- lm(yvalue ~ xvalue + peakNo, gdata)
            sm <- summary(mlinear)
            lmnames <- rownames(sm$coefficients)[2:4]
            slope <- round(sm$coefficients[2:4, 1], 6)
            pvalue <- sm$coefficients[2:4, 4]
            signif <- ifelse(sm$coefficients[2:4,4] < 0.05, "significant", "non-significant")
            Rsq <- round(sm$r.squared, 2)
            
          mm  <- c(paste(lmnames, ": slope", slope, "with a", signif, "pvalue of",
                  pvalue),paste("and an Rsquared of", Rsq))
          
          #paste(mm[1], mm[2], mm[3], mm[4], mm[5], sep = "\n")
          paste(mm, collapse = "\n")
          }else{
          mlinear <- lm(yvalue ~ xvalue, gdata)
          sm <- summary(mlinear)
          slope <- round(sm$coefficients[2, 1], 6)
          pvalue <- sm$coefficients[2, 4]
          signif <- ifelse(sm$coefficients[2,4] < 0.05, "significant", "non-significant")
          Rsq <- round(sm$r.squared, 2)
          
          paste("linear model: slope", slope, "with a", signif, "pvalue of",
              pvalue,"and an Rsquared of", Rsq)
          }
      }else if (input$inTabset == "Drop"){
      
      if (input$combine) {
        m <- t.test(yvalue ~ mixEvent, ddata)
        tstat <- m$statistic
        dfree <- m$parameter
        pvalue <- m$p.value
       paste("t test: t ", tstat, "with", dfree, "df", "and a pavalue of", pvalue)
      }else{
        
        m <- t.test(yvalue ~ mixEvent, ddata[ddata$drop == "drop1",])
        tstat <- m$statistic
        dfree <- m$parameter
        pvalue <- m$p.value
        p1 <- paste("t test mix event 1: t ", tstat, "with", dfree, "df", "and a pvalue of", pvalue)
        m <- t.test(yvalue ~ mixEvent, ddata[ddata$drop == "drop2",])
        tstat <- m$statistic
        dfree <- m$parameter
        pvalue <- m$p.value
        p2 <- paste("t test mix event 2: t ", tstat, "with", dfree, "df", "and a pvalue of", pvalue)
        paste(p1, p2, sep = "\n")
        }
      } else if(input$inTabset == "Sim"){
       sdata
        paste("The simulation has", sdata$pops[1], 
              "populations, with", sdata$popSize[1],
              "individuals in each population, and", 
              sdata$nSnps[1], "snps.")
    }
    
    
  })
}
# run app ---------------------------
shinyApp(ui = ui, server = server)

