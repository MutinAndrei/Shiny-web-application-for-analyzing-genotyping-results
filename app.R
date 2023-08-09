# needed packages
library(DT)
library(thematic)
library(bslib)
library(XML)
library(plyr)
library(dplyr)
library(rBLAST)
library(ggplot2)
library(tidyr)
library(xlsx)
library(shinyalert)
library(Biostrings)
library(shinycssloaders)
library(ggtext)
library(glue)
library(rclipboard)
library(BiocManager)
library(shinyBS)
library(shinyjs)
library(stringr)
library(seqinr)
library(forcats)
library(purrr)
library(lubridate)
library(readr)
library(tibble)
library(sf)
library(mapview)
library(leaflet)


custom_db <- c("COI")

options(repos = BiocManager::repositories())


appCSS <- "
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"

jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'


#The Blast code was taken from "https://github.com/ScientistJake/Shiny_BLAST"

ui <- fluidPage(useShinyjs(),
                inlineCSS(appCSS),
                
                # Loading message
                div(
                  id = "loading-content",
                  h1("Loading...")
                ),
                
                # The main app code goes here
                hidden(
                  div(
                    id = "app-content",
                  )
                ),
  theme = bs_theme(bg = "#0b3d91", fg = "white", primary = "#FCC780",
                                 base_font = font_google("Space Mono"),
                                 code_font = font_google("Space Mono")),
                tagList(
                  tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js"),
                    tags$style(type='text/css', "#button { vertical-align: middle; height: 70px; width: 120%; font-size: 50px;}"),
                    tags$script(HTML(jscode))
                  )
                ),
                
  
  shinyjs::useShinyjs(),
  id = "side-panel",
                #This block gives us all the inputs:
                
                  titlePanel(h1("speCOIdent")),
  sidebarLayout(sidebarPanel(possition = "right",
                             br(),
                             useShinyalert(force = TRUE),  # Set up shinyalert
                             actionButton("preview", "Instructions for using this app"),
                             br(),
                             br(),
                             helpText("Input your nucleotide sequence"),
  tagAppendAttributes(textAreaInput('query', 'Input sequence:', value = "", placeholder = "", width = "800px", height="300px"),
                      `data-proxy-click` = "blast"),
                             actionButton("reset_input", "Clear 'Input sequence'"),
                             br(),
                             br(),
                             helpText("Choose parametrs for identify"),
                             br(),
                             div(style="display:inline-block",
                                 selectInput("eval", "e-value:", choices=c(0.001,1e-4,1e-5,1e-10,1e-200), width="150px")),
                             div(style="display:inline-block",
                                 selectInput("id", "% ID:", choices=c(95, 97, 98, 99, 99.9, 100), width="150px")),
                             br(),
                             helpText("Click on the button 'IDENTIFY!' and wait for the result"),
                             br(),
                             br(),
                             actionButton("blast", "IDENTIFY!"),
                             downloadButton(outputId = "download_filtered",
                                             label = "Download Data Base")),
  
                              mainPanel(width = 8,
                                        h5("Species distribution map:"),
                                        leafletOutput('map', height = 775),
                                        h5("Determination results:"),
                                        rclipboardSetup(), 
                               DT::dataTableOutput("blastResults"),
                               p("", tableOutput("clicked") ),
                               verbatimTextOutput("alignment")),
                                 ),
  
                 )
  



thematic::thematic_shiny()

#The Blast code was taken from "https://github.com/ScientistJake/Shiny_BLAST"
server <- function(input, output, session){
  
  observeEvent(input$reset_input, {
    shinyjs::reset("side-panel")
  })
  
  observeEvent(input$uploadFilesBtn, {
    # When the button is clicked, wrap the code in a call to `withBusyIndicatorServer()`
    withBusyIndicatorServer("blast", {
      Sys.sleep(1)})})
  
  observeEvent(input$blast, {
    cat(input$text, "\n")
  })
  
  session$onSessionEnded(stopApp)
  
  # Simulate work being done for 1 second
  Sys.sleep(1)
  
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
  
 
  
  #page 1
  
  
  custom_db <- read.fasta("database_specoident999.fa")
  custom_db_path <- c("database_specoident99.blastdb")
  
  
 
  
  output$download_filtered <- downloadHandler(
    filename = "Data base SpeCOIdenT.fa",
    content = function(file){
      write.fasta(custom_db, names=names(custom_db), file)
    }
  )
  
  blastresults <- eventReactive(input$blast, {
    
    #gather input and set up temp file
    query <- input$query
    tmp <- tempfile(fileext = ".fa")

    
    #this makes sure the fasta is formatted properly
    if (startsWith(query, ">")){
      writeLines(query, tmp)
    } else {
      writeLines(paste0(">Query\n",query), tmp)
    }
    
    #calls the blast
    data <- system(paste0(" blastn "," -query ",tmp," -db ",custom_db_path," -evalue ",input$eval, " -perc_identity ", input$id,
                          " -outfmt 5 -max_hsps 1 -max_target_seqs 100 "), intern = T)
    xmlParse(data)
  }, ignoreNULL= T)
  
  
  coord <- read.csv("table_lang_and_latit.csv")
  #Now to parse the results...
  parsedresults <- reactive({
    if (is.null(blastresults())){}
    else {
      xmltop = xmlRoot(blastresults())
      
      #the first chunk is for multi-fastas
      results <- xpathApply(blastresults(), '//Iteration',function(row){
        query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
        Species <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        Hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
        bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
        Hsp_identity <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_identity') %>% sapply(., xmlValue)
        Hsp_align <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_align-len') %>% sapply(., xmlValue)
        eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
        Percent_identity <- signif((as.numeric(Hsp_identity) / as.numeric(Hsp_align)) * 100,digits=4)
        cbind(Species,Percent_identity, Hit_length)
      })
      #this ensures that NAs get added for no hits
      results <- rbind.fill(lapply(results,function(y){as.data.frame((y),stringsAsFactors=FALSE)}))
    }
  })
  
  str(parsedresults)
  parsedresults2 <- reactive({merge(as.data.frame(parsedresults()), as.data.frame(coord))})
  
  
    
    
    output$map <- renderLeaflet({if (is.null(blastresults())){
    } else {
      mapview(as.data.frame(parsedresults2()), xcol = "longitude", ycol = "latitude", crs = 4269, grid = FALSE)@map
    }})
    
  #makes the datatable
  output$blastResults <- DT::renderDataTable({
    datatable(if (is.null(blastresults())){
    } else {
      parsedresults() 
    }, options = list(pageLength = 5,
                      dom = 't'), rownames= FALSE)
    
  }, selection="single")
  
  parsedresults1 <- reactive({
    if (is.null(blastresults())){}
    else {
      xmltop = xmlRoot(blastresults())
      
      #the first chunk is for multi-fastas
      results <- xpathApply(blastresults(), '//Iteration',function(row){
        query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
        query_length <- getNodeSet(row, 'Iteration_query-len') %>% sapply(., xmlValue)
        hit_IDs <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
        bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
        Hsp_identity <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_identity') %>% sapply(., xmlValue)
        Hsp_align <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_align-len') %>% sapply(., xmlValue)
        eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
        perc_id <- signif((as.numeric(Hsp_identity) / as.numeric(Hsp_align)) * 100,digits=4 )
        gaps <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_gaps') %>% sapply(., xmlValue)
        cbind(query_ID,hit_IDs,hit_length,bitscore,Hsp_identity,Hsp_align,eval,perc_id, query_length, gaps)
      })
      #this ensures that NAs get added for no hits
      results <-  rbind.fill(lapply(results,function(y){as.data.frame((y),stringsAsFactors=FALSE)}))
    }
  })
  


  
  observeEvent(input$blast, {
     if (str_detect(parsedresults()[1, "Species"], ".eastern._18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (eastern)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".west._18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (west)"),
        easyClose = TRUE,
        footer = NULL,
        size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".southern._18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (southern)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".angara._18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (angara)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "vittatus_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus vittatus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "cyaneus_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus cyaneus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "cruentus_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus cruentus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "inconspicuous_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus inconspicuous</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "viridulus_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus viridulus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "maackii_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus maackii</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "testaceus_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus testaceus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "viridis_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus viridis</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "tchernykhi_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus tchernykhi</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "obtusatus_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus obtusatus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "etingovae_18S") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the 18S (small subunit of ribosomal DNA) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus etingovae</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".angara._COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (angara)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".eastern._COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (eastern)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".west._COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (west)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".southern._COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (southern)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_cyaneus_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus cyaneus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_vittatus_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus vitatus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_viridulus_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus viridulus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_inconspicuous_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus inconspicuous</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_viridis_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus viridis</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_cruentus_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("This sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus cruentus</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_maacki_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("his sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus maacki</i></b>"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], ".angara2._COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("his sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus verrucosus</i></b> (angara2)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else if (str_detect(parsedresults()[1, "Species"], "_vittatus2_COI") == TRUE) {
      showModal(modalDialog(
        title = "This sequence is the COI (cytochrome oxidase subunit I) gene",
        HTML("his sequence, with more than 95% identity, belongs to the species: <b><i>Eulimnogammarus vittatus</i></b> (2)"),
        easyClose = TRUE,
        footer = NULL, size = "xl"
      ))
    } else {
      showModal(modalDialog(
        title = "No matches found!",
        HTML("If your % ID is greater than 95 try lowering it"),
        easyClose = TRUE,
        footer = NULL, size = "m"
      ))
  }
})


  
  
  #this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmltop = xmlRoot(blastresults())
      clicked = input$blastResults_rows_selected
      tableout <- data.frame(parsedresults1()[clicked,])
      
      tableout <- t(tableout)
      names(tableout) <- c("")
      rownames(tableout) <- c("Query ID","Hit ID", "Hit length", "Bit Score", "Hsp identity", "Hsp align", 
                              "e-value", "% ID", "Query length", "Gaps")
      colnames(tableout) <- NULL
      data.frame(tableout)
    }
  },rownames =T,colnames =F)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmltop = xmlRoot(blastresults())
      
      clicked = input$blastResults_rows_selected
      
      #loop over the xml to get the alignments
      align <- xpathApply(blastresults(), '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top,mid,bottom)
      })
      
      #split the alignments every 40 carachters to get a "wrapped look"
      alignx <- do.call("cbind", align)
      splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3,clicked]),",")
      
      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind(paste0("Q-",splits[[1]][i],"\n"),paste0("M-",splits[[2]][i],"\n"),paste0("H-",splits[[3]][i],"\n"))
      })
      unlist(split_out)
    }
  })
  47
  #Data for plot and table
  

  
  #Loads the expression data and converts the data frame
  
  

  #make user guide 
  
  observeEvent(input$preview, {
    # Show a modal when the button is pressed
    showModal(modalDialog(HTML("<h3><P align='center'>This application was created to identify different morphologically similar species based on comparison of the nucleotide sequences of the COI gene<br><br><P align='center'></h3>",
                          "<span style='font-size: 22px'><P align='justify'>1. First, you need to copy the sequence of interest to the 'Input sequence' window (the sequence can be either in .fasta format or just a set of nucleotides).<br><br>
2. Click the button 'IDENTIFY' button.<br><br>
3. After a short delay, the results of the definition will appear. 
You will be presented with a dialog box containing information about the animal species definition, 
and if there is no match, the program will inform you about it. To close the dialog box, 
press 'Esc' or click on the space around the dialog box with your mouse. 
You will then be presented with a table with the top five search results.
By clicking on them you will be able to see the expanded search results 
(by default you will be presented with the species to which your COI sequence, 
percent sequence identity and length).<br><br><P align='justify'><span style='font-size: 22px'>
<P align='center'>For all questions and suggestions please contact us by e-mail:<br><br>
<span style='font-size: 30px'><b>andreimutin97@gmail.com</b><P align='center'></span>
"), easyClose = TRUE, footer = NULL, size = "xl"))
  })


  
}
      
      
      
shinyApp(ui = ui, server = server)