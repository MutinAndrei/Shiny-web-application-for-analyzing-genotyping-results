# Shiny-web-application-for-analyzing-genotyping-results
This web application was created to analyze genotyping of cryptic amphipod species of the Baikal region using two marker genes 18S and COI. However, the code used in this application can be applied not only to Baikal amphipods, but also to all other species of living organisms. The application was written in R language using the shiny package. The blast code was written based on https://github.com/ScientistJake/Shiny_BLAST.
# Required packages:
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

# Customizing Specoident:

You can use this application with your local databases and live detection coordinates. 

To do this, you need to load your database into the application, file app.R line 181.

Load the coordinates of the live detection in line 214.

Change the search pattern in the conditions that produce popups in lines 286 through 496.

Also, this application is located on our server:

http://bioinformatics.isu.ru:3838/SPecoident/ 
