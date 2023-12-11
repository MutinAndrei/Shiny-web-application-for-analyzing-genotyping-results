# Shiny-web-application-for-analyzing-genotyping-results
This web application was created to analyze genotyping of cryptic amphipod species of the Baikal region using two marker genes 18S and COI. However, the code used in this application can be applied not only to Baikal amphipods, but also to all other species of living organisms. The application was written in R language using the shiny package. The blast code was written based on https://github.com/ScientistJake/Shiny_BLAST.
# Required packages:
library(**DT**) Version: **0.29**

library(**thematic**) Version: **0.1.3**

library(**bslib**) Version: **0.5.1**

library(**XML**) Version: **3.99-0.14**

library(**plyr**) Version: **1.8.8**

library(**dplyr**) Version: **2.3.4**

library(**rBLAST**) Version: **0.99.2**

library(**ggplot2**) Version: **3.4.3**

library(**tidyr**) Version: **1.3.0**

library(**xlsx**) Version: **0.6.5**

library(**shinyalert**) Version: **3.0.0**

library(**Biostrings**) Version: **2.68.1**

library(**shinycssloaders**) Version: **1.0.0**

library(**ggtext**) Version: **0.1.2**

library(**glue**) Version: **1.6.2**

library(**rclipboard**) Version: **0.2.0**

library(**BiocManager**) Version: **1.30.22**

library(**shinyBS**) Version: **0.61.1**

library(**shinyjs**) Version: **2.1.0**

library(**stringr**) Version: **1.5.0**

library(**seqinr**) Version: **4.2-30**

library(**forcats**) Version: **1.0.0**

library(**purrr**) Version: **1.0.2**

library(**lubridate**) Version: **1.9.2**

library(**readr**) Version: **2.1.4**

library(**tibble**) Version: **3.2.1**

library(**sf**) Version: **1.0-14**

library(**mapview**) Version: **2.11.2**

library(**leaflet**) Version: **2.2.0**

# Customizing Specoident:

You can use this application with your local databases and live detection coordinates. 

To do this, you need to load your database into the application, file app.R line 181.

Load the coordinates of the live detection in line 214.

Change the search pattern in the conditions that produce popups in lines 286 through 496.

Also, this application is located on our server:

http://bioinformatics.isu.ru:3838/speCOIdent/ 
