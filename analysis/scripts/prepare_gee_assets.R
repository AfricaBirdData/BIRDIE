library(rgee)
library(dplyr)
library(stars)
library(raster)

# Import the raw Python Earth Engine module
ee <- reticulate::import("ee")

# Initialize rgee with an account that is authorised for project ee-seecatuct
# When prompted make sure you select this account and select the appropriate
# project
rgee::ee_check()
rgee::ee_Initialize(user = "birdie", project = 'ee-seecatuct', drive = TRUE)

# Create the folder using the raw backend dictionary structure
# If this folder was already created we would get an error that we cannot
# overwrite
ee$data$createAsset(list(type = "FOLDER"), "projects/ee-seecatuct/assets/myassets")

# Then we need to upload the national wetland map to asset home

