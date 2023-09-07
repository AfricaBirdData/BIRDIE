# This script is used to export the list of species codes used by BIRDIE to a text file

sp_codes <- unique(BIRDIE::waterbirds$SppRef)
sp_codes <- sp_codes[!is.na(sp_codes)]

# If we want one species per line run the following
writeLines(as.character(sp_codes), con = "analysis/scripts/linux/spp_codes.txt")

# If we want species codes separated by a space run this
cat(as.character(sp_codes), file = "analysis/scripts/linux/spp_codes.txt")
