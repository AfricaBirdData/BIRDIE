# Pattern to make

ptt <- "here/occu_fit_spOccupancy_2011_here_ZA.rds"

spp <- sort(BIRDIE::waterbirds$SppRef)

out_ptt <- vector(length = length(spp))

for(i in seq_along(out_ptt)){
    out_ptt[i] <- gsub("here", spp[i], ptt)
}

writeLines(out_ptt, con = "path_list.txt")
