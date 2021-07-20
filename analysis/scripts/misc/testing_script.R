visit.data = visit_data
site.data = site_data
names = list(visit = list(site = "Pentad", obs = "detc"),
             site = list(site = "Pentad", coords = c("lon", "lat")))

so.data = so_data

site <- so.data$site
visit <- so.data$visit
site$site.idx <- factor(site[, attr(site, "site")])
visit$site.idx <- factor(visit[, attr(visit, "site")], levels = levels(site$site.idx))
xy <- as.matrix(site[, attr(site, "coords")])
X <- as.matrix(model.matrix(occupancy.model, site))
W_vb <- as.matrix(model.matrix(detection.model, visit))
Y <- matrix(visit[, attr(visit, "obs")], ncol = 1)
z <- as.matrix(ifelse(table(visit$site.idx, visit[, attr(visit,
                                                         "obs")])[, "1"] > 0, 1, 0), ncol = 1)
nvisits <- as.numeric(table(visit$site.idx))
n.obs <- length(nvisits[which(nvisits > 0)])
n.site <- nrow(X)
nvisits <- nvisits[nvisits > 0]
unsurveyed_index <- ((n.obs + 1):n.site)
siteids <- matrix(rep(1:n.obs, nvisits), ncol = 1)
ysum <- tapply(so.data$visit$PAdata, so.data$visit$SiteName,
               sum)
alpha_m <- matrix(prior$mu.d, ncol = 1)
beta_m <- matrix(prior$mu.o, ncol = 1)




# Splines -----------------------------------------------------------------

library(mgcv)

# Add also month

occudata <- occudata %>%
    mutate(month = lubridate::month(StartDate))

spl_bs <- cSplineDes(x = sort(occudata$month),
                     knots = 1:12,
                     ord = 4, derivs = 0)

head(spl_bs)

plot(sort(occudata$month), spl_bs[,1], type = "l")
for (i in 2:ncol(spl_bs)) lines(sort(occudata$month), spl_bs[,i],col=i)



# New code ----------------------------------------------------------------

pa_dat <- getOccVisitData(region_type = "province", region = "North West",
                          species = 6, .path = "analysis/data")

saveRDS(pa_dat, "analysis/out_nosync/pa_dat_6_nw.rds")

pa_dat <- readRDS("analysis/out_nosync/pa_dat_6_nw.rds")

future::plan("multisession", workers = 3)

progressr::with_progress(
    covt_dat <- addOccVisitCovts(visit_data = pa_dat$visits,
                                 sites = pa_dat$sites,
                                 covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                                 covts_dir = "analysis/out_nosync/",
                                 file_fix = c("terraClim_", "_03_19_nw"))
)

future::plan("sequential")

saveRDS(covt_dat, "analysis/out_nosync/pa_dat_6_wcovts_nw.rds")

visit_data = pa_dat$visits
sites = pa_dat$sites
covts = c("prcp", "tmax", "tmin", "aet", "pet")
covts_dir = "analysis/out_nosync/"
file_fix = c("terraClim_", "_03_19_nw")

.sub_by = yearmonth[i]
.visit_data = sub_data[[i]]
.sites = sites
.covt = covt
.name = covts[i]
.p = p

visits = pa_dat$visits
sites = pa_dat$sites
covt = "prcp"
covts_dir = "analysis/out_nosync/"
file_fix = c("terraClim_", "_03_19_nw")
