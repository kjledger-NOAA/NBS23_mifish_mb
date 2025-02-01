
root_dir = R'(C:\Users\James.Thorson\Desktop\Git\NBS23_mifish_mb)'
data_dir = file.path( root_dir, "outputs" )
#setwd( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2024-10 -- Ledger eDNA tinyVAST)' )

Version = c("Tweedie", "cloglog")[1]
Date = Sys.Date()
date_dir = file.path( root_dir, paste0(Date,"_",Version) )
  dir.create(date_dir)

library(tinyVAST)
library(fmesher)
library(sf)

# Subset to top-3
species_to_use = c("Gadus_chalcogrammus", "Oncorhynchus_keta", "Clupea_pallasii")

################
# Read and merge data sets
################

#
CSV1 = read.csv( file.path(data_dir, "nbs23_taxon_reads.csv") )
CSV1$species = gsub( x = CSV1$species, pattern = " ", replace = "_", fixed = TRUE )
CSV1$taxon = ifelse( CSV1[,'species'] %in% species_to_use, CSV1[,'species'], "Other" )
CSV1 = subset( CSV1, depth == 10 )
# Add stationID
CSV1$station_ID = paste0( CSV1$longitude, "_", CSV1$latitude )
if( FALSE ){
  # Top species
  cumsum(sort(tapply( CSV1[,'reads'], INDEX=CSV1[,'species'], FUN=sum), dec = TRUE))  / sum(CSV1[,'reads'])

  # only 2023
  table(CSV1[,'collection_year'])

  #
  plot( x=CSV1[,'longitude'], y=CSV1[,'latitude'] )
}

#
CSV2 = read.csv( file.path(data_dir, "nbs23_basis_catch.csv") )
CSV2$SCIENTIFICNAME = ifelse( CSV2$SCIENTIFICNAME == "Clupea pallasi", "Clupea_pallasii", CSV2$SCIENTIFICNAME )
CSV2$SCIENTIFICNAME = gsub( x = CSV2$SCIENTIFICNAME, pattern = " ", replace = "_", fixed = TRUE )
CSV2$taxon = ifelse( CSV2[,'SCIENTIFICNAME'] %in% species_to_use, CSV2[,'SCIENTIFICNAME'], "Other" )
if( FALSE ){
  # Top species
  cumsum(sort(tapply( CSV2[,'TOTALCATCHWT'], INDEX=CSV2[,'SCIENTIFICNAME'], FUN=sum), dec = TRUE))  / sum(CSV2[,'TOTALCATCHWT'])

  # only 2023
  table(CSV2[,'year'])

  #
  plot( x=CSV2[,'EQ_LONGITUDE'], y=CSV2[,'EQ_LATITUDE'] )
}

#
#DF = expand.grid( "sample_ID"=CSV2[,'sample_ID'], "taxon"=colnames(CSV2)[-(1:3)] )
#DF$extraction_ID = CSV2[match(DF$sample_ID,CSV2$sample_ID),'extraction_ID']
DF1 = expand.grid( "Sample"=unique(CSV1[,'station_ID']), "taxon"=unique(CSV1$taxon) )
DF1 = cbind( DF1, CSV1[match(DF1$Sample,CSV1$station_ID),c("longitude","latitude")] )
Tmp = tapply( CSV1$reads, 
              INDEX = list(
                factor(CSV1$station_ID, unique(CSV1$station_ID)), 
                factor(CSV1$taxon, unique(CSV1$taxon) )
              ), 
              FUN = sum )
Tmp = ifelse( is.na(Tmp), 0, Tmp )
DF1$response = as.vector(Tmp) / 1e3

#
DF2 = expand.grid( "Sample"=unique(CSV2[,'STATIONID']), "taxon"=unique(CSV2$taxon) )
DF2 = cbind( DF2, 
             "longitude" = CSV2[match(DF2$Sample,CSV2$STATIONID),"EQ_LONGITUDE"],
             "latitude" = CSV2[match(DF2$Sample,CSV2$STATIONID),"EQ_LATITUDE"] )
Tmp = tapply( CSV2$TOTALCATCHWT, 
              INDEX = list(
                factor(CSV2[,'STATIONID'], unique(CSV2[,'STATIONID'])), 
                factor(CSV2$taxon, unique(CSV2$taxon) )
              ), 
              FUN = sum )
Tmp = ifelse( is.na(Tmp), 0, Tmp )
DF2$response = as.vector(Tmp) / 1e3

#
DF = rbind(
  data.frame( DF1[,-1], "gear"="eDNA", "Sample"=DF1[,1] ),
  data.frame( DF2[,-1], "gear"="trawl", "Sample"=factor(DF2[,1]) )
)

################
# Run model
################

# Start recording indices
B_mcz = array( NA, dim = c(3,length(species_to_use),2),
              dimnames = list("model" = c("joint","trawl","eDNA"),
                              "species" = species_to_use,
                              "stat" = c("Est","SE")) )

# Set trawl as reference gear
DF$gear = relevel( factor(DF$gear), ref = "trawl" )

# Descriptive
tapply( DF$response, INDEX=list(DF$taxon, DF$gear), FUN=\(x){mean(x>0)} )

# thin level .. 1 = trawl or first eDNA sample
DF$thin_level = ifelse( DF$gear=="eDNA", as.numeric(DF$Sample), 1 )  
DF$thin_level = relevel( factor(DF$thin_level), ref = 1 )

# make mesh
mesh = fm_mesh_2d( DF[,c('longitude','latitude')], 
                   cutoff = 0.1 )
mesh = fm_refine(mesh)

# Formula
formula =
  # species-specific intercept
  response ~ 0 + taxon +
  # composition-data .. control for differences in total read_count across taxa for a given sample_ID
  s(thin_level, bs = "re" )
  #thin_level +
  # catchability ratio for species in eDNA relative to trawl
  gear:taxon

# Separate spatial SD
sem = "
  Gadus_chalcogrammus <-> Gadus_chalcogrammus, SD # sd_cod
  Oncorhynchus_keta <-> Oncorhynchus_keta, SD # sd_chum
  Clupea_pallasii <-> Clupea_pallasii, SD # sd_herring
  #Other <-> Other, sd_Other
"

# Drop other category for now ... 100% encounter rate is annoying to incorporate
data = droplevels(subset( DF, taxon != "Other" ))

# Add error specification
#data$error_level = data$taxon
data$error_level = interaction(data$gear, data$taxon)
if( Version == "Tweedie" ){
  family = list(
    #"Gadus_chalcogrammus" = tweedie(link="log"),
    #"Oncorhynchus_keta" = tweedie(link="log"),
    #"Clupea_pallasii" = tweedie(link="log")
    "trawl.Gadus_chalcogrammus" = tweedie(link="log"),
    "trawl.Oncorhynchus_keta" = tweedie(link="log"),
    "trawl.Clupea_pallasii" = tweedie(link="log"),
    "eDNA.Gadus_chalcogrammus" = tweedie(link="log"),  # eDNA for pollock is 100% encounter
    "eDNA.Oncorhynchus_keta" = tweedie(link="log"),
    "eDNA.Clupea_pallasii" = tweedie(link="log")#,
    #"Other" = tweedie(link="log")
  )
}
if( Version == "cloglog" ){
  family = list(
    "Gadus_chalcogrammus" = binomial(link="cloglog"),
    "Oncorhynchus_keta" = binomial(link="cloglog"),
    "Clupea_pallasii" = binomial(link="cloglog")#,
    #"Other" = tweedie(link="log")
  )
  data$response = ifelse( data$response>0, 1, 0 )
}

# fit model
out = tinyVAST( data = data,
                formula = formula,
                spatial_graph = mesh,
                sem = sem,
                family = family,
                space_columns = c("longitude", "latitude"),
                variable_column = "taxon",
                distribution_column = "error_level",
                control = tinyVASTcontrol(
                  trace = 1,
                  calculate_deviance_explained = TRUE
                ) )

# Diagnostics:  simulate new data conditional on fixed and random effects
y_ir = replicate( n = 1000,
           expr = out$obj$simulate()$y_i )

# Plot using DHARMa
res = DHARMa::createDHARMa( simulatedResponse = y_ir,
                            observedResponse = data$response,
                            fittedPredictedResponse = fitted(out) )
plot(res)
plot(res, form = data$taxon)
plot(res, form = interaction(data$taxon, data$gear))

# Extract domain
locs = unique( data[,c("longitude", "latitude")] )
locs_sf = st_as_sf( locs, coords=c("longitude", "latitude") )
domain_sf = st_convex_hull(st_union(locs_sf))
pred_sf = st_make_grid( domain_sf, n=c(50,50) )
pred_sf = st_intersection(pred_sf,domain_sf)

# Record predicted densities
predDF = data.frame( st_coordinates(st_centroid(pred_sf)), "area" = st_area(pred_sf) )
colnames(predDF) = c("longitude","latitude", "area")
# Fill in dummie values .. should replace with something better (!)
predDF$gear = "trawl"
predDF$thin_level = 1

# Loop through species clumsily
# For Tweedie, use g^{-1}(p) as proportional to read-count density
if( Version == "Tweedie" ){
  predDF$taxon = "Gadus_chalcogrammus"
  predDF$error_level = paste("trawl.Gadus_chalcogrammus")
  predDF$log_pollock = predict(out, newdata=predDF, what="mu_g" )
  B_mcz["joint","Gadus_chalcogrammus",c("Est","SE")] = integrate_output( out, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]

  predDF$taxon = "Oncorhynchus_keta"
  predDF$error_level = paste("trawl.Oncorhynchus_keta")
  predDF$log_chum = predict(out, newdata=predDF, what="mu_g" )
  B_mcz["joint","Oncorhynchus_keta",c("Est","SE")] = integrate_output( out, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]

  predDF$taxon = "Clupea_pallasii"
  predDF$error_level = paste("trawl.Gadus_chalcogrammus")
  predDF$log_herring = predict(out, newdata=predDF, what="mu_g")
  B_mcz["joint","Clupea_pallasii",c("Est","SE")] = integrate_output( out, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]
}

# plot on map
#pred_sf = st_as_sf( predDF[,c("longitude","latitude","log_cod","log_chum","log_herring")], coords=c("longitude","latitude") )
plot_sf = st_sf( pred_sf, predDF[,c("log_pollock","log_chum","log_herring")] )
png( file=file.path(date_dir,"Density.png"), width=4, height=10, unit="in", res=200 )
  plot(plot_sf, pch=20, cex=2, border=NA )   #
dev.off()

################
# Compare with only trawl
################

data_trawl = droplevels(subset( data, gear=="trawl"))
data_trawl$error_level = data_trawl$taxon
if( Version == "Tweedie" ){
  family = list(
    "Gadus_chalcogrammus" = tweedie(link="log"),
    "Oncorhynchus_keta" = tweedie(link="log"),
    "Clupea_pallasii" = tweedie(link="log")
    #"Other" = tweedie(link="log")
  )
}

# Formula
formula =
  # species-specific intercept
  response ~ 0 + taxon
 
# fit model
out_trawl = tinyVAST( data = data_trawl,
                formula = formula,
                spatial_graph = mesh,
                sem = sem,
                family = family,
                space_columns = c("longitude", "latitude"),
                variable_column = "taxon",
                distribution_column = "error_level",
                control = tinyVASTcontrol(
                  trace = 1,
                  calculate_deviance_explained = FALSE
                ) )

# Loop through species clumsily
# For Tweedie, use g^{-1}(p) as proportional to read-count density
if( Version == "Tweedie" ){
  predDF$taxon = "Gadus_chalcogrammus"
  predDF$error_level = "Gadus_chalcogrammus"
  predDF$log_pollock = predict(out_trawl, newdata=predDF, what="mu_g" )
  B_mcz["trawl","Gadus_chalcogrammus",c("Est","SE")] = integrate_output( out_trawl, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]

  predDF$taxon = "Oncorhynchus_keta"
  predDF$error_level = paste("Oncorhynchus_keta")
  predDF$log_chum = predict(out_trawl, newdata=predDF, what="mu_g" )
  B_mcz["trawl","Oncorhynchus_keta",c("Est","SE")] = integrate_output( out_trawl, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]

  predDF$taxon = "Clupea_pallasii"
  predDF$error_level = paste("Gadus_chalcogrammus")
  predDF$log_herring = predict(out_trawl, newdata=predDF, what="mu_g")
  B_mcz["trawl","Clupea_pallasii",c("Est","SE")] = integrate_output( out_trawl, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]
}


# plot on map
plot_sf = st_sf( pred_sf, predDF[,c("log_pollock","log_chum","log_herring")] )
png( file=file.path(date_dir,"Density_onlytrawl.png"), width=4, height=10, unit="in", res=200 )
  plot(plot_sf, pch=20, cex=2, border=NA )   #
dev.off()

################
# Compare with only eDNA
################

data_eDNA = droplevels(subset( data, gear=="eDNA"))
data_eDNA$error_level = data_eDNA$taxon

# Formula
formula =
  # species-specific intercept
  response ~ 0 + taxon + s(thin_level, bs = "re")


# fit model
out_eDNA = tinyVAST( data = data_eDNA,
                formula = formula,
                spatial_graph = mesh,
                sem = sem,
                family = family,
                space_columns = c("longitude", "latitude"),
                variable_column = "taxon",
                distribution_column = "error_level",
                control = tinyVASTcontrol(
                  trace = 1,
                  calculate_deviance_explained = FALSE
                ) )
 
# Loop through species clumsily
# For Tweedie, use g^{-1}(p) as proportional to read-count density
if( Version == "Tweedie" ){
  predDF$taxon = "Gadus_chalcogrammus"
  predDF$error_level = paste("Gadus_chalcogrammus")
  predDF$log_pollock = predict(out_eDNA, newdata=predDF, what="mu_g" )
  B_mcz["eDNA","Gadus_chalcogrammus",c("Est","SE")] = integrate_output( out_eDNA, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]

  predDF$taxon = "Oncorhynchus_keta"
  predDF$error_level = paste("Oncorhynchus_keta")
  predDF$log_chum = predict(out_eDNA, newdata=predDF, what="mu_g" )
  B_mcz["eDNA","Oncorhynchus_keta",c("Est","SE")] = integrate_output( out_eDNA, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]

  predDF$taxon = "Clupea_pallasii"
  predDF$error_level = paste("Gadus_chalcogrammus")
  predDF$log_herring = predict(out_eDNA, newdata=predDF, what="mu_g")
  B_mcz["eDNA","Clupea_pallasii",c("Est","SE")] = integrate_output( out_eDNA, predDF, area = predDF$area )[c('Est. (bias.correct)','Std. Error')]
}

# plot on map
plot_sf = st_sf( pred_sf, predDF[,c("log_pollock","log_chum","log_herring")] )
png( file=file.path(date_dir,"Density_onlyeDNA.png"), width=4, height=10, unit="in", res=200 )
  plot(plot_sf, pch=20, cex=2, border=NA )   #
dev.off()

##########
# Plot proportions across methods
##########

Y = array(rnorm(40), dim = c(10,2,2))
P = sweep( Y, STAT = rowSums(Y[,,1]), MARGIN=1:2, FUN = "/")

P_mcz = sweep( B_mcz, STAT = rowSums(B_mcz[,,'Est']), MARGIN = 1:2, FUN = "/")

DF = expand.grid( dimnames(P_mcz[,,'Est']) )
DF$Est = as.vector(P_mcz[,,'Est'])
DF$upper = DF$Est + 1.96 * as.vector(P_mcz[,,'SE'])
DF$lower = DF$Est - 1.96 * as.vector(P_mcz[,,'SE'])

library(ggplot2)
ggplot(data=DF, aes(x=interaction(model), y=Est, color=model)) +
  geom_point( position=position_dodge(0.9) ) +
  geom_errorbar( aes(ymax=as.numeric(upper),ymin=as.numeric(lower)),
                 width=0.25, position=position_dodge(0.9)) +
  facet_grid( rows=vars(species) )    # , scales="free"
ggsave( file.path(date_dir,"Proportions.png"), width=4, height=4 )
