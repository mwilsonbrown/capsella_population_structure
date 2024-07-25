# Functions for Capsella population definitions
# Maya Wilson Brown
# March 11, 2024

####### Colors --------
# Roma (geographic regions)
region.cols <- c("E. Asia" = "#7E1700", "MENA/Medditeranean"= "#B58B31", "N. Europe"= "#C0E9C2", "US"= "#389CC6","C. Asia"="#023198")

# Hawaii (populations) -- chosen
pop.cols <- c("E_Asia" = "#B2F2FD", "MENA"= "#996330", "N_Europe"= "#7FC55F", "NYC"= "#585380")

###### MAP BASES ---------
states <- ne_states(returnclass = "sf", country="United States of America")

east <- ne_countries(returnclass = "sf", continent = c("Europe", "Asia", "Africa"))

whole <- ne_countries(returnclass = "sf")

# map of New York and New Jersey study area
ny_map <- ggplot()+ geom_sf(data= states, fill="gray93")+ 
  coord_sf(xlim = c(-74.5, -73),  ylim = c(40.5, 41), expand = FALSE) + 
  theme_light() + 
  theme(axis.text.x= element_text(size = 5), axis.text.y= element_text(size = 5))

# map of eurasia
eurasia <- ggplot()+ geom_sf(data= east, fill="gray93")+ 
  coord_sf(xlim = c(-30, 150.9),  ylim = c(22, 75), expand = FALSE) + 
  theme_light() + 
  theme(axis.text.x= element_text(size = 5), axis.text.y= element_text(size = 5))
# map of the world
world <- ggplot()+ geom_sf(data= whole, fill="gray93")+ 
  theme_light() + 
  theme(axis.text.x= element_text(size = 5), axis.text.y= element_text(size = 5))

######## DATA ANALYSIS FUNCTIONS-------
# Select K
plotK <- function(cv_error_log){
  cv.plot <- plot(x=cv_error_log$K_groups, y =cv_error_log$V2, ylab = "cross validation error rate", xlab = "K groups")
  return(cv.plot)
}

# K is the number of groupings for which ADMIXTURE has the initial lowest cross validation error rate
selectK <- function(cv_error_log = cv){
  slopes <- (lead(cv_error_log$V2)-cv_error_log$V2)/(lead(1:nrow(cv_error_log))-c(1:nrow(cv_error_log))) # slopes
  n.slopes <- slopes[slopes < 0] #only negative slopes
  K = which.max(n.slopes) # shallowest negative slope is biggest number (slope closest to zero)
  
  # plot
  cv.plot <- plot(cv_error_log$V2, ylab = "cross validation error rate", xlab = "K groups")
  
  # save returns as list ro capture rather than print K
  #return(as.list(K, cv.plot))
  return(K) #just prints K
  return(cv.plot)
}

#### Combine ancestry information with location data

#dyplr needed for rows_patch; stringr needed for str_remove and str_split; tidygeocoder needed for geocode

ancestry_location <- function(ancestry_data, whole_genome_sequencing){
  # make new column with the trailing S[number] removed from Panko samples
  ancestry_data$sample_name <- str_remove_all(ancestry_data$vcf_sample_name, "_S[0-9]+")
  
  # combine the sample info with ancestry
  anc_info <- left_join(ancestry_data, whole_genome_sequencing)
  # select those without latitude and longitude
  no_locs <- unique(anc_info[is.na(anc_info$latitude) & is.na(anc_info$longitude), c("sample_name","citation", "note")])
  #get city and state
  no_locs$loc_name <- str_split_i(no_locs$note, pattern = ";", i=1)
  #further split since I am having trouble passing directly to geocode
  no_locs$city <- str_split_i(no_locs$loc_name, pattern = ",", i=1)
  no_locs$country <-  str_split_i(no_locs$loc_name, pattern = ",", i=2)
  
  # takes about 1 second per query
  new_locs <- no_locs %>%
    geocode(city = city,
            country = country,
            lat = "latitude",
            long = "longitude")
  
  new_locs <- new_locs[c("sample_name", "latitude", "longitude")]
  # update with new location information; matches by sample_name
  anc_info <- rows_patch(anc_info, new_locs)
  return(anc_info)
}

# wanting to add a column indicating if the location was approximated; untested
ancestry_location2 <- function(ancestry_data, whole_genome_sequencing){
  # make new column with the trailing S[number] removed from Panko samples
  ancestry_data$sample_name <- str_remove_all(ancestry_data$vcf_sample_name, "_S[0-9]+")
  
  # combine the sample info with ancestry
  anc_info <- left_join(ancestry_data, whole_genome_sequencing)
  # select those without latitude and longitude
  no_locs <- unique(anc_info[is.na(anc_info$latitude) & is.na(anc_info$longitude), c("sample_name","citation", "note")])
  #get city and state
  no_locs$loc_name <- str_split_i(no_locs$note, pattern = ";", i=1)
  #further split since I am having trouble passing directly to geocode
  no_locs$city <- str_split_i(no_locs$loc_name, pattern = ",", i=1)
  no_locs$country <-  str_split_i(no_locs$loc_name, pattern = ",", i=2)
  
  # takes about 1 second per query
  new_locs <- no_locs %>%
    geocode(city = city,
            country = country,
            lat = "latitude",
            long = "longitude")
  
  new_locs <- new_locs[c("sample_name", "latitude", "longitude")]
  # add columns for location approximation indicator
  new_locs$approx_location <- "yes"
  
  # make new column for non-interpolated data
  anc_info$approx_location <- "no"
  
  # update with new location information; matches by sample_name
  anc_info <- rows_patch(anc_info, new_locs)
  return(anc_info)
}

## Assign majority ancestry group
assign_ancestry <- function(ancestry_proportion_list, list_element_name){
  # Print selected K and set K
  K = as.numeric(str_extract(list_element_name, "[0-9]+")) # extract only the number from list name and treat as number
  print(paste0("Selected K groups is ", K))
  
  # select dataframe associated with K
  props <- ancestry_proportion_list[[list_element_name]]
  # for each row (apply,1) select the column name in which the maximum value (of ancestry proportion) is found; replace the column indicator 'V' with 'pop'
  group <- str_replace(colnames(props)[apply(props[1:K],1, which.max)], "V", "pop")
  # add sample names back
  df <- as.data.frame(cbind(props, group))
  return(df)
}

######## PLOTTING FUNCTIONS
# Ancestry Bars Plot--ordered by geographic region
ancestry_cols <- function(ancestry_proportion_list, list_element_name, ggplot_opt = NULL){
  # Print selected K and set K
  K = as.numeric(str_extract(list_element_name, "[0-9]+")) # extract only the number from list name and treat as number
  print(paste0("Selected K groups is ", K))
  
  # select dataframe associated with K
  ancestry_dat <- ancestry_proportion_list[[list_element_name]]
  # change column names to vector of various K populations and vcf names
  colnames(ancestry_dat) <- str_replace(colnames(ancestry_dat), "V", "pop")
  
  # order region as a factor with specific levels, then arrange them by region
  ancestry_dat$region <- factor(ancestry_dat$region, levels = c("US","N. Europe", "MENA/Mediterranean", "C. Asia", "E. Asia"))
  df2 <- ancestry_dat %>% arrange(region)
  
  # Set the levels of the vcf sample names based on region order
  ind.lvl <- df2$vcf_sample_name #get level order that I want (by region)
  
  # then pivot longer
  df3 <- pivot_longer(df2, 
                      cols = colnames(df2)[1:K], 
                      names_to = "ancestry", 
                      values_to = "proportion")
  # set vcf_sample_name to factor with explicit levels
  df3$vcf_sample_name <- factor(df3$vcf_sample_name, levels = ind.lvl)
  
  #create map (named vector) between vcf_sample_name and sample_name
  sn <- df3$sample_name
  names(sn) <- df3$vcf_sample_name
  
  ## Plotting
  if(!is.null(ggplot_opt)) {
    plot <- ggplot() +
      geom_col(data=df3, aes(x=vcf_sample_name, y = proportion, fill=ancestry)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
            legend.position="bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      scale_x_discrete(labels=sn) +
      ggplot_opt
  } else {
    plot <- ggplot() +
      geom_col(data=df3, aes(x=vcf_sample_name, y = proportion, fill=ancestry)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
            legend.position="bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      scale_x_discrete(labels=sn)
  }
  return(list(plot, df3))
  
}

