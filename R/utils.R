

DKI_group_samples <- function(df, group_vars=c(),
                              column_sample_id=NA_character_,
                              column_samples_per_group=NA_character_){

  if(is.na(column_sample_id)){
    stop("column_sample_id, the name of the column containing sample id must be specified!")
  }

  group_vars_sample <- c(group_vars, column_sample_id, column_samples_per_group)
  group_vars <- c(group_vars, column_samples_per_group)

  df_samples <- df %>%
    group_by(across(all_of(group_vars_sample))) %>%
    summarise(n_spec=n(), .groups="drop")

  # we don't know if existing sample numbers are consecutive within each site
  # sometimes one could be missing
  # we will introduce some that we are sure are consecutive (they will be dropped later)

  if(length(group_vars)>0){
    df_samples <- df_samples %>%
      group_by(across(all_of(group_vars)))
  }

  df_samples <- df_samples %>%
    arrange(!!column_sample_id) %>%
    mutate(n_samples=n(), xsample_id=row_number()) %>%
    ungroup()

  df_samples <- df_samples %>%
    rowwise() %>%
    mutate(grp_id=.haps_group(xsample_id, n=n_samples, maxsize=!!sym(column_samples_per_group))) %>%
    ungroup() %>%
    select(all_of(c(group_vars_sample, "grp_id")))

  df <- df %>%
    left_join(df_samples, by=group_vars_sample)

  return(df)

}

.GetSpeciesIDs<-function(searchtext, verbose=F){
  # ------ get  AphiaID from search text using fuzzy match ---------
  #Build the URL to get the data from
  
  searchtext2 <- gsub(" ","%20",searchtext)
  searchtext2 =paste0("&scientificnames%5B%5D=",searchtext2)
  searchtext2 <- paste0(searchtext2, collapse="")
  
  url<-sprintf("https://marinespecies.org/rest/AphiaRecordsByMatchNames?marine_only=true&extant_only=true&match_authority=true%s",searchtext2)  
  
  x<-http_status(GET(url))
  AphiaRecord <- data.frame(AphiaID=NULL)
  if(x$reason!="OK"){
    if(verbose==T){
      cat(paste0(searchtext,": ",x$reason,"\n"))
    }
    return(AphiaRecord)
  }
  AphiaRecord <- fromJSON(url) 

  if(length(AphiaRecord)>0){
    names(AphiaRecord) <-  searchtext
      
    AphiaRecord <-  AphiaRecord %>%
      bind_rows(.id="species") %>%
      filter(status=="accepted")
  }else{
    AphiaRecord <- data.frame(AphiaID=NULL)
    }
  return(AphiaRecord)
}


match_species <- function(df, column_species="Artsnavn", 
                          column_species_corrected="species",
                          verbose=T){
  
  # function trying matches in GBIF and WoRMS
  
  # Get unique species names from the observation data
  column_species_matched <- "species_match"
  
  species <- df %>%
    distinct(!!as.name(column_species)) 
  
  # Use the `.fix_name()` function from [utils.R](R/utils.R). 
  # This will  do some simple corrections to species names 
  # e.g. replace any *indet.* in  species names and remove 
  # any double spaces in names.
  
  species <- species %>% 
    mutate(!!column_species_corrected := 
             .fix_name(!!as.name(column_species)))  %>%
    mutate(!!column_species_corrected := 
             gsub(" sp.","",!!as.name(column_species_corrected)))
  
  
  # start with gbif - it's faster
  species_list <- species %>%
    pull(column_species_corrected)
  
  if(verbose==T){
    cli::cli_inform(paste0("Checking ", length(species_list), " names in GBIF"))
  }
  
  
  res_gbif <- rgbif::name_backbone_checklist(species_list)
  res_unmatched <- res_gbif %>%
    filter(is.na(usageKey))
  
  if(verbose==T){
    nmissing <- res_gbif %>%
      filter(is.na(usageKey))%>%
      nrow()
    nfound <- length(species_list) - nmissing
    cli::cli_inform(paste0("...matched ", nfound, ""))
  }
  
  # check worms 
  species_list <- res_unmatched$verbatim_name %>% 
    unique()
  
  nmax <- 40
  n <- seq_along(species_list)
  
  if(verbose==T){
    cli::cli_inform(paste0("Checking ", nrow(res_unmatched), " names in WoRMS"))
  }
  species_split <- split(species_list, ceiling(n/nmax))
  
  res_worms <- species_split %>%
    purrr::map(.GetSpeciesIDs, .progress = T) %>%
    bind_rows()
  
  if(verbose==T){
    nfound <- res_worms %>%
      nrow()
    cli::cli_inform(paste0("...matched ", nfound, ""))
  }
  
  res_gbif <- res_gbif %>%
    mutate(rank=tolower(rank)) %>%
    select(verbatim_index,
           !!as.name(column_species_corrected):=verbatim_name,
           !!as.name(column_species_matched):=canonicalName,
           rank,
           kingdom,
           phylum,
           family,
           genus,
           class,
           order,
           id=usageKey) %>%
    mutate(source="gbif")
  
  res_worms <- res_worms %>%
    select(!!as.name(column_species_corrected):=species, 
           !!as.name(column_species_matched):=scientificname,
           rank,
           kingdom,
           phylum,
           family,
           genus,
           class,
           order,
           id=AphiaID) %>%
    mutate(source="worms")
  
  res_worms <- res_worms %>%
    left_join(res_gbif %>% 
                select(!!as.name(column_species_corrected), verbatim_index),
              by=c(column_species_corrected))
  
  res <- res_gbif %>%
    filter(!verbatim_index %in% res_worms$verbatim_index) %>%
    bind_rows(res_worms) %>%
    arrange(verbatim_index) %>%
    select(-verbatim_index)
  
  
  res <- species %>%
    left_join(res, by=c(column_species_corrected),
              relationship = "many-to-many") %>%
    select(-!!as.name(column_species_corrected))
  
  names(res)[names(res)==column_species_matched] <- column_species_corrected 
  
  # if(verbose==T){
  #   nfound <- res %>%
  #     filter(!is.na(id))
  #     nrow()
  #   ncorrected <- res %>%
  #       filter(!is.na(id)) %>%
  #     filter(!!as.name(column_species) != !!as.name(column_species_corrected)) %>%
  #     nrow()
  #     
  #   cli::cli_inform(paste0("Returning ", nrow(res), " species"))
  #   cli::cli_inform(paste0("matched: ", found, " species"))
  #   cli::cli_inform(paste0("corrected: ", ncorrected, " species"))
  # }
  # 
  
  return(res)
}






.GetSpeciesIDFuzzy<-function(searchtext, verbose=F){
  # ------ get  AphiaID from search text using fuzzy match ---------
  #Build the URL to get the data from
  searchtext2 <- gsub(" ","%20",searchtext)
  url<-sprintf("http://marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]=%s&marine_only=false",searchtext2)
  x<-http_status(GET(url))
  ID<-""
  if(x$reason!="OK"){
    if(verbose==T){
      cat(paste0(searchtext,": ",x$reason,"\n"))
    }
    return(ID)
  }
  AphiaRecord <- fromJSON(url) 
  if(length(AphiaRecord)>0){
    AphiaRecord <- as.data.frame(AphiaRecord[[1]])
    if(nrow(filter(AphiaRecord,status=="Accepted"))>0){
      AphiaRecord <- AphiaRecord %>% 
        filter(status=="Accepted")
    }
    if(!is.na(AphiaRecord$scientificname[1])){
      ID<-AphiaRecord$AphiaID[1]
      if(verbose==T){
        cat(paste0(searchtext,": AphiaID=",ID,"\n"))
        }
    }
  }
  return(ID)
  
}


.GetSpeciesNameWoRMS<-function(searchtext, verbose=F){
  # ---------- get the AphiaID from the search text ---------------
  #Build the URL to get the data from
  
  if(grepl("\\/", searchtext)==TRUE){
    # invalid character
    return(df)
  }
  
  searchtext <- URLencode(searchtext)
  url<-sprintf("https://marinespecies.org/rest/AphiaIDByName/%s?marine_only=false",searchtext)
  #browser()
  x<-http_status(GET(url))
  if(x$reason!="OK"){
    if(verbose==T){
      cat(paste0(searchtext,": ",
                 x$reason," (trying fuzzy search instead)\n"))
    }
    
    AphiaID <- .GetSpeciesIDFuzzy(searchtext, verbose=verbose)
    if(AphiaID==""){
      return(NA_character_)
    }
  }else{
    #Get the AphiaID
    AphiaID <- fromJSON(url)
    if(verbose==T){
      cat(paste0(searchtext,": AphiaID=",AphiaID))
    }
  }
  
  # ---------- get the Aphia record from the AphiaID ---------------
  
  url<-sprintf("http://marinespecies.org/rest/AphiaRecordByAphiaID/%d",AphiaID)
  url
  AphiaRecord <- fromJSON(url)
  
  
  validID<-AphiaRecord$valid_AphiaID
  if(is.null(validID)){
    correct_name <- NA_character_
  }else{
    if(validID != AphiaID){
      if(verbose==T){
        cat(paste0(" (Using AphiaID=",validID,")\n"))
      }
      AphiaIDorig <- AphiaID
      AphiaID <- validID
      AphiaRecordOrig<-AphiaRecord
      
      # get the correct record
      url<-sprintf("http://marinespecies.org/rest/AphiaRecordByAphiaID/%d",AphiaID)
      AphiaRecord <- fromJSON(url)
    }else{
      if(verbose==T){
        cat(paste0(" (Valid ID)\n"))
      }
    }
  }
  
  correct_name <- AphiaRecord$scientificname
  return(correct_name)
}

#' functions for matching species name with names in a list
#' ignoring whether or not sp. is used in the AMBI list or the observation
#' species name

.species_match_single<- function(species, species_match, exact=T){

  species <- ifelse(is.na(species),"",species)
  res <- species_match[species_match==species]
  
  if(length(res)==0){
    if(exact!=T){
      species_match2 <- stringr::str_remove_all(species_match, " sp.")
      species2 <- stringr::str_remove_all(species, " sp\\.")
      species2 <- stringr::str_remove_all(species2, " spp\\.")
      res <- species_match[species_match2==species2]
      if(length(res)==0){
        res <- NA_character_
      }
    }else{
      res <- NA_character_
    }
  }
  if(length(res)>1){
    n_res <- length(res)
    msg <- paste0(res, collapse = ", ")
    msg <- paste0(species, " had ", n_res, " matches: ", msg)
    cli::cli_warn(msg)
    res <- res[1]
  }
  return(res)
}

.species_match <- function(df_obs, df_species,
                           var_species_obs="species",
                           var_species="species",
                           exact=T){
  
  species_list <- df_species %>%
    pull(var_species)
  
  species_match <- df_obs %>%
    pull(var_species_obs) %>%
    lapply(.species_match_single, species_list, exact) %>%
    unlist()
  
  df_obs$species_match <- species_match
  
  return(df_obs)
}

.fix_name <- function(species){
  # if the name includes a comma, extract only the part before the comma
  # e.g. "Polychaeta, havbørsteorme" -> "Polychaeta"
  species <- stringr::str_remove(species, ",.+")
  species <- stringr::str_remove(species, "\\(\\w+\\)")
  species <- stringr::str_to_sentence(species)
  species <- stringr::str_replace_all(species, "  ", " ")
  species <- stringr::str_replace_all(species, " indet.", "")
  species <- stringr::str_remove(species, "_(\\w+_)")
  
  return(species)
  }

.valid_species <- function(species){
  species <- ifelse(species=="Atylus swammerdami",
                    "Nototropis swammerdamei",species)
  species <- ifelse(species=="Callinassa subterranea",
                    "Callianassa subterranea",species)
  species <- ifelse(species=="Cheirocratus sundevalli",
                    "Cheirocratus sundevallii",species)
  species <- ifelse(species=="Corophium bonelli",
                    "Corophium bonnellii",species)
  species <- ifelse(species=="Diastyloides biplicata",
                    "Diastyloides biplicatus",species)
  
  species <- ifelse(species=="Pseudocuma longicornis",
                    "Pseudocuma longicorne",species)
  species <- ifelse(species=="Priapulidae sp.",
                    "Priapulida sp.",species)
  species <- ifelse(species=="Terebellida sp.",
                    "Terebellidae sp.",species)
  
  # here AMBI is not correct - it should be Priapulidae
  
  species <- ifelse(species=="",
                    "",species)
  species <- ifelse(species=="",
                    "",species)
  return(species)
}


# auxiliary function to group samples by sample number

.haps_group <- function(id, n, maxsize=7){
  if(n<1){
    stop("n must be greater than or equal to 1")
  }
  if(id<1){
    stop("id must be greater than or equal to 1")
  }
  if(id>n){
    stop("id must be less than or equal to n")
  }
  n_grps <- ceiling(n / maxsize)
  sizes <- rep(maxsize, n_grps)
  diff <- sum(sizes) - n
  while(diff > n_grps){
    sizes <- sizes - 1
    diff <- sum(sizes) - n
  }
  if(diff>0){
    sizes[1:diff] <- sizes[1:diff] - 1
  }
  grps <- 1:n_grps
  grps <- rep(grps, rep, times=sizes)
  grp <- grps[id]
  return(grp)
}

