

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


.GetSpeciesName<-function(searchtext, verbose=F){
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


.fix_name <- function(species){
  species <- stringr::str_to_sentence(species)
  species <- stringr::str_replace_all(species, "  ", " ")
  suffix <- ifelse(stringr::str_detect(species," "),""," sp.")
  species <- ifelse(is.na(species),NA,paste0(species, suffix))
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

