
# myDF<-myInnerDataFrame
# Variable<-"Value"
# Response<-"LongitudinalVariable"
# Subsets<-"CategoricalVariable"
# IDs<-"LinkVariable"



PairedWilcox_withNAs<-function(myDF, Variable, Response, Subsets=NULL, IDs){
  require(utils)

  myInnerDF<-myDF %>%
    dplyr::rename(Response = !!sym(Response),
                Variable = !!sym(Variable),
                IDs = !!sym(IDs))

  subs <- myDF %>%
    pull(Response) %>%
    unique()

  myPossiblePairs <-
    combn(subs,2, simplify = F) %>%

      map(function(comb){
        nlev <-
          myInnerDF %>%
          dplyr::select(IDs, Variable, Response) %>%
          dplyr::filter(Response %in% comb &
                          !is.na(Variable)) %>%
          pull(Response) %>%
          unique() %>%
          length()

        if(nlev < 2) {
          return(NULL)
        } else {
          return(comb)
        }
      }) %>%
        purrr::discard(is.null)



  myPossiblePairs %>%
    purrr::set_names() %>%
    purrr::map_dfr(function(set){

      if(is.null(Subsets)){
        try(
          myTestDF <-
            myInnerDF %>%

            dplyr::filter(Response %in% set) %>%
            dplyr::mutate(Variable = as.numeric(Variable),
                          Response = as.character(Response)) %>%
            tibble::as_tibble() %>%
            dplyr::arrange(IDs, Response) %>%
            dplyr::group_by(IDs) %>%
            dplyr::filter(n()==2) %>%
            dplyr::ungroup() %>%
            rstatix::wilcox_test(.,Variable ~ Response, paired = T, p.adjust.method = "BH") %>%
            rstatix::add_y_position() %>%
            dplyr::select(-groups)
        )

      } else {
        try(
          myTestDF <-
            myInnerDF %>%

            dplyr::filter(Response %in% set,
                          !is.na(Variable)) %>%
            dplyr::mutate(Variable = as.numeric(Variable),
                   Response = as.character(Response)) %>%
            tibble::as_tibble() %>%
            dplyr::arrange(IDs, Response) %>%
            dplyr::group_by(IDs,!!sym(Subsets)) %>%
            dplyr::filter(n()==2) %>%
            dplyr::ungroup(IDs) %>%
            rstatix::wilcox_test(.,Variable ~ Response, paired = T, p.adjust.method = "none")  %>%
            rstatix::add_y_position() %>%
            dplyr::select(-groups)
        )

      }


    }) %>% mutate(.y. = paste(Variable))

  }

# test<-PairedWilcox_withNAs(myDF = myInnerDataFrame,Subsets = "CategoricalVariable" ,Variable="Value", Response="LongitudinalVariable", IDs="LinkVariable")
# test<-test

# colnames(test)[1] <- "LinkVariable"
