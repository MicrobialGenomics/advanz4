get_group_comparisons <- function(dat,
                                  link_var,
                                  long_var,
                                  comps,
                                  cat_vector,
                                  num_vector,
                                  type = "categorical"){

  stopifnot("argument 'type' must have one of the following values:
            'categorical' or 'longitudinal"=type %in% c("categorical", "longitudinal"))


  if(type == "categorical"){
    cat_vector %>%
      purrr::set_names() %>%
      purrr::map(function(cv){

        num_vector %>%
          purrr::set_names() %>%
          purrr::map(function(nv){
            df <-
              dat %>%
              dplyr::select(cat_var = !!sym(cv),
                            num_var = !!sym(nv),
                            long_var = !!sym(long_var),
                            link_var = !!sym(link_var)) %>%
              dplyr::filter(cat_var %in% comps) %>%
              dplyr::filter(!is.na(num_var)) %>%
              dplyr::group_by(long_var)

            wilc <- df %>%
              wilcox_test(formula = num_var ~ cat_var, conf.level = .95, exact = T,detailed = T)

            stats <- df %>%
              dplyr::group_by(cat_var,long_var, .drop = T) %>%
              summarise(median = median(num_var, na.rm = T),
                        upper_ci = quantile(num_var, .75, type = 8, na.rm = T),
                        lower_ci = quantile(num_var,.25, type = 8, na.rm = T),
                        iqr = IQR(num_var, type = 8, na.rm = T))
            return(list(stats = stats, test = wilc))
          })
      })
  }else{
    cat_vector %>%
      purrr::set_names() %>%
      purrr::map(function(cv){
        num_vector %>%
          purrr::set_names() %>%
          purrr::map(function(nv){
            df <-
              dat %>%
              dplyr::select(cat_var = !!sym(cv),
                            num_var = !!sym(nv),
                            long_var = !!sym(long_var),
                            link_var = !!sym(link_var)) %>%
              dplyr::group_by(link_var) %>%
              dplyr::filter(long_var %in% comps)


            stats <- df %>%
              dplyr::summarise(change = num_var[long_var == comps[2]] - num_var[long_var == comps[1]]) %>%
              dplyr::right_join(df, by="link_var") %>%
              dplyr::group_by(cat_var) %>%
              dplyr::filter(!is.na(change)) %>%
              # dplyr::filter(!is.na(change)) %>%
              dplyr::select(link_var, long_var, cat_var,num_var ,change) %>%
              # dplyr::select(link_var, change, num_var) %>%
              unique() %>%
              summarise(medianChange=median(change, na.rm = T),
                        median1 = median(num_var[long_var == comps[1]], na.rm = T),
                        iqr1 = IQR(num_var[long_var == comps[1]], type = 8,na.rm = T),
                        median2 = median(num_var[long_var == comps[2]], na.rm = T),
                        iqr2 = IQR(num_var[long_var == comps[2]], type = 8,na.rm = T),
                        iqr_change = IQR(change, na.rm = T, type = 8))

            wilc <- df %>%
              dplyr::filter(n()==2) %>%
              dplyr::arrange(link_var, long_var) %>%
              dplyr::group_by(cat_var, .drop=T) %>%
              rstatix::wilcox_test(formula = num_var~long_var, paired = T, detailed = T)
            return(list(stats = stats, test = wilc))
          })

      })
  }

}
