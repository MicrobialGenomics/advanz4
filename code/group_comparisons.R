PairedWilcox_withNAs<-function(myDF, Variable, Response, Subsets=NULL, IDs){

  # This function gets the matched paired samples in a longitudinal study and removes those groups that don't have all the members.
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


    })

}



get_group_comparisons <- function(dat,
                                  link_var,
                                  long_var,
                                  comps,
                                  cat_vector,
                                  num_vector,
                                  graph_coords = F,
                                  type = "categorical"){

  # This is a function to obtain the stats of a certain conmparison, be it categorical or longitudinal. It returns pvals and statistical info (mean+IQR).

  stopifnot("argument 'type' must have one of the following values:
            'categorical' or 'longitudinal'"=type %in% c("categorical", "longitudinal"))


  if(type == "categorical"){ # Doesn't need matching pairs, groups can have different Ns.
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
            if(graph_coords){
              wilc <- df %>%
                rstatix::wilcox_test(formula = num_var ~ cat_var, conf.level = .95, exact = T,detailed = T) %>%
                rstatix::add_xy_position(x = "long_var")
            } else {
              wilc <- df %>%
                wilcox_test(formula = num_var ~ cat_var, conf.level = .95, exact = T,detailed = T)
            }

            stats <- df %>%
              dplyr::group_by(cat_var,long_var, .drop = T) %>%
              summarise(median = median(num_var, na.rm = T),
                        upper_ci = quantile(num_var, .75, type = 8, na.rm = T),
                        lower_ci = quantile(num_var,.25, type = 8, na.rm = T),
                        iqr = IQR(num_var, type = 8, na.rm = T))


            return(list(stats = stats, test = wilc))
          })
      })
  }else{ # Performs longitudinal test, keeps only those pairs of samples where both members are present
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

            if(graph_coords){
              wilc <- df %>%
                dplyr::filter(n()==2) %>%
                dplyr::arrange(link_var, long_var) %>%
                dplyr::group_by(cat_var, .drop=T) %>%
                rstatix::wilcox_test(formula = num_var~long_var, paired = T, detailed = T) %>%
                rstatix::add_xy_position(x = "long_var", group = "cat_var")

            } else {
              wilc <- df %>%
                dplyr::filter(n()==2) %>%
                dplyr::arrange(link_var, long_var) %>%
                dplyr::group_by(cat_var, .drop=T) %>%
                rstatix::wilcox_test(formula = num_var~long_var, paired = T, detailed = T) %>%
                rstatix::add_xy_position(x = "long_var", group = "cat_var")

            }
            return(list(stats = stats, test = wilc))
          })

      })
  }

}


get_comp_boxplots <-
  function(dat,
           cat_var,
           num_var,
           long_var,
           link_var,
           pal = "Accent",
           y_trans = c("identity")) {

    # set NULL transformation

      trans_fun <- function(x){do.call(y_trans, list(x))}




    df <-
      dat  %>%
      dplyr::select(SampleID,
                    link_var = !!sym(link_var),
                    long_var = !!sym(long_var),
                    cat_var = !!sym(cat_var),
                    num_var = !!sym(num_var)) %>%
      dplyr::filter(!is.na(cat_var))


    test_cat <-
      df %>%
      dplyr::group_by(long_var) %>%
      rstatix::wilcox_test(num_var ~ cat_var, p.adjust.method = "none") %>%
      rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
      rstatix::add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>%
      rstatix::add_xy_position("long_var", ) %>%
      mutate(y.trans = trans_fun(y.position))



    test_long <-
      df %>%
      PairedWilcox_withNAs(
        myDF = .,
        IDs = "link_var",
        Variable = "num_var",
        Response = "long_var",
        Subsets = "cat_var"
      ) %>%
      group_by(cat_var) %>%
      rstatix::adjust_pvalue(p.col = "p",
                             output.col = "p.adj",
                             method = "BH") %>%
      ungroup() %>%
      add_significance(p.col = "p.adj") %>%
      rstatix::add_x_position(x = c("long_var"), group = "cat_var") %>%
      mutate(y.trans = trans_fun(y.position))


    boxplot <-
      df %>%
      ggplot(., aes(x = as.factor(long_var), y = trans_fun(num_var))) +
      geom_boxplot(aes(fill = cat_var), position = "dodge", outlier.size = .6) +

      # coord_trans(y="sqrt")+
      # scale_y_continuous(labels = scales::trans_) +
      ggpubr::stat_pvalue_manual(
        data = test_cat,
        label = "p.adj.signif",
        hide.ns = T,
        size = 6,
        step.increase = 0.1,
        bracket.size = .75,
        y.position = "y.trans"
      ) +
      ggpubr::stat_pvalue_manual(
        data = test_long,
        label = "p.adj.signif",
        hide.ns = T,
        step.increase = 0.1,
        size = 6,
        bracket.size = .75,
        y.position = "y.trans"

      ) +
      theme_bw() +
      scale_fill_brewer(palette = pal) +
      labs(y = num_var, x = "weeks", fill = cat_var) +
      theme(
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.key.size = unit(1.3, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
      )

    return(list(test = list(cat = test_cat, long = test_long), plot = boxplot))


  }
    # Performs both categorical and longitudinal tests and places the in a boxplot.








