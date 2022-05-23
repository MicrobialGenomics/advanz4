create_LMM <- function(metadata, num_var, cat_var, long_var, link_var,breakpoints){

  myLMM <-metadata %>%
    as.data.frame() %>%
    dplyr::select(link_var = !!sym(link_var),
                  num_var = !!sym(num_var),
                  cat_var = !!sym(cat_var),
                  long_var = !!sym(long_var)) %>%
    dplyr::mutate(cat_var = as.factor(cat_var)) %>%
    myLMM_two.piece.GLMM(Datos = .,
                         Y = "num_var",
                         Time = "long_var",
                         Factor = "cat_var",
                         t = breakpoints,
                         ID = "link_var")


  myLMM
}


get_LMM_GraphParams <- function(Table){

  myList <- Table %>%
    pull(cat_var) %>%
    unique() %>%
    set_names() %>%
    map_dfr( ~ {

      model<-Table %>%
        dplyr::select(num_var, long_var, cat_var, link_var) %>%
        dplyr::filter(cat_var == .x) %>%
        lmerTest::lmer(formula = "num_var ~ long_var + (1|link_var)") %>%
        summary() %>%
        magrittr::extract2("coefficients")

      graphparams<- Table %>%
        summarise(xmin=min(long_var),
                  xmax=max(long_var),
                  ymin = model[1,1] + xmin*model[2,1],
                  ymax = model[1,1] + xmax*model[2,1]) %>%
        mutate(cat_var = .x) %>%
        relocate(cat_var)
    }) %>%
    as.data.frame()


  return(myList)
}


get_lmm_effects <- function(data, cat_vector, num_vector, long_var){
  myUnifiedModelList <- cat_vector %>%
    set_names() %>%

    map(function(cat_var) {

      num_vector %>%
        set_names() %>%

        map(function(num_var) {
          dat <-
            data %>%
            dplyr::select(SampleID,
                          link_var = link_var,
                          long_var = long_var,
                          cat_var = cat_var,
                          num_var = num_var)

          cat_levels <-
            dat %>%
            dplyr::select(link_var, cat_var) %>%
            unique() %>%
            dplyr::group_by(cat_var) %>%
            dplyr::summarize(n=n()) %>%
            dplyr::filter(n>1,
                          !is.na(cat_var)) %>%
            dplyr::pull(cat_var)


          dat <-
            dat %>%
            dplyr::filter(cat_var %in% cat_levels,
                          !is.na(cat_var))


          if (length(!is.na(cat_levels)) >= 2){

            myUnifiedModel <- dat %>%
              lmerTest::lmer(formula = "num_var ~ long_var * cat_var + (1|link_var)", data = .)

            myGraphParams <- dat %>%
              get_LMM_GraphParams() %>%
              mutate(across(contains("x"), ~ as.numeric(.x)))

            myPval <-
              anova(myUnifiedModel) %>%

              # magrittr::extract2("coefficients") %>%
              as.data.frame() %>%
              dplyr::slice(3L) %>%
              magrittr::set_names(c("SumSq", "MeanSq", "df", "DenDF", "Fval","pval")) %>%
              dplyr::mutate(pval = round(pval, 3))

            myplot <- dat %>%
              ggplot(.,
                     aes(
                       x = long_var,
                       y = num_var,
                       group = link_var,
                       color = cat_var
                     )) +
              geom_point(alpha = 0.2) +
              geom_line(alpha = 0.2) +
              scale_x_continuous(breaks = unique(dat$long_var),
                                 limits = c(min(dat$long_var),
                                            max(dat$long_var))) +

              geom_segment(
                data = myGraphParams,
                aes(
                  x = xmin,
                  y = ymin,
                  xend = xmax,
                  yend = ymax,
                  color = cat_var
                ),
                inherit.aes = T,
                size = 2
              ) +
              scale_color_brewer(palette = "Set1") +
              labs(color = cat_var,
                   x = "weeks",
                   y=num_var) +
              geom_label(
                data = myPval,
                aes(
                  alpha = NULL,
                  label = paste("ANOVA p=",
                                pval,
                                sep = ""),
                  x = -Inf,
                  y = Inf
                ),
                color = "black",
                hjust = 0.01,
                vjust = 1
              ) +
              theme_bw()

            mySummary <-
              summary(myUnifiedModel) %>%
              magrittr::extract2("coefficients") %>%
              as.data.frame() %>%
              kableExtra::kable(., format = "markdown")

            return(list(plot = myplot, summary = mySummary))

          } else {return(NA)}
        })

    })

}




