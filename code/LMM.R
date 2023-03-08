#' Creates LMMs separated for each group.
#' @metadata data frame containing the variables of choice.
#' @cat_var a character containing the name of the categorical variables you want to split the model with.
#' @num_var character vector. Contains the names of the numerical variables (usually the Y variable)
#' @long_var a character, name of the Time or longitudinal variable. Alternatively, just the variable yu want to put on the X axis.
#' @link_var link variable eg. patient, mouse... it stands for the random variable.
#' @breakpoints In case the study has a breakpoint, like a discontinuation of the treatment and washout period. It will break the model at this point of the X axis and create a ner model onwards.

create_LMM <- function(data, num_var, cat_var, long_var, link_var,breakpoints = NULL){

  myLMM <-
    data %>%
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

  stats <- sapply(myLMM[[1]], function(x){x[2,c(1,5)]}) %>%
    t() %>%
    as.data.frame() %>%
    mutate(CategoricalVariable=rownames(.)) %>%
    setNames(c("slope","pval","cat_var")) %>%
    mutate(pval=round(pval,3),
           slope = round(slope, 3))




  plt <- myLMM[[2]] +
    geom_label(data = stats, inherit.aes = F,aes(label=paste("p= ",pval, ", slope= ", slope, sep=""),
                                x = 0, y = Inf),
               color = "black",
               hjust=0.01, vjust=1) +
    geom_smooth(aes(group=cat_var, fill=cat_var),method=glm ,alpha=0.4, linetype=0)+

    labs(x=long_var, title=num_var, fill=long_var)

  return(list(plot = plt, stats = stats))
}

#' Helper function to obtain the graphical parameters to build a beautiful LMM plot.
#' @data is the same data frame as the other functions.

get_LMM_GraphParams <- function(data, t) {

  # Summary of how the old function gets the breakpoints:
  # It creates a new variable (T1) which is basically

  data %>%
    pull(cat_var) %>%
    unique() %>%
    set_names() %>%
    map_dfr( ~ {
      model <- data %>%
        dplyr::select(num_var, long_var, cat_var, link_var) %>%
        dplyr::filter(cat_var == .x) %>%
        lmerTest::lmer(formula = "num_var ~ long_var + (1|link_var)") %>%
        summary() %>%
        magrittr::extract2("coefficients")

      data %>%
        summarise(
          xmin = min(long_var, na.rm = T),
          xmax = max(long_var, na.rm = T),
          ymin = model[1, 1] + xmin * model[2, 1],
          ymax = model[1, 1] + xmax * model[2, 1]
        ) %>%
        mutate(cat_var = .x) %>%
        relocate(cat_var)
    }) %>%
    as.data.frame()

}


#' A wrapper to perform unified lmms models which test for group effect in the slope. Test performed by ANOVA. It returns a nested list organized by cat_var/num_var and 2 slots: plot and statistics
#'
#' @data A data_frame/tibble object containing all the variables you want to test.
#' @cat_vector a character vector containing the names of the categorical variables you want to use.
#' @num_vector character vector. Contains the names of the numerical variables (usually the Y variable)
#' @long_var a character, name of the Time or longitudinal variable. Alternatively, just the variable yu want to put on the X axis.
#' @link_var link variable eg. patient, mouse... it stands for the random variable.
#'
get_lmm_effects <- function(data, cat_vector, num_vector, long_var, link_var, title = T){



   myUnifiedModelList <- cat_vector %>%
    set_names() %>%

    map(function(cat_var) {
      print(cat_var)
      num_vector %>%
        set_names() %>%

        map(function(num_var) {
          print(num_var)
          dat <-
            data %>%
            dplyr::select(SampleID,
                          link_var = !!sym(link_var),
                          long_var = !!sym(long_var),
                          cat_var = !!sym(cat_var),
                          num_var = !!sym(num_var))

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
              stats::setNames(c("estimate", "st.err", "df", "tval","p"))
              kableExtra::kable(., format = "markdown")

            if(title == T){
              myplot <-
                myplot +
                labs(title = num_var)
            } else {
              myplot <- myplot
            }


            return(list(plot = myplot, summary = mySummary))

          } else {return(NA)}
        })

    })
  return(myUnifiedModelList)
}




