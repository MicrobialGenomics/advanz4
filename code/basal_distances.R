distance_from_basal <- function(dist_mat, metadata, link_var, long_var, cat_var, basal = 0) {

  basal <-
    dplyr::pull(metadata, !!sym(long_var)) %>%
    min()

  dist_df <-
    as.matrix(dist_mat) %>%
    tibble::as_tibble(rownames = "SampleID") %>%
    tidyr::pivot_longer(-1, names_to = "basal_id", values_to = "distance")

  dat <-
    metadata %>%
    dplyr::select(SampleID, !!link_var, !!long_var, !!cat_var) %>%
    dplyr::group_by(!!sym(link_var)) %>%
    dplyr::filter(any(!!sym(long_var) == basal) & n() > 1) %>%
    dplyr::mutate(basal_id = SampleID[!!sym(long_var) == basal], .before = 2) %>%
    dplyr::arrange(!!sym(link_var), !!sym(long_var)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(dist_df, by = c("SampleID", "basal_id")) %>%
    dplyr::filter(!!sym(long_var) != basal)

  plt <-
    ggplot(dat, aes(
      !!sym(long_var),
      distance,
      group = !!sym(link_var),
      colour = !!sym(cat_var),
      fill = !!sym(cat_var)
    )) +
    facet_grid(cols = vars(!!sym(cat_var))) +
    geom_point() +
    geom_line() +
    stat_summary(aes(group = !!sym(cat_var)), fun = median, geom = "line", size = 2) +
    geom_smooth(aes(group = !!sym(cat_var)), method = "glm", formula = 'y ~ x', colour = "black") +
    theme_bw()

  plt
}
