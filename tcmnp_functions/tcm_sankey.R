#' Sankey diagram
#'
#' @param data_sankey data.frame
#' @param text.size text size
#' @param text.position text position (kept for compatibility)
#' @param x.axis.text.size Location of x-axis labels in sankey plot
#' @param ... additional parameters
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggplot2 ggplot aes theme scale_fill_manual
#' @importFrom ggplot2 element_text element_blank
#' @importFrom ggplot2 position_nudge
#' @importFrom dplyr mutate
#' @importFrom ggsankey geom_sankey geom_sankey_text
#' @importFrom ggsankey theme_sankey
#' @importFrom cols4all c4a

tcm_sankey <- function(data_sankey,
                       text.size = 3,
                       text.position = 0,
                       x.axis.text.size = 14,
                       ...) {

  library(ggplot2)
  library(dplyr)
  library(ggsankey)
  library(cols4all)

  ## ---- rename columns ----
  colnames(data_sankey) <- c("Plant", "Phytochemical", "Target")

  ## ---- build sankey dataframe ----
  sankey_df <- do.call(
    rbind,
    apply(data_sankey, 1, function(x) {
      data.frame(
        x = names(x),
        node = x,
        next_x = dplyr::lead(names(x)),
        next_node = dplyr::lead(x),
        stringsAsFactors = FALSE
      )
    })
  ) %>%
    mutate(
      x = factor(x, names(data_sankey)),
      next_x = factor(next_x, names(data_sankey))
    )

  ## ---- colors ----
  df_long <- ggsankey::make_long(data_sankey, colnames(data_sankey))
  cols <- cols4all::c4a("rainbow_wh_rd", length(unique(df_long$node)))
  cols <- sample(cols)

  ## ---- plot ----
  p <- ggplot(
    sankey_df,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = node,
      label = node
    )
  ) +
    ggsankey::geom_sankey(
      flow.alpha = 0.5,
      node.color = NA,
      node.width = 0.08,
      show.legend = FALSE
    ) +
    ggsankey::geom_sankey_text(
      size = text.size,
      color = "black",
      hjust = -0.1,  # Negative value pushes text RIGHT of bars
      position = position_nudge(x = 0)
    ) +
    scale_fill_manual(values = cols) +
    theme_sankey(base_size = 18) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(
        size = x.axis.text.size,
        face = "plain",
        colour = "black"
      ),
      plot.margin = margin(10, 150, 10, 150)
    ) +
    coord_cartesian(clip = "off")

  return(p)
}