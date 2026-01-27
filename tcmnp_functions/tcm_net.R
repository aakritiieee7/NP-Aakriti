#' Network interaction diagram of herbs, ingredients, and targets
#'
#' @param network.data data frame
#' must contain herb, molecule, target three columns of data
#' @param node.color node color
#' see "RColorBrewer::display.brewer.all()"
#' @param node.size  node size
#' @param label.size label size
#' @param label.degree
#' the node degree is the number of connections that
#' the node has with the other nodes.
#' Nodes with connections greater than or
#' equal to degree will be displayed.
#' @param edge.color edge color
#' @param edge.width edge width
#' @param graph.layout etwork Diagram Layout:
#' @param graph.layout Network Diagram Layout:
#' "kk", "nicely", "circle", "sphere",
#' "bipartite", "star", "tree", "randomly",
#' "gem", "graphopt","lgl", "grid",
#' "mds", "sugiyama","fr"
#' @param rem.dis.inter remove single free unconnected nodes
#' @param label.repel label repel
#'
#' @return Network Diagram
#' @export
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 aes
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr count
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_fan
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph scale_edge_width
#' @importFrom ggraph geom_edge_link0
#' @importFrom ggraph theme_graph
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph graph_from_data_frame
#' @examples
#' \dontrun{
#' data("xfbdf", package = "TCMNP")
#' network.data <- xfbdf %>%
#'   dplyr::select(herb, molecule, target) %>%
#'   sample_n(100, replace = FALSE) %>%
#'   as.data.frame()
#' tcm_net(network.data, node.color= "Spectral",
#' label.degree = 0, rem.dis.inter = TRUE,graph.layout = "fr",
#' label.size = 3)
#' }
tcm_net <- function(network.data,
                    node.size = c(2, 8),
                    label.size = 3.5,
                    label.repel = TRUE,
                    edge.color = "grey80",
                    edge.width = c(0.3, 1.5),
                    graph.layout = "fr",
                    label.degree = 0,
                    max.overlaps=10,
                    rem.dis.inter = FALSE){

  # -----------------------------
  # Data cleaning
  # -----------------------------
  network.data <- as.data.frame(network.data)

  network.data$herb     <- stringr::str_trim(as.character(network.data$herb))
  network.data$molecule <- stringr::str_trim(as.character(network.data$molecule))
  network.data$target   <- stringr::str_trim(as.character(network.data$target))

  network.data <- unique(network.data)

  # -----------------------------
  # Edge list with FIXED weights
  # -----------------------------
  links <- rbind(
    network.data %>%
      dplyr::select(herb, molecule) %>%
      dplyr::rename(from = herb, to = molecule) %>%
      dplyr::mutate(weight = 0.5),

    network.data %>%
      dplyr::select(molecule, target) %>%
      dplyr::rename(from = molecule, to = target) %>%
      dplyr::mutate(weight = 0.2)
  ) %>%
    dplyr::distinct()

  # -----------------------------
  # Nodes
  # -----------------------------
  nodes <- data.frame(
    name = unique(c(links$from, links$to)),
    stringsAsFactors = FALSE
  )

  net <- igraph::graph_from_data_frame(
    d = links,
    vertices = nodes,
    directed = FALSE
  )

  # -----------------------------
  # Node attributes
  # -----------------------------
  igraph::V(net)$degree <- igraph::degree(net)

  igraph::V(net)$class <- ifelse(
    igraph::V(net)$name %in% network.data$herb, "Herb",
    ifelse(igraph::V(net)$name %in% network.data$molecule,
           "Molecule", "Target")
  )

  # -----------------------------
# Plot
# -----------------------------
p <- ggraph::ggraph(net, layout = graph.layout) +

  # Edges
  ggraph::geom_edge_link0(
    aes(edge_linewidth = weight),
    edge_colour = edge.color,
    alpha = 0.35
  ) +

  # Nodes
  ggraph::geom_node_point(
    aes(size = degree, shape = class, color = class),
    alpha = 1
  ) +

  # Labels (ALL nodes)
  ggraph::geom_node_text(
    aes(label = name, color = class),
    size = label.size,
    repel = label.repel,
    max.overlaps = 10,
    show.legend = FALSE
  ) +

  # Class colors
  ggplot2::scale_color_manual(
    values = c(
      "Herb"     = "#009E73",
      "Molecule" = "#0072B2",
      "Target"   = "#D55E00"
    )
  ) +

  # Node size (NO degree legend)
  ggplot2::scale_size_continuous(
    range = node.size,
    guide = "none"
  ) +

  # Edge width (NO legend)
  ggraph::scale_edge_width(
    range = edge.width,
    guide = "none"
  ) +

  # Base theme
  ggraph::theme_graph() +

  # Legend styling + placement
  ggplot2::theme(
    legend.justification = c("left", "top"),

    legend.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    legend.box.background = ggplot2::element_rect(
      fill = "grey80",
      colour = NA
    ),

    legend.title = ggplot2::element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      size = 14,
      face = "bold"
    ),

    legend.key.size = grid::unit(1.4, "cm"),
    legend.spacing.y = grid::unit(0.6, "cm"),

    plot.margin = ggplot2::margin(20, 150, 20, 20)
  ) +

  # Allow drawing outside panel
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::guides(
  color = ggplot2::guide_legend(
    override.aes = list(size = 8) 
  ),
  shape = ggplot2::guide_legend(
    override.aes = list(size = 22)
  )
)

return(p)}
  