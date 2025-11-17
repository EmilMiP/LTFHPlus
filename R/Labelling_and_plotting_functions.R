utils::globalVariables("id1")
utils::globalVariables("id2")
utils::globalVariables("lab")
utils::globalVariables("gen.x")
utils::globalVariables("gen.y")
utils::globalVariables("k")
utils::globalVariables("N")
utils::globalVariables("half")
utils::globalVariables("x")
utils::globalVariables("xend")
utils::globalVariables("y")
utils::globalVariables("yend")
utils::globalVariables("reldegree")
utils::globalVariables("reported_label")
utils::globalVariables("id2.y")
utils::globalVariables("R")
utils::globalVariables("r")
utils::globalVariables("contained")
utils::globalVariables("path")
utils::globalVariables("shortest")
utils::globalVariables("pair")
utils::globalVariables("short_path")
utils::globalVariables("generations")
utils::globalVariables("fid")



#' Plot the (Average) Number of Relatives per Proband by Relationship Type
#'
#' @description
#' Produces a structured visualisation of the (average) number of each relationship
#' type per proband, based on the labelled pairwise relationship data from
#' [label_relatives()]. The plot arranges relationship types according to
#' generational distance and degree of relatedness, providing an intuitive overview
#' of kinship structure within the study sample.
#'
#' @param labelled_relations A tibble or data frame containing pairwise relationship
#'   labels and their associated metadata. Must include the following columns:
#'   \describe{
#'     \item{id1}{Identifier for the first individual (typically the proband).}
#'     \item{id2}{Identifier for the second individual (the relative).}
#'     \item{gen.x, gen.y}{Number of generations separating each individual from their most recent common ancestor.}
#'     \item{k}{Kinship coefficient between the two individuals.}
#'     \item{lab}{Relationship label assigned by [label_relatives()].}
#'   }
#' @param proband_vec A vector of identifiers for probands. The function restricts
#'    the analysis to pairs where \code{id1} is included in this vector.
#' @param reported_info Chose which information is reported on the figure.
#'    \describe{
#'    \item{total}{shows the total number of relatives of each type across all probands.}
#'    \item{average}{shows the average number of relatives of each type per proband.}
#'    \item{both}{(default) shows both total and average numbers on the plot.}
#'    }
#'
#' @details
#' If any relationship types in the input are not recognised in the predefined
#' mapping (e.g., rare or complex kinships), these are aggregated and shown as 
#' `"Other"`.
#'
#' @return
#' A \pkg{ggplot2} object showing the (mean) number of relatives per proband for each
#' relationship type. The plot can be further modified using standard
#' \pkg{ggplot2} functions (e.g., \code{+ theme()} or \code{+ labs()}).

#' @seealso
#' [label_relatives()] for generating the relationship labels used as input.
#'
#' @examples
#' # see vignette on identifying and labelling relatives
#' 
#' @importFrom ggplot2 annotate geom_label geom_segment geom_vline ggplot labs scale_fill_discrete theme_classic coord_equal scale_y_continuous theme element_text element_blank aes theme_bw
#' @importFrom dplyr full_join join_by
#' @export
Relation_per_proband_plot = function(labelled_relations, proband_vec, reported_info = "both") {
  if (!(reported_info %in% c("both", "total", "average")) ) {
    warning("reported_info must be one of 'both', 'total', or 'average'. Setting to default: 'both'")
    reported_info = "both"
  }
  
  # creating a padding_label for the boxes:
  padding_label = paste0(paste0(rep(" ", 14), collapse = ""), "\n  \n")
  
  # restrict labelled_relations to only the relations related to the probands
  # in proband_vec
  labelled_relations = labelled_relations %>%
    filter(id1 %in% proband_vec) %>%
    group_by(lab, gen.x, gen.y, k) %>%
    summarise(N = n(), .groups = "drop") %>%
    mutate(
      half = grepl("H", lab),
      y = gen.x - gen.y,
      x = pmin(gen.x, gen.y) * ifelse(half, -1, 1)
    )
  
  # number of probands
  nb_proband = length(proband_vec)
  
  # create dummy data to fill out the plot when no relationships of that type exist
  # center column labels:
  center = tibble(
    x = 0,
    y = c(5:1, 0, -(1:5)),
    lab = c(paste0(4:2, "GP"),"GP", "P", "Proband", "Ch", "GCh", paste0(2:4, "GCh")),
    k = 2^{-c(5:1, 0, 1:5)}/2
  )
  
  # first column, right:
  
  first_col = tibble(
    x = 1,
    y = c(4:0, -(1:4)),
    lab = c(paste0(3:2, "GPib"),"GPib", "Pib", "S", "Nib", "GNib", paste0(2:3, "GNib")),
    k = 2^{-c(5:1, 2:5)} /2
  )
  
  # second column, right:
  second_col = tibble(
    x = 2,
    y = c(3:0, -(1:3)),
    lab = c(paste0("1C", 3:1, "R"), "1C", paste0("1C", 1:3, "R")),
    k = 2^{-c(6:3, 4:6)} / 2
  )
  
  # third column, right:
  third_col = tibble(
    x = 3,
    y = c(2:0, -(1:2)),
    lab = c(paste0("2C", 2:1, "R"), "2C", paste0("2C", 1:2, "R")),
    k = 2^{-c(7:5, 6:7)} / 2
  )
  # fourth column, right:
  fourth_col = tibble(
    x = 4,
    y = 0,
    lab = "3C",
    k = 2^{-8} / 2
  )
  
  # combining all columns and creating the mirrored (half) relationships
  dummy_relations = bind_rows(center, first_col, second_col, third_col, fourth_col) %>%
    bind_rows(
      filter(., x != 0) %>%
        mutate(x = -x,
               lab = paste0("H", lab),
               k = k /2)
      
    ) %>%
    # add degree of relatedness
    mutate(reldegree = -log2(2*k))
  
  # Does any relationship type exist in the labelled_relations data that is not
  # present in the prespecified relations?
  
  # are there any labels in the input data, that does not appear in the dummy relations?
  if ( any( !(labelled_relations$lab %in% dummy_relations$lab) ) ) {
    # setting these unidentified relations to "Other"
    labelled_relations = labelled_relations %>%
      filter(lab %in% dummy_relations$lab) %>%
      bind_rows(tibble(
        lab = "Other",
        x = 4,
        y = 5,
        k = NA,
        N = labelled_relations %>% filter( !(lab %in% dummy_relations$lab) ) %>% pull(N) %>% sum)
      )
    
    # add "Other" to the plotted squares in dummy_relations
    dummy_relations = dummy_relations %>%
      bind_rows(tibble(
        x = 4,
        y = 5,
        lab = "Other",
        k = NA
      ))
  }
  
  # combining the dummy labeled data with the provided labelled data:
  plot_data = full_join(labelled_relations %>%
                          # calculate degree of relatedness
                          mutate(reldegree = -log2(2*k)),
                        dummy_relations,
                        by = join_by(lab, k, y, x, reldegree)) %>%
    # if N is not observed, set to 0
    mutate(N = if_else( is.na(N), 0L, N),
           avg_pr_proband = round(N/nb_proband,3),
           reldegree = if_else(reldegree == 0, NA, reldegree),
           reported_label = case_when(
             reported_info == "total"   ~ paste0(lab, "\n", N),
             reported_info == "average" ~ paste0(lab, "\n", avg_pr_proband),
             reported_info == "both"    ~ paste0(lab, "\n", N, "\n(", avg_pr_proband, ")")
           ))
  
  # initialise plot
  ggplot()+
    # add dashed vertical lines + full / half labels
    geom_vline(xintercept = -.5,linetype="dashed")+
    annotate(geom = "text",x = -2.5,y=4,label="Half")+
    geom_vline(xintercept = .5,linetype="dashed")+
    annotate(geom = "text",x = 2.5,y=4,label="Full")+theme_bw()+
    # add kinship structure lines
    geom_segment(data=bind_rows(
      tibble( # vertical lines - center + below generation 0
        x=    -3:3,
        xend= -3:3,
        y=    c(-(2:5), -(4:2)),
        yend= c(rep(0, 3), 5, rep(0, 3))),
      tibble( # diagonal lines - non-double linked lines
        x=    c(-3, -4:-1,2:4, 3),
        xend= c(0,0,0,0,0,1,1,1, 1),
        y=    c(2, rep(0, 7), 2),
        yend= c(5:1,1:4)),
      tibble( # diagonal extra lines
        x=    .9,
        xend= 0,
        y=    0:4,
        yend= 1:5-0.1),
      tibble( # diagonal extra lines
        x=    1.1,
        xend= 0,
        y=    0:4,
        yend= 1:5+0.1)),
      aes(x=x, y=y, xend=xend, yend=yend), linewidth=1.3) +
    # add and color boxes for each relationship
    geom_label(data = plot_data,
               aes(x=x, y=y, fill=factor(reldegree)),
               label = padding_label,
               size=3, color="black") +
    # color proband box white
    annotate(geom  = "label", x=0, y=0, fill="white",
             label = padding_label,
             color = "black", size=3) +
    # add proband box text
    geom_label(data=plot_data,
               aes(x=x,
                   y=y,
                   label = reported_label),
               fill=NA,label.size=NA,size=3) +
    # theme
    theme_classic() +
    # remove ticks on x-axis
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    #
    coord_equal() +
    labs(x = "", y = "Generation", fill = "Degree of \nRelatedness",
         title = "Average Relationship per Proband") +
    scale_fill_discrete(na.value = "grey90", na.translate = FALSE) +
    scale_y_continuous(breaks = c(-4, -2, 2, 4))
}



#' Label Pairwise Relationships Based on Generational Distance and Kinship Coefficient
#'
#' @description
#' Assigns standard pedigree relationship labels (e.g., *Parent*, *Child*, *Sibling*, *Grandparent*, *Cousin*)
#' to all pairs of individuals based on their generational distances (`gen.x`, `gen.y`)
#' and kinship coefficients (`k`), typically produced by [get_generations()].
#'
#' @param tbl A tibble or data frame containing at least the following columns:
#' \describe{
#'   \item{fid}{Column with family identifier (typically the proband's id).}
#'   \item{id1}{Identifier for the first individual.}
#'   \item{id2}{Identifier for the second individual.}
#'   \item{gen.x}{Number of generations between `id1` and their most recent common ancestor with `id2`.}
#'   \item{gen.y}{Number of generations between `id2` and their most recent common ancestor with `id1`.}
#'   \item{k}{Estimated kinship coefficient between the two individuals.}
#' }
#'
#' @details
#' This function derives descriptive relationship labels using generational differences
#' and kinship patterns. The labels are written in a short-hand notation, an explaination of a subset is given below:
#' \itemize{
#'   \item *P* - Parent
#'   \item *Ch* - Child
#'   \item *S* - Sibling
#'   \item *GP* - Grandparent
#'   \item *Pib* - "Pibling" (parental sibling; aunt/uncle)
#'   \item *Nib* - "Nibling" (sibling's child; niece/nephew)
#'   \item *GCh* - Grandchild
#'   \item *GPib* - Grandpibling (grandparent's sibling)
#'   \item *GNib* - Grandnibling (sibling's grandchild)
#'   \item *C* - Cousin
#'   \item *1C1R* - First Cousin Once Removed
#'   \item *2C2R* - Second Cousin Twice Removed
#'   \item *H* prefix - Half relationships (e.g., *HS* for Half-Sibling)
#' }
#'
#' @return
#' A tibble with the following columns:
#' \describe{
#'   \item{fid}{Column with family identifier (typically the proband's id).}
#'   \item{id1}{Identifier for the first individual.}
#'   \item{id2}{Identifier for the second individual.}
#'   \item{gen.x}{Generational distance for \code{id1}.}
#'   \item{gen.y}{Generational distance for \code{id2}.}
#'   \item{k}{Kinship coefficient between the two individuals.}
#'   \item{lab}{Assigned relationship label (e.g., `"S"`, `"P"`, `"1C"`, `"H1C"`, `"2GP"`, etc.).}
#' }
#'
#' @seealso
#' [get_generations()] for computing the generational and kinship inputs used by this function.
#'
#' @examples
#' # see vignette on identifying and labelling relatives
#'
#' @importFrom dplyr if_else select mutate filter 
#' 
#' @export
label_relatives <- function(tbl) {
  tbl %>%
    select(fid, id1, id2, gen.x, gen.y, k) %>%
    mutate(
      R = k * 2,
      r = if_else(gen.x == 0 | gen.y == 0,
                  0.5^(gen.x + gen.y),
                  0.5^(gen.x + gen.y) * 2),
      half = (r / 2 == R),
      lab = paste0(pmin(gen.x, gen.y) - 1, "C", abs(gen.y - gen.x), "R"),
      lab = if_else(pmin(gen.x, gen.y) - 1 == 0, paste0(gen.x - 2, "GPib"), lab),
      lab = if_else(lab == "0GPib", "Pib", lab),
      lab = if_else(lab == "1GPib", "GPib", lab),
      lab = if_else(lab == "-1GPib", paste0(gen.y - 2, "GNib"), lab),
      lab = if_else(lab == "0GNib", "Nib", lab),
      lab = if_else(lab == "1GNib", "GNib", lab),
      lab = if_else(gen.y == 0 & gen.x > 2, paste0(gen.x - 1, "GP"), lab),
      lab = if_else(gen.y == 0 & gen.x == 2, "GP", lab),
      lab = if_else(gen.y == 0 & gen.x == 1, "P", lab),
      lab = if_else(gen.y == 1 & gen.x == 0, "Ch", lab),
      lab = if_else(gen.y == 2 & gen.x == 0, "GCh", lab),
      lab = if_else(gen.y > 2 & gen.x == 0, paste0(gen.y - 1, "GCh"), lab),
      lab = if_else(gen.y == gen.x & !(gen.x == 1), paste0(gen.x - 1, "C"), lab),
      lab = if_else(gen.y == gen.x & gen.x == 1, "S", lab),
      lab = if_else(gen.y == gen.x & gen.x == 0, "Proband", lab),
      lab = if_else(half, paste0("H", lab), lab)
    ) %>%
    filter(r == R | half) %>%
    select(fid, id1, id2, gen.x, gen.y, k, lab)
}


#' Compute Generational Distances and Kinship Coefficients from a Family Graph
#'
#' @description
#' Calculates generational distances and kinship coefficients between all pairs of individuals
#' represented in a directed family graph. The function identifies shortest paths between
#' individuals, accounts for common ancestors, and derives kinship coefficients based on
#' the number of generations separating each pair.
#'
#' @param fam_graph An \link[igraph]{igraph} object representing the family structure,
#' where directed edges indicate parent–child relationships (from parent to child). See get_family_graphs().
#'
#' @return
#' A tibble with one row per unique pair of individuals and the following columns:
#' \describe{
#'   \item{id1}{Identifier for the first individual.}
#'   \item{id2}{Identifier for the second individual.}
#'   \item{k}{Estimated kinship coefficient between \code{id1} and \code{id2}.}
#'   \item{gen.x}{Mean number of generations separating \code{id1} from the common ancestor.}
#'   \item{gen.y}{Mean number of generations separating \code{id2} from the common ancestor.}
#' }
#' The family graph is centered on a proband (typically the same id as fid),
#' all relations for the proband can be found by selecting only relations with the
#' proband's id in the column id1.
#'
#' @examples
#' # see vignette on identifying and labelling relatives
#'
#' @importFrom dplyr anti_join inner_join  
#' @importFrom utils tail
#' @importFrom igraph delete_edges E V all_shortest_paths which_mutual vertex_attr
#' @export
get_generations <- function(fam_graph) {
  ### function originally developed by Morten Krebs and added here with his permission:
  ### https://github.com/BioPsyk/PAFGRS/blob/master/R/kinship_sparse.R
  ### Minor changes made by Emil Pedersen so it accept a graph input directly
  
  # ids - intended to include everyone
  ids = igraph::V(fam_graph)$name
  
  #remove any duplicate edges in the graph, if present
  fam_graph = delete_edges(fam_graph, E(fam_graph)[which_mutual(fam_graph)])
  
  # measure all distances between family members
  anc <- Reduce("rbind.data.frame", lapply(vertex_attr(fam_graph)[[1]],
                                           function(i) {
                                             p <- all_shortest_paths(fam_graph, from = i, mode = "in")$res
                                             if (length(p) != 0)
                                               cbind.data.frame(id1 = as.numeric(i), id2 = as.numeric(names(sapply(p,
                                                                                                                   tail, n = 1))), gen = sapply(p, length) -
                                                                  1, path = I(lapply(p, function(x) as.numeric(names(x)))))
                                           }))
  # removes distances to self
  anc <- anc[!anc$id1 == anc$id2, ]
  # getting distances in the other direction
  dec <- setNames(anc, c("id2", "id1", "gen", "path"))
  
  # only one person was in the family graph:
  # returning just self-relations
  if (nrow(anc) == 0) {
    return(tibble(id1=as.character(ids),id2=as.character(ids),k=0.5,gen.x=0,gen.y=0)) # only self-relations
  }
  
  # combining all distances to one object
  other = inner_join(anc, dec, by = c("id2" = "id1"), relationship = "many-to-many")
  
  # modifying anc and dec for later use
  anc$gen.x = 0
  dec$gen.y = 0
  names(anc)[3] <- "gen.y"
  names(dec)[3] <- "gen.x"
  
  # formatting other and removing any duplicate pairs
  other = other %>%
    select(-id2) %>%
    rename(id2 = id2.y) %>%
    filter(!(id1 == id2)) %>%
    mutate(pair = paste(id1, id2, sep = "_")) %>%
    anti_join(
      dec %>% mutate(pair = paste(id1, id2, sep = "_")) %>% select(pair),
      by = "pair"
    ) %>%
    anti_join(
      dec %>% mutate(pair = paste(id2, id1, sep = "_")) %>% select(pair),
      by = "pair"
    ) %>%
    select(-pair)
  
  # construct full paths, with no double counting of shared ancestor
  other$path <- Map(function(px, py) c(px, rev(py)[-1]), other$path.x, other$path.y)
  
  # identify shortest paths per id1-id2 pair
  other = other %>%
    group_by(id1, id2) %>%
    mutate(shortest = (gen.x + gen.y) == min(gen.x + gen.y)) %>%
    ungroup()
  
  # get list of shortest paths per id1-id2 pair
  other_shortest = other %>%
    filter(shortest) %>%
    group_by(id1, id2) %>%
    summarise(short_path = list(path), .groups = "drop")
  
  # join shortest paths back to other
  other = left_join(other, other_shortest, by = c("id1", "id2"))
  
  # for non-shortest paths, check if any shortest path is fully contained within
  other_not_shortest = other %>% filter(!shortest) %>%
    mutate(contained = purrr::map2_lgl(path, short_path, function(p, short_list) {
      # Return TRUE if any shortest path is fully contained within p
      any(purrr::map_lgl(short_list, function(s) all(s %in% p)))
    }))
  
  # keep only shortest paths initially
  other = other %>% filter(shortest)
  
  # add non-shortest paths that are not contained within any shortest path
  not_contained = other_not_shortest %>% filter(!contained)
  
  # for these, keep only the shortest per pair
  if (nrow(not_contained) > 0) {
    not_contained = not_contained %>%
      group_by(id1, id2) %>%
      filter((gen.x + gen.y) == min(gen.x + gen.y)) %>%   # keep shortest per pair
      ungroup() %>%
      select(-contained)
    
    other = bind_rows(other, not_contained)
  }
  
  # calculate kinship coefficients from paths
  # an set gen.x and gen.y
  k <- bind_rows(
    anc %>% select(id1, id2, gen.y = gen.x, gen.x = gen.y, path),
    dec %>% select(id1, id2, gen.y = gen.x, gen.x = gen.y, path),
    other %>% select(id1, id2, gen.x, gen.y, path)
  ) %>%
    mutate(r = 0.5^gen.x * 0.5^gen.y) %>%
    group_by(id1, id2) %>%
    summarise(
      k = sum(r) / 2,
      gen.x = mean(gen.x),
      gen.y = mean(gen.y),
      .groups = "drop"
    )
  # Check types of ids and convert if necessary
  type_mismatch <- !(
    typeof(ids) == typeof(k$id1) && typeof(ids) == typeof(k$id2)
  )
  # If there is a type mismatch, convert k$id1 and k$id2 to the type of ids
  if (type_mismatch) {
    target_type <- typeof(ids)
    if (target_type == "character") {
      k$id1 <- as.character(k$id1)
      k$id2 <- as.character(k$id2)
    } else if (target_type %in% c("double", "integer")) {
      k$id1 <- as.numeric(k$id1)
      k$id2 <- as.numeric(k$id2)
    } else {
      stop("Unsupported type for id columns")
    }
  }
  # add self-kinship entries
  k = k %>% bind_rows(
    tibble(id1 = ids, id2 = ids, k = 0.5, gen.x = 0, gen.y = 0)
  )
  # return result
  return(k)
}


#' Compute and Label Pairwise Relationships Across Multiple Family Graphs
#'
#' @description
#' Applies [get_generations()] and [label_relatives()] to a tibble of
#' family graph objects from [get_family_graphs()], returning a unified table of labelled pairwise relationships
#' for all individuals across the specified families.
#' 
#' @param family_graphs A tibble containing family-specific graph objects 
#'  from [get_family_graphs()] (typically of class `igraph`). Each row should correspond to a
#'   distinct family, with one column containing the graph object and another
#'   containing the family identifier (typically the proband's id).
#' @param fid A character string specifying the name of the column in 
#'    `family_graphs` that holds the family identifiers. Defaults to `"fid"`.
#' @param family_id_vec An optional character or numeric vector specifying
#'   which families to process. If `NULL` (default), the function will process
#'   all families in `family_graphs`.
#' @param fam_graph_col A character string specifying the name of the column
#'   in `family_graphs` that contains the family graph objects. Defaults to `"fam_graph"`.
#'   
#' @return
#' A tibble containing all labelled pairwise relationships across the specified
#' families, with columns:
#' \describe{
#'   \item{fid}{Family identifier (typically the proband's id).}
#'   \item{id1, id2}{Identifiers for the two individuals being compared.}
#'   \item{gen.x, gen.y}{Number of generations separating each individual
#'         from their most recent common ancestor.}
#'   \item{k}{Kinship coefficient between the pair.}
#'   \item{lab}{Relationship label (e.g., `"S"`, `"P"`, `"GP"`, `"1C"`, `"HNib"`).}
#' }
#'
#' @seealso
#' [get_generations()] for computing generational distances and kinship coefficients,  
#' [label_relatives()] for labelling relationships based on generational patterns.
#'
#' @examples
#' # See vignette
#' 
#' @export
get_relations = function(family_graphs, fid = "fid", family_id_vec = NULL, fam_graph_col = "fam_graph") {
  # get generations for each family graph in family_id_vec
  if ( is.null(family_id_vec) ) {
    res = family_graphs %>%
      mutate(generations = lapply(!!as.symbol(fam_graph_col),
                               function(cur_fam_graph) get_generations(cur_fam_graph))) %>%
      select(!!as.symbol(fid), generations)

  } else {
    res = family_graphs %>% 
      filter(!!as.symbol(fid) %in% family_id_vec) %>%
      mutate(generations = lapply(!!as.symbol(fam_graph_col),
                               function(cur_fam_graph) get_generations(cur_fam_graph))) %>%
      select(!!as.symbol(fid), generations)   

  }
  # concatenate all results and label relatives
  res %>% tidyr::unnest(generations) %>% label_relatives()
}
