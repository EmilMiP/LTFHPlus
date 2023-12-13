utils::globalVariables("cip")
utils::globalVariables("all_combinations")
utils::globalVariables("children")
utils::globalVariables("comb")
utils::globalVariables("distances")
utils::globalVariables("from")
utils::globalVariables("nchildren")
utils::globalVariables("parent_id")
utils::globalVariables("ph")
utils::globalVariables("to")

#' Attempts to convert the list entry input format to a long format
#' 
#' @param family a tibble with two entries, family id and personal id. personal id should end in "_role", if a role column is not present.
#' @param threshs thresholds, with a personal id (without role) as well as the lower and upper thresholds
#' @param personal_id_col column name that holds the personal id
#' @param role_col column name that holds the role
#' 
#' @return returns a format similar to \code{prepare_LTFHPlus_input}, which is used by \code{estimate_liability}
#'
#' @export 

convert_format = function(family, threshs, personal_id_col = "pid", role_col = NULL) {
  # standardising input -----------------------------------------------------
  
  #are there any list columns in family (list entry format)?
  which_list_columns = which(sapply(family, is.list))
  if (length(which_list_columns) > 0) { 
    #updating fam to get rid of list columns
    family = tidyr::unnest(family, cols = names(which_list_columns))
    
    ###  checking if role is present in ID or separate column
    # if "_" is present, a role will be there too.
    if (any(stringr::str_detect(family[[personal_id_col]], "_"))) { #if true, extract role
      #split pid_role in two:
      family[[role_col]] = strsplit(family[[personal_id_col]], "_(?=[^_]+$)", perl=TRUE) %>% sapply(., function(x) x[2])
      family[[personal_id_col]] = strsplit(family[[personal_id_col]], "_(?=[^_]+$)", perl=TRUE) %>% sapply(., function(x) x[1])
      
      cat("We've tried converting from list entries to a long format internally. If you see this print, please run prepare_LTFHPlus_input and use the .tbl input going forward! \n")
    } else {
      if (is.null(role_col)) stop("Please provide family roles for each family member. e.g. father(f), mother(m), siblings (s1-s9), etc.") 
      stop("We weren't able to convert data input automatically. Please use prepare_LTFHPlus_input and use the .tbl input!\n")
    }
  }
  .tbl = left_join(family, threshs, by = personal_id_col)
  return(.tbl)
}

#' Prepares input for \code{estimate_liability}
#' 
#' @param family contains family and personal ids and role with a family.
#' @param CIP tibble with population representative cumulative incidence proportions has the interpretation of "proportion of people that has experienced the trait subset by \code{CIP_columns}. 
#' @param CIP_columns the columns the CIPs are subset by, e.g. CIPs by birth_year, sex. 
#' @param status_col Column that contains the status of each family member 
#' @param use_fixed_case_thr Should the threshold be fixed for cases? Can be used if CIPs are detailed, e.g. stratified by birth_year and sex.
#' @param fam_id_col Column that contains the family ID
#' @param personal_id_col Column that contains the personal ID
#' @param role_col Column that cotnains the role of each individual
#' 
#' @export

prepare_LTFHPlus_input = function(family, CIP, 
                                  CIP_columns = c("sex", "birth_year", "age"), 
                                  status_col = "status", 
                                  use_fixed_case_thr = F, 
                                  fam_id_col = "fam_id", 
                                  personal_id_col = "pid",
                                  role_col = "role") {

  

# Merging CIPs and assigning thresholds -----------------------------------

  family = dplyr::left_join(family, CIP, by = CIP_columns) %>% 
    dplyr::mutate(
      lower = ifelse(!!as.symbol(status_col) == 1, stats::qnorm(cip, lower.tail = F), -Inf),
      upper = ifelse(!!as.symbol(status_col) == 1, ifelse(use_fixed_case_thr, stats::qnorm(cip, lower.tail = F), Inf), stats::qnorm(cip, lower.tail = F))
    )
  

# returning formatted input -----------------------------------------------

  dplyr::select(family, !!as.symbol(fam_id_col), !!as.symbol(personal_id_col), !!as.symbol(role_col), lower, upper)
}



#' Construct graph from register information
#' 
#' \code{prepare_graph} constructs a graph based on mother, father, and offspring links. 
#'
#' @param .tbl tibble with columns icol, fcol, mcol. Additional columns will be attributes in the constructed graph.
#' @param icol column name of column with proband ids.
#' @param fcol column name of column with father ids.
#' @param mcol column name of column with mother ids.
#' @param missingID_patterns string of missing values in the ID columns. Multiple values can be used, but must be separated by "|". Defaults to "0".
#' 
#' @return An igraph object. A (directed) graph object based on the links provided in .tbl. 
#' 
#' @export
#' @importFrom dplyr %>% rename relocate mutate filter group_by summarise select bind_rows pull


prepare_graph = function(.tbl, fcol, mcol, icol, missingID_patterns = "^0$") {
  
  # formatting .tbl from trio info to graph compatible input
  prep = .tbl %>% 
    # making from column
    tidyr::pivot_longer(cols = c(!!as.symbol(fcol), !!as.symbol(mcol)), 
                        values_to = "from") %>% 
    # renaming id to "to"
    rename(to = !!as.symbol(icol)) %>% 
    # reloacting to and from columns to first two columns
    relocate(from, to) %>% # directed graph -> order is important! 
    # replacing "0"s with NA to ease later computation and reflect true data set
    # with missing / unknown links
    mutate(to = ifelse(str_detect(to, missingID_patterns), NA, to),
           from = ifelse(str_detect(from, missingID_patterns), NA, from)) 
  
  # remove connections with unknown links, i.e. only known links / edges
  parent_links = prep %>% 
    filter(!is.na(to), !is.na(from)) %>% 
    # ensure from and to are character vectors
    mutate(from = as.character(from),
           to = as.character(to))
  
  # we need to add a direct link between (full) siblings
  sibling_links = .tbl %>% 
    filter(str_detect(!!as.symbol(fcol), missingID_patterns, negate = TRUE), 
           str_detect(!!as.symbol(mcol), missingID_patterns, negate = TRUE)) %>% 
    mutate(parent_id = purrr::map2_chr(
      .x = !!as.symbol(fcol),
      .y = !!as.symbol(mcol),
      ~ paste0(sort(c(.x, .y)), collapse = "_")
    )) %>% 
    group_by(parent_id) %>% 
    summarise(children = list(!!as.symbol(icol))) %>% 
    mutate(nchildren = sapply(children, length)) %>% 
    filter(nchildren > 1) %>% 
    mutate(all_combinations = purrr::map(.x = children, ~ get_all_combs(.x))) %>% 
    # avoiding many duplicate rows by selecting only columns to unnest
    select(all_combinations) %>% 
    tidyr::unnest(cols = c(all_combinations)) %>% 
    mutate(ph = str_split(all_combinations, "_"),
           from = sapply(ph, function(x) x[1]),
           to = sapply(ph, function(x) x[2])) %>% 
    select(from, to)
  
  # combining sibling and parent links
  if (nrow(sibling_links) > 0) {
    graph_input = bind_rows(parent_links, sibling_links)
  } else {
    graph_input = parent_links
  }
  
  
  graph = igraph::graph_from_data_frame(d = graph_input)
  
  # isolating potential solo nodes: where to or from column is NA
  # extracting only node names
  solo = prep %>% 
    mutate(comb = purrr::map2_chr(.x = from, 
                                  .y = to,
                                  ~ paste0(sort(c(.x, .y)), collapse = "_"))) %>% 
    filter(str_detect(comb, "_", negate = TRUE)) %>% 
    pull(comb)
  
  # all linked nodes; no NAs in to or from,
  # just a list of node names
  duos = graph_input %>% 
    filter(!is.na(from) & !is.na(to)) %>% 
    select(from, to) %>% unlist() %>% unique()
  
  # which nodes appear only in solo nodes, but not in duo nodes?
  # meaning they have no links
  solo_points = setdiff(solo, as.character(duos))
  
  igraph::add.vertices(graph, nv = length(solo_points), name = solo_points) 
}

