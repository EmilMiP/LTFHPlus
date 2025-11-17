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
utils::globalVariables("event_age")
utils::globalVariables("cip_pred")
utils::globalVariables("thr")
utils::globalVariables("attrs")


#' Attempts to convert the list entry input format to a long format
#' 
#' @param family a tibble with two entries, family id and personal id. personal id should end in "_role", if a role column is not present.
#' @param threshs thresholds, with a personal id (without role) as well as the lower and upper thresholds
#' @param personal_id_col column name that holds the personal id
#' @param role_col column name that holds the role
#' 
#' @return returns a format similar to \code{prepare_LTFHPlus_input}, which is used by \code{estimate_liability}
#'
#' @examples
#' family <- data.frame(
#' fam_id = c(1, 1, 1, 1),
#' pid = c(1, 2, 3, 4),
#' role = c("o", "m", "f", "pgf")
#' )
#' 
#' threshs <- data.frame(
#'   pid = c(1, 2, 3, 4),
#'   lower = c(-Inf, -Inf, 0.8, 0.7),
#'   upper = c(0.8, 0.8, 0.8, 0.7)
#' )
#' 
#' convert_format(family, threshs)
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
      
      warning("We've tried converting from list entries to a long format internally. If you see this print, please run prepare_LTFHPlus_input and use the .tbl input going forward! \n")
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
#' @param .tbl contains family and personal ids and role with a family.
#' @param CIP tibble with population representative cumulative incidence proportions. CIP values should be merged by \code{CIP_columns}. 
#' @param age_col name of column with age
#' @param aoo_col name of column with age of onset
#' @param CIP_merge_columns The columns the CIPs are subset by, e.g. CIPs by birth_year, sex. 
#' @param CIP_cip_col name of column with CIP values
#' @param status_col Column that contains the status of each family member 
#' @param use_fixed_case_thr Should the threshold be fixed for cases? Can be used if CIPs are detailed, e.g. stratified by birth_year and sex.
#' @param fam_id_col Column that contains the family ID
#' @param personal_id_col Column that contains the personal ID
#' @param interpolation type of interpolation, defaults to NULL.
#' @param bst.params list of parameters to pass on to xgboost
#' @param min_CIP_value minimum cip value to allow, too low values may lead to numerical instabilities.
#' @param xgboost_itr Number of iterations to run xgboost for. 
#' 
#' 
#' @importFrom stats qnorm predict
#' @importFrom dplyr all_of mutate select %>% left_join group_by ungroup arrange
#' 
#' @return tibble formatted for \code{estimate_liability}
#' 
#' @examples
#' tbl = data.frame(
#'   fam_id = c(1, 1, 1, 1),
#'   pid = c(1, 2, 3, 4),
#'   role = c("o", "m", "f", "pgf"),
#'   sex = c(1, 0, 1, 1),
#'   status = c(0, 0, 1, 1),
#'   age = c(22, 42, 48, 78),
#'   birth_year = 2023 - c(22, 42, 48, 78),
#'   aoo = c(NA, NA, 43, 45))
#' 
#' cip = data.frame(
#'   age = c(22, 42, 43, 45, 48, 78),
#'   birth_year = c(2001, 1981, 1975, 1945, 1975, 1945),
#'   sex = c(1, 0, 1, 1, 1, 1),
#'   cip = c(0.1, 0.2, 0.3, 0.3, 0.3, 0.4))
#' 
#' prepare_LTFHPlus_input(.tbl = tbl,
#'                        CIP = cip, 
#'                        age_col = "age",
#'                        aoo_col = "aoo",
#'                        interpolation = NA)
#' 
#' @export
prepare_LTFHPlus_input = function(.tbl, 
                                  CIP, 
                                  age_col,
                                  aoo_col,
                                  CIP_merge_columns = c("sex", "birth_year", "age"), 
                                  CIP_cip_col = "cip", 
                                  status_col = "status", 
                                  use_fixed_case_thr = FALSE, 
                                  fam_id_col = "fam_id", 
                                  personal_id_col = "pid",
                                  interpolation = NULL,
                                  bst.params = list(
                                    max_depth = 10,
                                    base_score = 0,
                                    nthread = 4,
                                    min_child_weight = 10
                                  ),
                                  min_CIP_value = 1e-5,
                                  xgboost_itr = 50
                                  ) {
  #### checks for the presence of all columns in .tbl and CIP goes here ####
  #### checks for the presence of all columns in .tbl and CIP goes here ####
  
  
  # interpolation with xgboost or merge on raw values?
  if (!is.na(interpolation) && interpolation != "xgboost") stop("Invalid choice of interpolation method. Must be NULL or xgboost.")

  # interpolate CIP values based on xgb
  if (is.na(interpolation)) {
    
    # merge on raw values
    
    # Merging CIPs and assigning thresholds -----------------------------------
    .tbl = .tbl %>% dplyr::mutate(event_age = pmin(!!as.symbol(age_col), !!as.symbol(aoo_col), na.rm = T),
                                  !!as.symbol(age_col) := event_age) %>% 
      dplyr::select(-event_age) %>% 
      dplyr::left_join(CIP, by = CIP_merge_columns) %>% 
      dplyr::mutate(thr = qnorm(!!as.symbol(CIP_cip_col), lower.tail = FALSE),
                    lower = ifelse(!!as.symbol(status_col) == 1, thr, -Inf),
                    upper = ifelse(!!as.symbol(status_col) == 1, 
                                   ifelse(use_fixed_case_thr, thr, Inf),
                                   thr))
  
    
    if (any(is.na(.tbl$lower)) | any(is.na(.tbl$upper))) {
      warning(paste0("There are ", sum(is.na(select(.tbl, lower, upper))), " NA values in the upper and lower thresholds. \n Do the age and age of onset values match the ages given in the CIPs?"))
    }
    
    # returning formatted input -----------------------------------------------
    return(.tbl)
    

  } else if (interpolation == "xgboost") { 
    
    
    # extract cip values
    y = CIP[[CIP_cip_col]]
    
    # force the remaining values in cur_cip to be a matrix
    X = as.matrix(select(CIP, all_of(c(CIP_merge_columns, age_col)), -!!as.symbol(CIP_cip_col)))
    # train xgboost
    xgb = xgboost::xgboost(X, y, nrounds = xgboost_itr, params = bst.params)
    
    
    # get the predicted CIP value based on the merge columns.
    .tbl$cip_pred = .tbl %>% 
      mutate(age = pmin(!!as.symbol(age_col), !!as.symbol(aoo_col), na.rm = TRUE)) %>% 
      select(all_of(c(CIP_merge_columns, age_col))) %>% 
      as.matrix() %>% 
      predict(xgb,.) %>% 
      pmax(min_CIP_value)
    
    .tbl %>% 
      mutate(event_age = pmin(!!as.symbol(age_col), !!as.symbol(aoo_col), na.rm = TRUE)) %>% 
      group_by(across(all_of(setdiff(CIP_merge_columns, c(age_col, aoo_col))))) %>% 
      arrange(event_age) %>% 
      mutate(cip_pred = cummax(cip_pred)) %>% 
      ungroup() %>% 
      mutate(thr = qnorm(cip_pred, lower.tail = FALSE),
             lower = ifelse(!!as.symbol(status_col), thr, -Inf),
             upper = ifelse(!!as.symbol(status_col), 
                            ifelse(use_fixed_case_thr, thr, Inf),
                            thr)) #%>% select(!!as.symbol(fam_id_col), !!as.symbol(personal_id_col), !!as.symbol(role_col), lower, upper)
  } else {
    stop("unsupported interpolation method. Please use xgboost or NA.")
  }
}



#' Construct graph from register information
#'
#' \code{prepare_graph} constructs a graph based on mother, father, and offspring links.
#'
#' @param .tbl tibble with columns icol, fcol, mcol. Additional columns will be attributes in the constructed graph.
#' @param icol column name of column with proband ids.
#' @param fcol column name of column with father ids.
#' @param mcol column name of column with mother ids.
#' @param node_attributes tibble with icol and any additional information, such as sex, lower threshold, and upper threshold. Used to assign attributes to each node in the graph, e.g. lower and upper thresholds to individuals in the graph.
#' @param missingID_patterns string of missing values in the ID columns. Multiple values can be used, but must be separated by "|". Defaults to "^0$". OBS: "0" is NOT enough, since it relies on regex.
#'
#' @return An igraph object. A (directed) graph object based on the links provided in .tbl, potentially with provided attributes stored for each node.
#'
#' @importFrom dplyr %>% rename relocate mutate filter group_by summarise select bind_rows pull
#'
#' @examples
#' fam <- data.frame(
#'   id = c("pid", "mom", "dad", "pgf"),
#'   dadcol = c("dad", 0, "pgf", 0),
#'   momcol = c("mom", 0, 0, 0))
#'
#' thresholds <- data.frame(
#'   id = c("pid", "mom", "dad", "pgf"),
#'   lower = c(-Inf, -Inf, 0.8, 0.7),
#'   upper = c(0.8, 0.8, 0.8, 0.7))
#'
#' prepare_graph(fam, icol = "id", fcol = "dadcol", mcol = "momcol", node_attributes = thresholds)
#'
#' @export
prepare_graph = function(.tbl, icol, fcol, mcol, node_attributes = NA, missingID_patterns = "^0$") {
  
  # helper boolean to check if node_attributes is provided, since an offered node_attributes may have NA values in one or more
  # entries, hence the check on class of the object. NAs are logical and will lead to a FALSE.
  attachAttributes = any(class(node_attributes) %in% c("data.frame", "tibble", "matrix", "data.table", "tbl_df", "tbl"))
  
  # formatting .tbl from trio info to graph compatible input
  prep = .tbl %>%
    # making from column
    tidyr::pivot_longer(cols = c(!!as.symbol(fcol), !!as.symbol(mcol)),
                        values_to = "from") %>%
    # renaming id to "to"
    rename(to = !!as.symbol(icol)) %>%
    # reloacting to and from columns to first two columns
    select(from, to) %>% # directed graph -> order is important!
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
    select(!!as.symbol(icol), !!as.symbol(fcol), !!as.symbol(mcol)) %>%
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
  
  # extract unique list of ids of individuals in graph input; graph_input has only 2 columns.
  present_ids = unlist(graph_input) %>% unique()
  
  # if node_attributes is provided, we will attach the node information to the graph.
  if (attachAttributes) {
    graph = igraph::graph_from_data_frame(d = graph_input,
                                          #use unique id list to attach threshold info
                                          vertices = filter(node_attributes, !!as.symbol(icol) %in% present_ids))
  } else {
    graph = igraph::graph_from_data_frame(d = graph_input)
  }
  
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
  # only run the below code if solo_points has any solo points to add.
  if (length(solo_points) > 0) {
    if (attachAttributes) {
      graph = igraph::add.vertices(graph, nv = length(solo_points), name = solo_points, attr = filter(node_attributes, !!as.symbol(icol) %in% solo_points))
    } else {
      graph = igraph::add.vertices(graph, nv = length(solo_points), name = solo_points)
    }
  }
  return(graph)
}

#' Automatically identify family members of degree n
#'
#' This function identifies individuals ndegree-steps away from the proband in the population graph.
#'
#' @param pop_graph Population graph from prepare_graph()
#' @param ndegree Number of steps away from proband to include
#' @param proband_vec Vector of proband ids to create family graphs for. Must be strings.
#' @param fid Column name of proband ids in the output.
#' @param fam_graph_col Column name of family graphs in the output.
#' @param mindist Minimum distance from proband to exclude in the graph (experimental, untested), defaults to 0, passed directly to make_neighborhood_graph.
#' @param mode Type of distance measure in the graph (experimental, untested), defaults to "all", passed directly to make_neighborhood_graph.
#'
#' @returns Tibble with two columns, family ids (fid) and family graphs (fam_graph_col).
#' @export
#'
#' @examples
#' # See Vignettes.
get_family_graphs = function(pop_graph, ndegree, proband_vec, fid = "fid", fam_graph_col = "fam_graph", mindist = 0, mode = "all") {
  ## TODO: mindist > 0, will get_covmat still work?
  
  # later computations with igraph seem to internally convert ids to strings
  # we will enforce string format for IDs here.
  if (typeof(proband_vec) != "character") {
    warning("get_family_graphs: proband_vec is not a character vector. Converting to character.")
    proband_vec = as.character(proband_vec)
  }
  
  
  tibble(
    !!as.symbol(fid) := proband_vec,
    !!as.symbol(fam_graph_col) := igraph::make_neighborhood_graph(pop_graph, order = ndegree, mindist = mindist, nodes = proband_vec, mode = mode)
  )
}



#' Attach attributes to a family graphs
#'
#' This function attaches attributes to family graphs, such as lower and upper thresholds, for each family member. This allows for a user-friendly way to attach personalised thresholds and other per-family specific attributes to the family graphs.
#'
#' @param cur_fam_graph An igraph object (neighbourhood graph around a proband) with family members up to degree n.
#' @param cur_proband Current proband id (center of the neighbourhood graph).
#' @param pid Column name of personal id (within a family).
#' @param attr_tbl Tibble with family id and attributes for each family member.
#' @param attr_names Names of attributes to be assigned to each node (family member) in the graph.
#' @param proband_cols_to_censor Which columns should be made uninformative for the proband? Defaults to NA. Used to exclude proband's information for prediction with, e.g. c("lower", "upper").
#'
#' @returns igraph object (neighbourhood graph around a proband) with updated attributes for each node in the graph.
#'
#' @export

attach_attributes = function(cur_fam_graph, cur_proband, pid, attr_tbl, attr_names, proband_cols_to_censor = NA) {
  # get node names
  graph_vertex_names = igraph::vertex_attr(cur_fam_graph)$name
  
  # which nodes are present in thresholds?
  to_keep_indx = which(graph_vertex_names %in% attr_tbl[[pid]])
  # get names of present nodes
  to_keep = graph_vertex_names[to_keep_indx]
  
  # order attributes after graph
  attr_tbl_matched = attr_tbl %>% slice(match(to_keep, !!as.symbol(pid)))
  
  # any attr_names columns not in attr_tbl?
  if (any(!(attr_names %in% colnames(attr_tbl_matched)))) {
    warning(paste0("attach_attributes: Not all attributes are present in attr_tbl!\n",
                   " Missing attributes: ", paste(setdiff(attr_names, colnames(attr_tbl_matched)), collapse = ", ")))
  }
  
  # only include columns of attr_tbl that are in attr_names
  for (attr in attr_names) {
    cur_fam_graph <- igraph::set_vertex_attr(graph = cur_fam_graph, name = attr, value = attr_tbl_matched[[attr]])
  }
  
  # censor proband thresholds if requested
  if (is.character(proband_cols_to_censor) & any(!is.na(proband_cols_to_censor))) {
    # are any cols to censor not present?
    if (any(!(proband_cols_to_censor %in% colnames(attr_tbl)))) {
      warning(paste0("attach_attributes: Not all proband_cols_to_censor are present in attr_tbl!\n",
                     " Missing columns: ", paste(setdiff(proband_cols_to_censor, colnames(attr_tbl)), collapse = ", ")))
    }
    
    for (attr in proband_cols_to_censor) {
      cur_fam_graph = igraph::set_vertex_attr(
        graph = cur_fam_graph,
        index = cur_proband,
        name = attr,
        value = case_when(
          str_detect(attr, "lower") ~ -Inf,
          str_detect(attr, "upper") ~ Inf,
          TRUE ~ NA
        )
      )
    }
  }
  
  return(cur_fam_graph)
}


#' Wrapper to attach attributes to family graphs
#'
#' This function can attach attributes to family graphs, such as lower and upper thresholds, for each family member. This allows for personalised thresholds and other per-family specific attributes.
#' This function wraps around attach_attributes to ease the process of attaching attributes to family graphs in the standard format.
#'
#' @param family_graphs tibble with family ids and family graphs
#' @param fam_attr tibble with attributes for each family member
#' @param fam_graph_col column name of family graphs in family_graphs. defailts to "fam_graph"
#' @param attached_fam_graph_col column name of the updated family graphs with attached attributes. defaults to "masked_fam_graph".
#' @param fid column name of family id. Typically contains the name of the proband that a family graph is centred on. defaults to "fid".
#' @param pid personal identifier for each individual in a family. Allows for multiple instances of the same individual across families. Defaults to "pid".
#' @param cols_to_attach columns to attach to the family graphs from fam_attr, typically lower and upper thresholds. Mixture input also requires K_i and K_pop.
#' @param proband_cols_to_censor Should proband's upper and lower thresholds be made uninformative? Defaults to TRUE. Used to exclude proband's information for prediction.
#'
#' @returns tibble with family ids and an updated family graph with attached attributes. If lower and upper thresholds are specified, the input is ready for estimate_liability().
#'
#' @export
#'
#' @examples
#' # See Vignettes.
familywise_attach_attributes = function(family_graphs,
                                        fam_attr,
                                        fam_graph_col = "fam_graph",
                                        attached_fam_graph_col = "masked_fam_graph",
                                        fid = "fid",
                                        pid = "pid",
                                        cols_to_attach = c("lower", "upper"),
                                        proband_cols_to_censor = NA) {
  fam_attr %>%
    tidyr::nest(attrs = -fid) %>%
    left_join(family_graphs, ., by = fid) %>%
    mutate(
      !!as.symbol(attached_fam_graph_col) := purrr::pmap(
        .l = list(!!as.symbol(fid), !!as.symbol(fam_graph_col), attrs),
        ~ attach_attributes(
          cur_fam_graph = ..2,
          attr_tbl = ..3, cur_proband = ..1,
          pid = pid, attr_names = cols_to_attach,
          proband_cols_to_censor  = proband_cols_to_censor))) %>%
    select(-attrs, -!!as.symbol(fam_graph_col))
}