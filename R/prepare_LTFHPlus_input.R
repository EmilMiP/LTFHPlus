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
#' @export

prepare_LTFHPlus_input = function(.tbl, 
                                  CIP, 
                                  age_col,
                                  aoo_col,
                                  CIP_merge_columns = c("sex", "birth_year", "age"), 
                                  CIP_cip_col = "cip", 
                                  status_col = "status", 
                                  use_fixed_case_thr = F, 
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
      mutate(age = pmin(!!as.symbol(age_col), !!as.symbol(aoo_col), na.rm = T)) %>% 
      select(all_of(c(CIP_merge_columns, age_col))) %>% 
      as.matrix() %>% 
      predict(xgb,.) %>% 
      pmax(min_CIP_value)
    
    .tbl %>% 
      mutate(event_age = pmin(!!as.symbol(age_col), !!as.symbol(aoo_col), na.rm = T)) %>% 
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
#' @param thresholds tibble with icol, lower_col and upper_col. Used to assign lower and upper thresholds to individuals in the graph as attributes.
#' @param lower_col Column name of column with proband's lower threshold.
#' @param upper_col Column name of column with proband's upper threshold.
#' @param missingID_patterns string of missing values in the ID columns. Multiple values can be used, but must be separated by "|". Defaults to "^0$".
#' 
#' @return An igraph object. A (directed) graph object based on the links provided in .tbl with the lower and upper thresholds stored as attributes. 
#' 
#' @export
#' @importFrom dplyr %>% rename relocate mutate filter group_by summarise select bind_rows pull
#' 
#' @export


prepare_graph = function(.tbl, icol, fcol, mcol, thresholds, lower_col = "lower", upper_col = "upper", missingID_patterns = "^0$") {
 
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
  
  graph = igraph::graph_from_data_frame(d = graph_input, 
                                        #use unique id list to attach threshold info
                                        vertices = filter(thresholds, !!as.symbol(icol) %in% present_ids))
  
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
     graph = igraph::add.vertices(graph, nv = length(solo_points), name = solo_points, attr = filter(thresholds, !!as.symbol(icol) %in% solo_points)) 
  }
  return(graph)
}

