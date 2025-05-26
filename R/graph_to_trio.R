utils::globalVariables("M")
utils::globalVariables("F")
utils::globalVariables("sortedID")

#' Convert from igraph to trio information
#' 
#' This function converts an igraph object to a trio information format.
#'
#' @param graph An igraph graph object.
#' @param id Column of proband id. Defaults to id.
#' @param dadid Column of father id. Defaults to dadid.
#' @param momid Column of mother id. Defaults to momid.
#' @param sex Column of sex in igraph attributes. Defaults to sex.
#' @param fixParents Logical. If TRUE, the kinship2's fixParents will be run on the trio information before returning. Defaults to TRUE.
#'
#' @returns A tibble with trio information.
#' 
#' @details The sex column is required in the igraph attributes. The sex information is used to determine who is the mother and father in the trio.
#' 
#' @importFrom igraph E vertex.attributes
#' @importFrom stats setNames
#' @importFrom dplyr select filter left_join rename mutate bind_rows %>%
#' 
#' @export
#' @examples
#' if (FALSE) {
#' 
#' family = tribble(
#' ~id, ~momcol, ~dadcol,
#' "pid", "mom", "dad",
#' "sib", "mom", "dad",
#' "mhs", "mom", "dad2",
#' "phs", "mom2", "dad",
#' "mom", "mgm", "mgf",
#' "dad", "pgm", "pgf",
#' "dad2", "pgm2", "pgf2",
#' "paunt", "pgm", "pgf",
#' "pacousin", "paunt", "pauntH",
#' "hspaunt", "pgm", "newpgf",
#' "hspacousin", "hspaunt", "hspauntH",
#' "puncle", "pgm", "pgf",
#' "pucousin", "puncleW", "puncle",
#' "maunt", "mgm", "mgf",
#' "macousin", "maunt", "mauntH",
#' "hsmuncle", "newmgm", "mgf",
#' "hsmucousin", "hsmuncleW", "hsmuncle"
#' )
#' 
#' 
#' thrs =  tibble(
#'   id = family %>% select(1:3) %>% unlist() %>% unique(),
#'  sex = case_when(
#'    id %in% family$momcol ~ "F",
#'    id %in% family$dadcol ~ "M",
#'     TRUE ~ NA)) %>%
#'   mutate(sex = sapply(sex, function(x) ifelse(is.na(x), 
#'   sample(c("M", "F"), 1), x)))
#' graph = prepare_graph(.tbl = family, 
#' icol = "id", fcol = "dadcol", mcol = "momcol", node_attributes = thrs)
#' graph_to_trio(graph)
#' }
#' 
graph_to_trio = function(graph, id = "id", dadid = "dadid", momid = "momid", sex = "sex", fixParents = TRUE) {
  graph_attr = vertex.attributes(graph) %>% as_tibble() %>% select(name, sex)
  if (typeof(graph_attr[[sex]]) %in% c("integer", "numeric", "double")) {
    graph_attr[[sex]] = c("F", "M")[graph_attr[[sex]] + 1] # assuming sex is coded as 0 for female and 1 for male
    # now, we can work solely with characters
  }
  # get end points of edges
  ph = str_split(attributes(E(graph))$vnames, "\\|")
  
  # split end points into to and from
  ph2 = lapply(setNames(1:2, c("from", "to")), function(i) {
    sapply(ph, function(x) x[i])
  }) %>%
    bind_cols()
  
  # identify any potential sibling links through duplicated (sorted) edges
  ph3 = ph2 %>% mutate(
    sortedID = purrr::map2_chr(.x = from, .y = to,
                               ~ paste(sort(c(.x, .y)), collapse = "_"))
  )
  sib_links = ph3$sortedID[duplicated(ph3$sortedID)]
  
  # remove sibling links (if present) 
  ph4 = ph3 %>%
    filter(!(sortedID %in% sib_links)) %>%
    # attaching sex info
    left_join(graph_attr, by = c("from" = "name"))  %>%
    select(-sortedID) %>%
    # formatting to pedigree format
    tidyr::pivot_wider(names_from = !!as.symbol(sex), values_from = "from") %>%
    # if either M or F does not exist, we will create them here:
    mutate("M" = if(!"M" %in% names(.)) NA else M,
           "F" = if(!"F" %in% names(.)) NA else F) %>%
    rename(!!as.symbol(id) := to,
           !!as.symbol(dadid) := "M",
           !!as.symbol(momid) := "F") %>%
    left_join(graph_attr, setNames(c("name"), id))
  
  
  # Identify individuals without their own row:
  to_add = with(ph4, setdiff(unique(c(dadid, momid)), id))
  
  # Add them to the pedigree with empty columns for dadid and momid
  trio = bind_rows(
    ph4,
    tibble(id = to_add, !!as.symbol(dadid) := "", !!as.symbol(momid) := "") %>%
      left_join(graph_attr, by = setNames(c("name"), id))
  ) %>% 
    filter(!is.na(id)) %>% 
    mutate(
      !!as.symbol(dadid) := ifelse(is.na(!!as.symbol(dadid)), "", !!as.symbol(dadid)),
      !!as.symbol(momid) := ifelse(is.na(!!as.symbol(momid)), "", !!as.symbol(momid))
    ) 
  
  if (fixParents) {
    # fixing parents coding to be suitable for kinship2 input
    if (is.character(trio[[id]])) { # when coded as characters
      # assuming "" is used to indicate missing/unknown values
      # coverting to NAs and performing fixing on NA values.
      trio = trio %>% 
        mutate(!!as.symbol(dadid) := ifelse(!!as.symbol(dadid) == "" & !!as.symbol(momid) != "", NA, !!as.symbol(dadid)),
               !!as.symbol(momid) := ifelse(!!as.symbol(dadid) != "" & !!as.symbol(momid) == "", NA, !!as.symbol(momid)))
      
      missingMaxItr = 1
      for (x in c(momid, dadid)) {
        if ( is.character(trio[[x]]) ) {
          missingIndx = which(is.na(trio[[x]]))
          trio[[x]][missingIndx] <- paste0("added_", missingMaxItr:(missingMaxItr - 1 + length(missingIndx)))
          missingMaxItr = length(missingIndx) + missingMaxItr
          sex_coding = is.character(trio[[sex]])
          
          trio = trio %>% bind_rows(
            tibble(!!as.symbol(id) := trio[[x]][missingIndx],
                   !!as.symbol(sex) := fixSexCoding(x = x, sex_coding = sex_coding, dadid = dadid, momid = momid), 
                   !!as.symbol(dadid) := "",
                   !!as.symbol(momid) := "")
          )
        }
      }
    } 
    
    if (is.integer(trio[id])) { # when coded as integers
      # assuming 0 is used to indicate missing/unknown values
      missingIndxm = which(trio[[momid]] == 0 & trio[[dadid]] != 0)
      max_id = max(trio[, c(id, momid, dadid)], na.rm = T)
      sex_coding = is.character(trio[[sex]])
      
      if (length(missingIndxm) > 0) {
        to_add_m = max_id + 1:length(missingIndxm)
        trio[[momid]][missingIndxm] <- to_add_m
        trio = trio %>% bind_rows(
          tibble(!!as.symbol(id) := to_add_m,
                 !!as.symbol(sex) := fixSexCoding(x = momid, sex_coding = sex_coding, dadid = dadid, momid = momid), 
                 !!as.symbol(dadid) := 0,
                 !!as.symbol(momid) := 0)
        )
      }
      missingIndxf = which(trio[[momid]] != 0 & trio[[dadid]] == 0)
      # now with added mothers
      max_id = max(trio[, c(id, momid, dadid)], na.rm = T)
      if (length(missingIndxf) > 0) {
        to_add_f = max_id + 1:length(missingIndxf)
        trio[[dadid]][missingIndxf] <- to_add_f
        trio = trio %>% bind_rows(
          tibble(!!as.symbol(id) := to_add_f,
                 !!as.symbol(sex) := fixSexCoding(x = dadid, sex_coding = sex_coding, dadid = dadid, momid = momid), 
                 !!as.symbol(dadid) := 0,
                 !!as.symbol(momid) := 0)
        )
      }
      
    }
    
    
  }
  return(trio)
}




#' Fixing sex coding in trio info
#'
#' Internal function used to assist in fixing sex coding separately from id coding type.
#'
#' @param x current row to check against
#' @param sex_coding logical. Is sex coded as character?
#' @param dadid column name of father ids
#' @param momid column name of mother ids
#'
#' @returns appropriate sex coding
#' 
#' @importFrom dplyr case_when
#' 
fixSexCoding = function(x, sex_coding = TRUE, dadid, momid) {
  if (sex_coding) {
    case_when(
      x == dadid ~ "M",
      x == momid ~ "F",
      sex_coding ~ "unknown",
      TRUE ~ NA)
  } else {
    case_when(
      x == dadid ~ 1,
      x == momid ~ 2,
      !sex_coding ~ 3,
      TRUE ~ NA)
  }
}
