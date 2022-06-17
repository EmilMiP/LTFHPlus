#' Check that the strings are valid for family members 
#'
#' \code{check_family_role} checks that the strings used for the relatedness times the
#' heritability for a pair of family members is valid
#'
#' @param s1 String representing the family member.
#' The string must be one of the following:
#' - g (Genetic component of full liability)
#' - o (Full liability)
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - mgp (Maternal grandparent)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - pgp (Paternal grandparent)
#' - s\[0-9\]* (Full siblings)
#' - c\[0-9\]* (children)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#' 
#' @return If s1 is a valid string from the mentioned list of strings no errors will be returned.
#'
#' @examples
#' check_family_role("g")
#' check_family_role("f")
#' check_family_role("s")
#'
#' \dontrun{
#' # This will result in errors:
#' check_family_role("a")
#' check_family_role(mhs)
#' }
#' 
#' @importFrom stringr str_detect
#' 
#' @export
check_family_role = function(s1) {
  # Checking that s1 are valid strings
  sapply(s1, function(s) {
    if (any(!(str_detect(s, "^[gomf]$") | str_detect(s, "^[mp]g[mf]$") | str_detect(s, "^[mp]gp[0-9]*") |
              str_detect(s, "^s[0-9]*") | str_detect(s, "^c[0-9]*") | str_detect(s, "^[mp]hs[0-9]*") | 
              str_detect(s, "^[mp]au[0-9]*")))) stop(paste0(s, " is not a valid string! Use a string from the following list: \n
        - m (Mother)\n
        - f (Father)\n
        - mgm (Maternal grandmother)\n
        - mgf (Maternal grandfather)\n
        - mgp[1-2]* (Maternal grandparent)\n
        - pgm (Paternal grandmother)\n
        - pgf (Paternal grandfather)\n
        - pgp[1-2]* (Paternal grandparent)\n
        - s[0-9]* (Full siblings)\n
        - c[0-9]* (children)\n
        - mhs[0-9]* (Half-siblings - maternal side)\n
        - phs[0-9]* (Half-siblings - paternal side)\n
        - mau[0-9]* (Aunts/Uncles - maternal side)\n
        - pau[0-9]* (Aunts/Uncles - paternal side)."))
  })
  NULL
}

