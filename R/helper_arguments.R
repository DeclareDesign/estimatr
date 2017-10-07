# Turn unquoted variables into strings
deparse_var <- function(var_name) {
  if (is.character(var_name)) {
    return(var_name)
  } else {
    return(deparse(var_name))
  }
}
