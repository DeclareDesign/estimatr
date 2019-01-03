

extract_terms <- function(x){
  x <- gsub("*", "", x, fixed=TRUE)
  x <- gsub(" ", "", x, fixed=TRUE)
  x <- gsub("\\n", " ", x)
  x <- gsub("\\t", " ", x)
  x <- gsub("-", "+-", x, fixed = TRUE)
  x <- strsplit(x, "+", fixed = TRUE)[[1]]
  x[x!=""]
}


extract_coefficients <- function( combination, coefnames){
 comb_vector <- extract_terms(combination)
 terms <- gsub("^[-\\ 0-9]+", "", comb_vector)
 i <- coefnames %in%  terms
 if( length(terms) > sum(i)) stop("Linear combination is")
 a <- mapply(function(x, y) unlist(strsplit(y, x, fixed=TRUE)),  terms, comb_vector) #split hypothesis
 a[a == ""] <- 1
 a[a == "-"] <- -1

 out <- coefnames
 names(out) <- out
 out[names(a)] <- a
 out[!i] <- 0
 out <- as.numeric(out)
 names(out) <- coefnames
 out
}

combined_se <- function(l, vcov){
  vcov_lb <-  t(l) %*%  vcov(model) %*% l
  sqrt(diag(vcov_lb))
}

combined_estimate <- function(l, coeffs){
  sum(l*coeffs)
}


combine_coefficients <- function(model, combination, level = 0.95, ...){
 coeffs <- coefficients(model)
 l <- sapply(combination,extract_coefficients, coefnames)
 estimate  <- apply(l, 2, combined_estimate, coeffs)
 std.error <- apply(l, 2, combined_se , coeffs)
 df <- model$df.residual
 tt <- estimate/std.error
 p.value <-  2 * pt(abs(tt), df)
 alpha <- (1 - level)/2
 alpha <- c(alpha, 1 - alpha)
 ci <- estimate + std.error %o% qt(alpha, df)
 conf.low  <-  ci[,1]
 conf.high <-  ci[,2]

 data.frame(estimate = estimate,
            std.error =  std.error,
            p.value = p.value ,
            conf.low = conf.low,
            conf.high = conf.high,
            df = df
            )
}
