combine_coefficients <- function(model, combined_coefficients, level = 0.95, ...){

 out <- car::linearHypothesis(model, combined_coefficients)

 n <- nrow(model$model)
 p <- nrow(out) - 1
 estimate  <- drop(attr(out, "value"))
 std.error <- sqrt(diag(attr(out, "vcov")))
 df <- model$df.residual
 tt <- estimate/std.error
 p.value <-  2 * pt(abs(tt), df)
 alpha <- (1 - level)/2
 alpha <- c(alpha, 1 - alpha)
 ci <- estimate + std.error %o% qt(alpha, df)
 conf.low  <-  ci[,1]
 conf.high <-  ci[,2]
 outcome <- rownames(attr(model$terms, "factors"))[1]

 data.frame(estimator_label = estimator_label,
            term = combined_coefficient,
            estimate = estimate,
            std.error =  std.error,
            p.value = p.value ,
            conf.low = conf.low,
            conf.high = conf.high,
            df = df,
            outcome =  outcome,
            estimand_label = estimand_label
            )
}
