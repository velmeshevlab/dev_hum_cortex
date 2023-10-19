function (x, listw, wc, alternative = "greater", randomisation = TRUE)
{
    zero.policy = TRUE
    adjust.n = TRUE
    na.action = stats::na.fail
    drop.EI2 = FALSE
    xname <- deparse(substitute(x))
    wname <- deparse(substitute(listw))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    n <- length(listw$neighbours)
    S02 <- wc$S0 * wc$S0
    model <- lm(exp ~ sex + region + as.numeric(PMI) + as.numeric(UMI) + chemistry, data = x)
    res <- spdep::lm.morantest(model, listw, zero.policy = zero.policy, alternative = alternative)
    statistic = as.numeric(res[1])
    names(statistic) <- "Moran I statistic standard deviate"
    PrI = as.numeric(res[2])
    vec <- c(res[3]$estimate[1], res[3]$estimate[2], res[3]$estimate[3])
    names(vec) <- c("Moran I statistic", "Expectation", "Variance")
    method <- paste("Moran I test under", ifelse(randomisation, 
        "randomisation", "normality"))
    res <- list(statistic = statistic, p.value = PrI, estimate = vec)
    if (!is.null(na.act)) 
        attr(res, "na.action") <- na.act
    class(res) <- "htest"
    res
}
