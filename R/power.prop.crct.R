#' Power Calculations for Clustered Two-Sample Test for Proportions
#'
#' Compute power of test, or determine parameters to obtain target power.
#'
#' The procedure uses uniroot to find the root of a discontinuous function so
#' some errors may pop up due to the given setup that causes the root-finding
#' procedure to fail. Also, since exact binomial tests are used we have
#' discontinuities in the function that we use to find the root of but despite
#' this the function is usually quite stable.
#'
#' @param n Number of observations
#' @param p0 Probability under the null
#' @param pa Probability under the alternative
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param alternative One- or two-sided test
#' @return Object of class \code{power.htest}, a list of the arguments
#' (including the computed one) augmented with method and note elements.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{binom.test}}
#' @keywords htest
#' @examples
#'
#'
#'
#' @export power.prop.crct
power.prop.crct <- function(n = NULL,
                            p1 = NULL,
                            p2 = NULL,
                            sig.level = 0.05,
                            power = NULL,
                            rho=NULL,
                            nclusters=NULL,
                            alternative = c("two.sided", "one.sided"), strict = FALSE,
                            tol = .Machine$double.eps^0.25) {



    power <- function(p2, p1=0.018, alpha=0.05, m=50000/2, rho=.1, n=50000/20) {




    p.body <- quote({
        DE <- mean(cs[group1])*length(group1) / (sum(cs[group1]/(1 + (cs[group1]-1)*rho)))
        pnorm(  sqrt( m/(DE) * (p1-p2)^2 / ((p1*(1-p1) + p2*(1-p2)))) - qnorm(1-alpha/2) )

        pnorm(  sqrt(
            m/(1 + (n-1)*rho) * (p1-p2)^2 / ((p1*(1-p1) + p2*(1-p2)))
            ) - qnorm(1-alpha/2) )
        #  (qnorm(1-alpha/2) + qnorm(power))*(*(1+(n-1)*rho)/(p1-p2)^2
    }


    })
    if (is.null(power))
        power <- eval(p.body)




    if (sum(sapply(list(n, p0, pa, power, sig.level), is.null)) != 1)
        stop("exactly one of 'n', 'p0', 'pa', 'power', and 'sig.level' must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >
        sig.level | sig.level > 1))
        stop("'sig.level' must be numeric in [0, 1]")

    alternative <- match.arg(alternative)


    NOTE <- NULL
    METHOD <- "One-sample exact binomial power calculation"
    structure(list(n = n, p0 = p0, pa = pa, sig.level = sig.level,
        power = power, alternative = alternative, note = NOTE,
        method = METHOD), class = "power.htest")
}


