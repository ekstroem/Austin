#' Power Calculations for Exact and Asymptotic McNemar Test in a 2 by 2 table
#'
#' Compute power of test, or determine parameters to obtain target power for
#' matched case-control studies.
#'
#' If psi is less than 1 then the two probabilities p_12 and p_21 are reversed.
#'
#' @param n Number of observations (number of pairs)
#' @param paid The probability that a case patient is not exposed and that the
#' corresponding control patient was exposed (specifying p_12 in the 2 x 2 table).
#' @param psi The relative probability that a control patient is not exposed and that the corresponding case patient was exposed compared to the probability that a case patient is not exposed and that the corresponding control patient was exposed (p12 / p21 in the 2x2 table). Also called the discordant proportion ratio
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param alternative One- or two-sided test
#' @param method Power calculations based on exact or asymptotic test. The default (normal) corresponds to an approximative test, "exact" is the unconditional exact test, while "cond.exact" is a conditional exact test (given fixed n).
#' @return Object of class \code{power.htest}, a list of the arguments
#' (including the computed one) augmented with method and note elements.
#' @note \code{uniroot} is used to solve power equation for unknowns, so you
#' may see errors from it, notably about inability to bracket the root when
#' invalid arguments are given.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{mcnemar.test}}
#' @references Duffy, S (1984). Asymptotic and Exact Power for the McNemar Test
#' and its Analogue with R Controls per Case
#'
#' Fagerland MW, Lydersen S, Laake P. (2013) The McNemar test for binary matched-pairs data: mid-p and asymptotic
#' are better than exact conditional. BMC Medical Research Methodology.
#' @keywords htest
#' @examples
#'
#' # Assume that pi_21 is 0.125 and we wish to detect an OR of 2.
#' # This implies that pi_12=0.25, and with alpha=0.05, and a power of 90% you get
#' power_mcnemar_test(n=NULL, paid=.125, psi=2, power=.9)
#' 
#' power_mcnemar_test(n=NULL, paid=.1, psi=2, power=.8, method="normal")
#' power_mcnemar_test(n=NULL, paid=.1, psi=2, power=.8)
#'
#'
#'
#' @export
power_mcnemar_test <- function(n = NULL, paid = NULL, psi = NULL, sig.level = 0.05, power = NULL,
                               alternative = c("two.sided", "one.sided"),
                               method = c("normal", "exact", "cond.exact")) {

    if (sum(sapply(list(n, paid, psi, power, sig.level), is.null)) != 1)
        stop("exactly one of 'n', 'paid', 'psi', 'power', and 'sig.level' must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
        stop("'sig.level' must be numeric in [0, 1]")
    if (any(paid <=0) || any(paid>=1)) {
        stop("paid is a probability and must be 0<paid<1")
    }
    if (any(psi <=0)) {
        stop("psi must be non-negative")
    }    
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    
    ## Fix if psi was specified to be less that 1
    if (psi<1) {
        paid <- paid*psi
        psi <- 1/psi
    }

    ## Conditional power (conditional on n)
    f <- function(n, paid, psi, sig.level, power) {
        bc <- ceiling(paid * n * (1+psi))
        pbinom(qbinom(0.025, size=bc, prob=0.5)-1, size=bc, prob=1/(1+psi)) + 1-pbinom(qbinom(0.975, size=bc, prob=0.5), size=bc, prob=1/(1+psi))
    }

    ## Unconditional power
    fexact <- function(n, paid, psi, sig.level, power, alt=alternative) {
        sum(dbinom(seq(n), size=n, prob=paid*(1+psi))*power_binom_test(seq(n), p0=.5, pa=1/(1+psi),
                                                                       power=power, sig.level=sig.level,
                                                                       alternative=ifelse(alt=="two.sided", "two.sided", "less"))$power)
    }

        
    if (method=="normal") {
        p.body <- quote( pnorm (
            (sqrt(n * paid) * (psi-1) - qnorm(sig.level/tside, lower.tail=FALSE)*sqrt(psi+1)) / sqrt((psi+1) - paid*(psi-1)^2)))
    } else if (method=="exact") { 
        p.body <- quote( fexact(n, paid, psi, sig.level, power) )
    } else {
        p.body <- quote( f(n, paid, psi, sig.level, power) )
    }


    if (is.null(power)) {
        power <- eval(p.body)
    } else if (is.null(n)) {
        n <- uniroot(function(n) eval(p.body) - power, c(ceiling(log(sig.level)/log(.5)), 1e+07))$root
    } else if (is.null(paid))
        paid <- uniroot(function(paid) eval(p.body) - power, c(0, 1-psi-1e-10))$root
    else if (is.null(psi))
        psi <- uniroot(function(psi) eval(p.body) - power, c(1e-10, 1/paid-1-1e-10))$root
    else if (is.null(sig.level))
        sig.level <- uniroot(function(sig.level) eval(p.body) -
            power, c(1e-10, 1 - 1e-10))$root
    else stop("internal error", domain = NA)
    NOTE <- "n is number of pairs"
    METHOD <- paste("McNemar paired comparison of proportions", ifelse(method=="normal", "approximate", "exact") ,"power calculation")
    structure(list(n = n, paid = paid, psi = psi, sig.level = sig.level,
        power = power, alternative = alternative, note = NOTE,
        method = METHOD), class = "power.htest")

}
