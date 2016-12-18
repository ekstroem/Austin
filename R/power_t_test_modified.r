#WORK IN PROGRESS!
#-------------------------------------------------------------------

#Styrkeberegninger for t-test, hvor der korrigeres for cluster randomized sampling

#Bygget p? funktionen power.t.test og 
#bogen "Allan Donner & Neil Klar: Design and Analysis of 
#Cluster Randomization Trials in Health Research" (Donner & Klark 2000)

#-----------------------------------------------------------------------
#F?lgende funktioner beregner standard-error p? estimatoren af 
#behandlingseffekten, n?r der korrigeres for intracluster korrelations 
#koefficienten. J?vnf?r Donner & Klark 2000, side 111-116.

#Hjælpe-funktion. J?vnf?r Donner & Klark 2000, side 112.
m_Ai <- function(cluster_sizes_group_i) {
  return(sum(cluster_sizes_group_i^2)/sum(cluster_sizes_group_i))
}

#Hj?lpe-funktion. J?vnf?r Donner & Klark 2000, side 115.
var_IF_group_i <- function(m_Ai,rho) {
  return(1+(m_Ai-1)*rho)
}

#Korrigeret standard-error. J?vnf?r Donner & Klark 2000, side 114.
SE_mod <- function(sd,var_IF_group_1,var_IF_group_2,group_1_size,group_2_size) {
  return (sd*sqrt((var_IF_group_1/group_1_size)+(var_IF_group_2/group_2_size)))
}
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#Styrken af t-test med cluster sampling

#Parametre:
#cluster_sizes_group_1: Vektor med antallet af observationer i hver cluster i gruppe 1
#cluster_sizes_group_2: Vektor med antallet af observationer i hver cluster i gruppe 2
#delta: Behandlingseffekt
#sd: Standardafvigelse
#rho: Intracluster korrelations koefficient
#sig.level: Signifikans-niveau
#alternative: Angiver om testen er et- eller to-sidet
#strict: Angiver om der ogs? skal medregnes sandsynlighed i den nedre del af halen i en to-sidet test

power.t.test.mod <- function (cluster_sizes_group_1=NULL, cluster_sizes_group_2=cluster_sizes_group_1, 
                              delta=NULL, sd = NULL, rho = NULL, power = NULL, sig.level = 0.05, 
                              alternative = c("two.sided","one.sided"),strict=FALSE, tol = .Machine$double.eps^0.25,
                              size_calculation=FALSE, num_clusters_pr_group=NULL) {
    
    if (size_calculation) {
      if (is.null(num_clusters_pr_group)) {
        stop("to calculate sample-size, you need to specify the number of clusters per group (parameter:'num_clusters_pr_group')")
      }
    }
    else if (sum(sapply(list(delta, sd, rho, power, sig.level), is.null)) != 1) {
      stop("exactly one of 'delta', 'sd', 'power', and 'sig.level' must be NULL, or size_calculation must be TRUE")
    }

    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    
    if (tside==2 && !is.null(delta)) {
      delta <- abs(delta)
    }
    
    if (size_calculation) {
      p.body <- quote({
        
        var_IF_group_1 <- var_IF_group_i(num_obs_pr_group/num_clusters_pr_group,rho)
        var_IF_group_2 <- var_IF_group_1
        group_1_size <- num_obs_pr_group
        group_2_size <- num_obs_pr_group
        this_SE_mod <- SE_mod(sd, var_IF_group_1, var_IF_group_2, group_1_size, group_2_size)
        this_ncp <- delta/this_SE_mod 
        
        df <- (2*num_clusters_pr_group - 2) 
        
        q <- qt(sig.level/tside, df, lower.tail = FALSE)
        
        if (strict && tside ==2) {
          pt(q, df, ncp = this_ncp, lower.tail = FALSE) + pt(-q, df,ncp=this_ncp,lower.tail=TRUE)
        }
        
        else {
          pt(q, df, ncp = this_ncp, lower.tail = FALSE) 
        }
        
      })
      
    }
    
    else {
      p.body <- quote( {
        
        var_IF_group_1 <- var_IF_group_i(m_Ai(cluster_sizes_group_1),rho)
        var_IF_group_2 <- var_IF_group_i(m_Ai(cluster_sizes_group_2),rho)
        group_1_size <- sum(cluster_sizes_group_1)
        group_2_size <- sum(cluster_sizes_group_2)
        this_SE_mod <- SE_mod(sd,var_IF_group_1,var_IF_group_2,group_1_size,group_2_size)
        this_ncp <- delta/this_SE_mod 
        
        df <- (length(cluster_sizes_group_1) + length(cluster_sizes_group_2) - 2) 

        q <- qt(sig.level/tside, df, lower.tail = FALSE)
        
        if (strict && tside ==2) {
          pt(q, df, ncp = this_ncp, lower.tail = FALSE) + pt(-q, df,ncp=this_ncp,lower.tail=TRUE)
        }
        
        else {
          pt(q, df, ncp = this_ncp, lower.tail = FALSE) 
        }
        
        } )
    }
    
    if (is.null(power)) {power <- eval(p.body)}
    
    else if (is.null(sd)) {
      sd <- uniroot(function(sd) eval(p.body) - power, delta * 
                      c(1e-07, 1e+07), tol = tol, extendInt = "downX")$root
    }
    
    else if (is.null(delta)) {
      delta <- uniroot(function(delta) eval(p.body) - power, 
                       sd * c(1e-07, 1e+07), tol = tol, extendInt = "upX")$root
    }
    
    else if (is.null(sig.level)) {
      sig.level <- uniroot(function(sig.level) eval(p.body) - 
                             power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
    }
    
    else if (is.null(rho)) {
      rho <- uniroot(function(rho) eval(p.body) - 
                          power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
    }
    
    else if (size_calculation) {
        num_obs_pr_group <- tryCatch(uniroot(function(num_obs_pr_group) eval(p.body) - power, 
                          c(1, 1e+10), tol = tol, extendInt = "upX")$root,
                          error = function(e){
                              if (e$message == "did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0") {
                                stop(paste("Error in call to the uniroot function. The error has the message:\n",
                                            "\"", e$message, "\" \n",
                                           "This error might be caused by the experiment needing more than a trillion observations to obtain the specified power\n",
                                           "Try changing the parameters of the power.t.test"))
                              }
                              else {
                                stop(paste("Error in call to the uniroot function. The error has the message:\n",
                                            "\"", e$message, "\" \n"))
                              }
                          }  
                          )
        n <- 2*ceiling(num_obs_pr_group)
    }
    
    else {
      stop("internal error", domain = NA)
    }
    
    METHOD <- paste(switch(tside, one.sample = "One-sample", two.sample = "Two-sample"), 
                    "t test power calculation")
    
    if (size_calculation) {
      NOTE <- "n is the number of observations needed in total"
      structure(list(n = n, clusters_pr_group = num_clusters_pr_group,
                     delta = delta, sd = sd, sig.level = sig.level,
                     power = power, alternative = alternative, rho = rho,
                     note = NOTE, method = METHOD), class = "power.htest")
    }
    
    else {
      structure(list(cluster_sizes_group_1=cluster_sizes_group_1, 
                   cluster_sizes_group_2 = cluster_sizes_group_2,
                   delta = delta, sd = sd, sig.level = sig.level,
                   power = power, alternative = alternative, rho = rho,
                   method = METHOD), class = "power.htest")
    }
      
}

#-----------------------------------------------------------------------




