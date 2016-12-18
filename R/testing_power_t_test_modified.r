#-------------------------------------------------------------------------
#Test om power.t.test.mod giver samme resultat som power.t.test n?r rho=0
#Der laves 100 test, hvor parametrene for hver test genereres tilf?ldigt indenfor visse betingelser

test_against_power.t.test <- function() {
  for (i in 1:100) {
    this_n <- sample(50:150,1)
    num_clus_1 <- sample(3:6,1)
    num_clus_2 <- sample(3:6,1)
    k_1 <- floor(this_n/num_clus_1)
    k_2 <- floor(this_n/num_clus_2)
    clus_size_1 <- sample(5:k_1,num_clus_1-1,replace=TRUE)
    last_one_1 <- (this_n - sum(clus_size_1))
    if (last_one_1 <= 0) {stop("Wrong logic, stupid! (1)")}
    clus_size_1 <- c(clus_size_1,last_one_1)
    clus_size_2 <- sample(5:k_2,num_clus_2-1,replace=TRUE)
    last_one_2 <- (this_n - sum(clus_size_2))
    if (last_one_2 <= 0) {stop("Wrong logic, stupid! (2)")}
    clus_size_2 <- c(clus_size_2,last_one_2)
    if ((!(sum(clus_size_1)==sum(clus_size_2))) || (!(sum(clus_size_1)==this_n))) {
      stop("Wrong logic, stupid! (3)")
    }
    this_delta <- runif(1,min=-3,max=3)
    this_sd <- runif(1,min=1,10)
    this_sig <- runif(1, min=0.01,max=0.1)
    find_alt <- sample(1:2,1)
    if (find_alt==1) {
      this_alt<-"one.sided"
      this_strict<-FALSE
    }
    if (find_alt==2) {
      this_alt<-"two.sided"
      find_strict <- sample(0:1,1)
      if (find_strict==0) { 
        this_strict<-FALSE
      }
      if (find_strict==1) {
        this_strict<-TRUE
      }
    }
    pow1 <- (power.t.test(n=this_n,delta=this_delta,sd=this_sd,sig.level=this_sig,power=NULL,type="two.sample",alternative=this_alt,strict=this_strict))$power
    pow2 <- power.t.test.mod(clus_size_1,clus_size_2,delta=this_delta,sd=this_sd,rho=0,sig.level=this_sig,alternative=this_alt,strict=this_strict)
    print("POWER from original:")
    print(pow1)
    print("POWER from modified:")
    print(pow2)
    if (!isTRUE(all.equal(pow1,pow2))) {
      print("n")
      print(this_n)
      print("clus_size_1")
      print(clus_size_1)
      print("clus_size_2")
      print(clus_size_2)
      print("this_delta")
      print(this_delta)
      print("this_sd")
      print(this_sd)
      print("this_sig")
      print(this_sig)
      print("this_alt")
      print(this_alt)
      stop("Something is HORRIBLY wrong!")
    }
    rm(clus_size_1,clus_size_2,find_alt,k_1,k_2,last_one_1,last_one_2,num_clus_1,num_clus_2,pow1,pow2,this_alt,this_delta,this_n,this_sd,this_sig)
  }
  print("HOOOOORAY!")
}
#Testen er succesfuld.
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#Power.t.test.mod testes mod simulerede data

#..................................................
#Vi starter med at lave funktioner, der kan simulere data med den ønskede varians-struktur
data_generate_group <- function(num_obs_clus, group_mean, sd_group, rho) {
  var_group <- sd_group^2
  var_among_clus <- rho*var_group #udregner variansen among clusters
  var_within_clus <- var_group - var_among_clus #udregner variansen within clusters
  clus_effect <- rnorm(n=length(num_obs_clus),mean=0,sd=sqrt(var_among_clus)) #finder cluster-effekten for hvert cluster vha. among-clus-var
  data<-vector(mode='numeric',length = 0)
  for (i in 1:length(num_obs_clus)) {
    x <- ( (rnorm(num_obs_clus[i], mean=group_mean, sd=sqrt(var_within_clus))) + clus_effect[i]) #Genererer normalfordelte m?linger med within-clus-var st?j og cluster-effekt
    data <- c(data,x)
  }
  return (data)
}
#Genererer data for en enkelt gruppe (kontrol eller behandling)
#Dataen genereres ud fra en mixed model med fast effekt for hver gruppe og tilf?ldig effekt for hvert cluster
#Parametre:
#num_obs_clus: vektor med antallet af m?linger i hvert cluster i gruppen
#group_mean: gruppens grand mean
#sd_group: den samlede standardafvigelse for gruppens m?linger
#rho: intra-cluster korrelations-koefficienten for gruppens clustre

data_generate_experiment <- function(num_obs_clus_group_1, num_obs_clus_group_2=num_obs_clus_group_1, 
                                     delta, sd = 1, rho) {
  respons_group_1 <- data_generate_group(num_obs_clus_group_1,group_mean=0,sd_group=sd,rho)
  respons_group_2 <- data_generate_group(num_obs_clus_group_2,group_mean=delta,sd_group=sd,rho)
  respons <- c(respons_group_1,respons_group_2)
  group <- as.factor(c(rep(1,sum(num_obs_clus_group_1)),rep(2,sum(num_obs_clus_group_2))))
  num_obs_clusters <- c(num_obs_clus_group_1,num_obs_clus_group_2)
  cluster <- vector(mode="numeric",length = 0)
  for (i in 1:length(num_obs_clusters)) {
    cluster <- c(cluster,rep(i,num_obs_clusters[i]))
  }
  cluster <- as.factor(cluster)
  data <- data.frame(group,cluster,respons)
  return (data)
}

#Generer data for kontrolgruppe og behandlingsgruppe
#Parametre:
#num_obs_clus_group_i: vektor med antallet af m?linger i hvert cluster i gruppe i
#delta: behandlingseffekten - alts? forskellen p? gruppernes grand means
#sd: den samlede standardafvigelse for hver af gruppernes m?linger
#rho: intra-cluster korrelations-koefficienten for gruppernes clustre
test_data_generate_experiment <- function(num_test,num_experiments) {
  test_results <- matrix(NA,nrow=num_test,ncol=2)
  colnames(test_results) <- c("sd-est-err.", "rho-est-err.")
  for (j in 1:num_test) {
    
    num_clus <- sample(5:10,1)
    num_obs_clus_group_1 <- sample(50:75,num_clus)
    num_obs_clus_group_2 <- sample(50:75,num_clus)
    delta <- runif(1,min=0,max=10)
    sd <- runif(1,min=1,max=10)
    rho <- runif(1,min=0,max=1)
    
    results <- matrix(data=NA,nrow=num_experiments,ncol=2)
    
    for (i in 1:num_experiments) {
      data_experiment <- data_generate_experiment(num_obs_clus_group_1,num_obs_clus_group_2, delta, sd, rho)
      mixed_model <- lmer(respons~group+(1|cluster),data=data_experiment)
      var_within_est <- (attr(VarCorr(mixed_model),"sc"))^2
      var_among_est <- (as.numeric(attr(VarCorr(mixed_model)$cluster,"stddev")))^2
      rho_est <- var_among_est/(var_among_est+var_within_est)
      sd_est <- sqrt(var_among_est + var_within_est)
      results[i,1] <- sd_est
      results[i,2] <- rho_est
    }
    
    mean_sd_est <- mean(results[,1])
    mean_rho_est <- mean(results[,2])
    
    test_results[j,1] <- (mean_sd_est - sd)
    test_results[j,2] <- (mean_rho_est - rho)
  }
  return (test_results)
}
#Tester om data_generate_experiment fungerer efter hensigten ved at se på, hvor godt parametrene estimeres af lmer
#...............................................................

#Herefter laver vi den cluster-sampling-korrigerede t-test, jf Donner og Klark side 112-115
#Vi bruger lmer fra biblioteket lme4 til at få estimater for variansen mellem clustrer og indenfor clustre.
#Ud fra disse kan vi beregne et estimat for intracluster correlations koefficienten rho.

lmer_rho_est <- function(data_experiment) {
  mixed_model <- lmer(respons~group+(1|cluster),data=data_experiment)
  var_within_est <- (attr(VarCorr(mixed_model),"sc"))^2
  var_among_est <- (as.numeric(attr(VarCorr(mixed_model)$cluster,"stddev")))^2
  var_total_est <- var_within_est + var_among_est
  rho_est <- var_among_est/(var_among_est+var_within_est)
  return(rho_est)
}

m_Ai_from_data <- function(data_experiment_group_i) {
  clus_sizes <- as.numeric(table(data_experiment_group_i$cluster))
  group_size <- length(data_experiment_group_i$respons)
  return ((sum(clus_sizes^2))/group_size)
}

var_IF_group_i_est <- function(data_experiment_group_i,rho_est) {
  this_m_Ai_from_data <- m_Ai_from_data(data_experiment_group_i)
  return (1+(this_m_Ai_from_data-1)*rho_est)
}

pooled_sd <- function(data_experiment) {
  group_1_respons <- subset(data_experiment,group==1)$respons
  group_2_respons <- subset(data_experiment,group==2)$respons
  size_group_1 <- length(group_1_respons)
  size_group_2 <- length(group_2_respons)
  pooled_sd <- ((size_group_1-1)*sd(group_1_respons)+(size_group_2-1)*sd(group_2_respons))/(size_group_1+size_group_2-2)
  return (pooled_sd)
}

SE_est_fact <- function(data_experiment) {
  
  data_group_1 <- subset(data_experiment,group==1)
  data_group_2 <- subset(data_experiment,group==2)
  size_group_1 <- length(data_group_1$respons)
  size_group_2 <- length(data_group_2$respons)
  
  rho_est <- lmer_rho_est(data_experiment)
  
  var_IF_group_1 <- var_IF_group_i_est(data_group_1,rho_est)
  var_IF_group_2 <- var_IF_group_i_est(data_group_2,rho_est)
  
  a <- var_IF_group_1/size_group_1
  b <- var_IF_group_2/size_group_2
  
  return (sqrt(a+b))
  
}

SE_est <- function(data_experiment) {
  S_p <- pooled_sd(data_experiment)
  SE_fact <- SE_est_fact(data_experiment)
  return (S_p*SE_fact) 
}

t_test_stat <- function(data_experiment) {
  group_1_mean <- mean(subset(data_experiment,group==1)$respons)
  group_2_mean <- mean(subset(data_experiment,group==2)$respons)
  this_SE_est <- SE_est(data_experiment)
  return ((group_2_mean-group_1_mean)/this_SE_est)
}

t_test_pval <- function(data_experiment, alt, strict) {
  k <- nlevels(data_experiment$cluster)
  t <- t_test_stat(data_experiment)
  t <- abs(t)
  if (alt == "two.sided" && strict) {
    p <- pt(t,df=k-2,lower.tail = FALSE) + pt(-t,df=k-2,lower.tail=TRUE)
  }
  else {
    p <- pt(t,df=k-2,lower.tail = FALSE)
  }
  return (p)
}

t_test_mod <- function(data_experiment, alt, strict) {
  p <- t_test_pval(data_experiment, alt, strict)
  if (p<0.05) {return (0)}
  else {return (1)}
}
#.........................................................................................

#Vi tester nu styrken beregnet af power.t.test.mod mod hvor mange procent af de simulerede eksperimenter, som 
#den korrigerede t-test forkaster nul-hypotesen for

test_t_test_mod <- function(num_test=100, arg_num_obs_clus_g1=NULL, arg_num_obs_clus_g2=NULL, arg_delta=NULL, 
                            arg_sd=NULL, arg_rho=NULL, arg_alt=NULL, arg_strict=NULL) {
  
  accepted <- 0
  
  for (j in 1:num_test) {
    
    num_clus <- sample(5:10,1)
    
    if (is.null(arg_num_obs_clus_g1)) {num_obs_clus_group_1 <- sample(50:75,num_clus)}
    else {num_obs_clus_group_1 <- arg_num_obs_clus_g1}
    if (is.null(arg_num_obs_clus_g1)) {num_obs_clus_group_2 <- sample(50:75,num_clus)}
    else {num_obs_clus_group_2 <- arg_num_obs_clus_g2}
    if (is.null(arg_sd)) {sd <- runif(1,min=0,max=10)}
    else {sd <- arg_sd}
    if (is.null(arg_rho)) {rho <- runif(1,min=0,max=1)}
    else {rho <- arg_rho}
    if (is.null(arg_alt)) {
      s <- sample(1:2,1)
      if (s==1) {alt <- "one.sided"}
      if (s==2) {alt <- "two.sided"}
    }
    else {alt <- arg_alt}
    if (is.null(arg_delta)) {
      if (alt == "one.sided") {delta <- runif(1,min=0,max=5)}
      if (alt == "two.sided") {delta <- runif(1,min=-5,max=5)}
    }
    else {delta <- arg_delta}
    if (is.null(arg_strict)) {
      s <- sample(0:1,1)
      if (s==0) {strict <- FALSE}
      if (s==1) {strict <- TRUE}
    }
    else {strict <- arg_strict}
    
    test_data <- data_generate_experiment(num_obs_clus_group_1,num_obs_clus_group_2,delta,sd,rho)
    
    accepted <- accepted + t_test_mod(test_data,alt,strict)
  }
  
  return (accepted/num_test)
  
}

test_against_sim_experiments <- function(num_test_1=10, num_test_2=100, arg_num_obs_clus_g1=NULL,
                                         arg_num_obs_clus_g2=NULL,arg_delta=NULL,arg_sd=NULL,arg_rho=NULL,
                                         arg_alt = NULL, arg_strict=NULL) {
  
  results <- matrix(NA,nrow = num_test_1+1,ncol = 3)
  colnames(results) <- c("pow-cal","percent-null-denied","difference")
  
  for (j in 1:num_test_1) {
    
    num_clus <- sample(5:10,1)
    
    if (is.null(arg_num_obs_clus_g1)) {num_obs_clus_group_1 <- sample(50:75,num_clus)}
    else {num_obs_clus_group_1 <- arg_num_obs_clus_g1}
    if (is.null(arg_num_obs_clus_g1)) {num_obs_clus_group_2 <- sample(50:75,num_clus)}
    else {num_obs_clus_group_2 <- arg_num_obs_clus_g2}
    if (is.null(arg_sd)) {sd <- runif(1,min=0,max=10)}
    else {sd <- arg_sd}
    if (is.null(arg_rho)) {rho <- runif(1,min=0,max=1)}
    else {rho <- arg_rho}
    if (is.null(arg_alt)) {
      s <- sample(1:2,1)
      if (s==1) {alt <- "one.sided"}
      if (s==2) {alt <- "two.sided"}
    }
    else {alt <- arg_alt}
    if (is.null(arg_delta)) {
      if (alt == "one.sided") {delta <- runif(1,min=0,max=5)}
      if (alt == "two.sided") {delta <- runif(1,min=-5,max=5)}
    }
    else {delta <- arg_delta}
    if (is.null(arg_strict)) {
      s <- sample(0:1,1)
      if (s==0) {strict <- FALSE}
      if (s==1) {strict <- TRUE}
    }
    else {strict <- arg_strict}
    
    perc_accept <- test_t_test_mod(num_test = num_test_2, arg_num_obs_clus_g1 = num_obs_clus_group_1,
                                   arg_num_obs_clus_g2 = num_obs_clus_group_2, arg_delta = delta, 
                                   arg_sd = sd, arg_rho = rho, arg_alt = alt, arg_strict = strict)
    perc_denied <- 1 - perc_accept
    
    pow <- power.t.test.mod(cluster_sizes_group_1=num_obs_clus_group_1,
                            cluster_sizes_group_2=num_obs_clus_group_2,
                            delta=delta,sd=sd,rho=rho,alternative = alt,
                            strict=strict)$power
    
    results[j,1] <- 100*pow
    results[j,2] <- 100*perc_denied
    results[j,3] <- 100*(pow - perc_denied)
    
  }
  
  mean_diff <- mean(results[,3],na.rm=TRUE)
  results[num_test_1+1,3] <- mean_diff
  
  return (results)
  
}

#-------------------------------------------------------------------------

