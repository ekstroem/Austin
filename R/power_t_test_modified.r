#WORK IN PROGRESS!

#Styrkeberegninger for t-test, hvor der korrigeres for cluster randomized sampling

#Bygget på funktionen power.t.test og
#bogen "Allan Donner & Neil Klar: Design and Analysis of
#Cluster Randomization Trials in Health Research" (Donner & Klark 2000)

#-----------------------------------------------------------------------
#Følgende funktioner beregner standard-error på estimatoren af
#behandlingseffekten, når der korrigeres for intracluster korrelations
#koefficienten. Jævnfør Donner & Klark 2000, side 111-116.

#Hjælpe-funktion. Jævnfør Donner & Klark 2000, side 112.
m_Ai <- function(num_obs_clus_group_i) {
  return(sum(num_obs_clus_group_i^2)/sum(num_obs_clus_group_i))
}

#Hjælpe-funktion. Jævnfør Donner & Klark 2000, side 115.
C_i <- function(m_Ai,rho) {
  return(1+(m_Ai-1)*rho)
}

#Korrigeret standard-error. Jævnfør Donner & Klark 2000, side 114.
SE_mod <- function(sd,C_1,C_2,M_1,M_2) {
  return (sd*sqrt((C_1/M_1)+(C_2/M_2)))
}
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#Styrken af t-test med cluster sampling

#Parametre:
#num_obs_clus_group_1: Vektor med antallet af observationer i hver cluster i gruppe 1
#num_obs_clus_group_2: Vektor med antallet af observationer i hver cluster i gruppe 2
#delta: Behandlingseffekt
#sd: Standardafvigelse
#rho: Intracluster korrelations koefficient
#sig.level: Signifikans-niveau
#alternative: Angiver om testen er et- eller to-sidet
#strict: Angiver om der også skal medregnes sandsynlighed i den nedre del af halen i en to-sidet test

power.t.test.mod <-
function (num_obs_clus_group_1, num_obs_clus_group_2, delta, sd = 1, rho, sig.level = 0.05,
          alternative = c("two.sided","one.sided"),strict=FALSE)
{

  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside==2) {delta <- abs(delta)}

  C_1 <- C_i(m_Ai(num_obs_clus_group_1),rho)
  C_2 <- C_i(m_Ai(num_obs_clus_group_2),rho)
  M_1 <- sum(num_obs_clus_group_1)
  M_2 <- sum(num_obs_clus_group_2)
  this_SE_mod <- SE_mod(sd,C_1,C_2,M_1,M_2)
  this_ncp <- delta/this_SE_mod #Non-centralitets-parameter til bestemmelse af fordelingen af t-statistikken under den alternative hypotese
  df <- (M_1 + M_2 - 2)

  q <- qt(sig.level/tside, df, lower.tail = FALSE)

  if (strict && tside==2) {
    power <- (pt(q, df, ncp = this_ncp, lower.tail = FALSE)
              +pt(-q, df,ncp=this_ncp,lower.tail=TRUE))
  }
  else {
    power <- pt(q, df, ncp = this_ncp, lower.tail = FALSE)
  }
  return (power)

}
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#Test om power.t.test.mod giver samme resultat som power.t.test når rho=0
#Der laves 100 test, hvor parametrene for hver test genereres tilfældigt indenfor visse betingelser

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
#---------------------------------------------------------------------------------------------
