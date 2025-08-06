#' Regression on Phase (RoP) test
#'
#' Performs the Regression-on-Phase test to determine whether the \emph{cis}
#' or \emph{trans} configuration of two genetic variants influences the
#' outcome. Supports both biallelic and multiallelic loci and handles
#' continuous (Gaussian) or binary (binomial) responses.
#'
#' @param Y Numeric vector of length n containing the outcome.
#'          For \code{family = "gaussian"} values are treated as continuous;
#'          for \code{family = "binomial"} they must be coded as \code{0/1}.
#' @param Hap1 Numeric matrix of dimension n by 2 giving the haplotype on
#'             one homologous chromosome at the two variants.  Entries are
#'             integers indexing the observed alleles; within each column the
#'             smallest integer is taken as the reference allele.
#' @param Hap2 Numeric matrix of dimension n by 2 giving the haplotype on
#'             the other homologous chromosome, structured identically to
#'             \code{Hap1}.
#' @param family Character string specifying the GLM family to use.
#'               Accepts \code{"gaussian"} for continuous outcomes
#'               or \code{"binomial"} for binary outcomes.
#'
#' @return A named vector with elements \code{p_cis} and \code{p_trans}
#'         containing the p-values for the cis- and trans-phase effects.
#' @importFrom stats glm
#' @importFrom lmtest lrtest
#' @examples
#' ##Example1: Continous Y, biallelic variants, Cis effect
#' N=1000
#' Hap1<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
#' Hap2<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
#' Cis<-Hap1[,1]*Hap1[,2]+Hap2[,1]*Hap2[,2]
#' Y<-Cis+rnorm(n=N)
#' RoP(Y,Hap1,Hap2,family="gaussian")
#'
#' ##Example2: Binary Y,biallelic variants, Trans effect
#' N=1000
#' Hap1<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
#' Hap2<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
#' Trans<-Hap1[,1]*Hap2[,2]+Hap2[,1]*Hap1[,2]
#' p<-exp(2*Trans)/(1+exp(2*Trans))
#' Y<-rbinom(n=1000,size = 1,prob=p)
#' RoP(Y,Hap1,Hap2,family="binomial")
#'
#' ##Example3: Continous Y, multi-allelic variants, Cis effect
#' N=1000
#' Hap1<-cbind(matrix(sample(c(0,1,2),N,prob=c(0.3,0.4,0.3),replace=TRUE),ncol=1),
#'             matrix(sample(c(0,1,2,3),N,prob=c(0.25,0.25,0.25,0.25),replace=TRUE),ncol=1))
#' Hap2<-cbind(matrix(sample(c(0,1,2),N,prob=c(0.3,0.4,0.3),replace=TRUE),ncol=1),
#'             matrix(sample(c(0,1,2,3),N,prob=c(0.25,0.25,0.25,0.25),replace=TRUE),ncol=1))
#'
#' Cis<-as.numeric((Hap1[,1]==1)&(Hap1[,2]==3))+as.numeric((Hap2[,1]==1)&(Hap2[,2]==3))
#' Y<-Cis+rnorm(n=N)
#' RoP(Y,Hap1,Hap2,family="gaussian")
#' @export
RoP<-function(Y, Hap1, Hap2,family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  fam <- get(family, mode = "function")()

  #Check data format
  if (!requireNamespace("lmtest", quietly = TRUE))
    stop("Package 'lmtest' is required.")

  if (!is.matrix(Hap1) || !is.matrix(Hap2) ||
      ncol(Hap1) != 2   || ncol(Hap2) != 2)
    stop("Hap1 and Hap2 must be n by 2 matrices (variant A in col 1, variant B in col 2).")

  if (!all(Hap1 == as.integer(Hap1)) || !all(Hap2 == as.integer(Hap2)))
    stop("`Hap1` and `Hap2` must contain integer allele codes (0, 1, 2,...).")

  n <- nrow(Hap1)
  if (length(Y) != n)
    stop("Length of Y must equal number of rows in Hap1 / Hap2.")

  alleA <- sort(unique(c(Hap1[,1], Hap2[,1])))
  alleB <- sort(unique(c(Hap1[,2], Hap2[,2])))
  if (length(alleA) < 2 || length(alleB) < 2)
    stop("Both variants should have at least one alternative allele.")
  refA <- alleA[1];  altA <- alleA[-1]
  refB <- alleB[1];  altB <- alleB[-1]

  #Create design matrix
  GA <- lapply(altA, function(a) (Hap1[,1] == a) + (Hap2[,1] == a))
  GB <- lapply(altB, function(b) (Hap1[,2] == b) + (Hap2[,2] == b))
  DA <- lapply(GA,  function(g) as.integer(g == 1))
  DB <- lapply(GB,  function(g) as.integer(g == 1))
  names(GA) <- paste0("GA_", altA)
  names(GB) <- paste0("GB_", altB)
  names(DA) <- paste0("DA_", altA)
  names(DB) <- paste0("DB_", altB)

  Cis   <- list();  Trans <- list()
  for (a in altA){
    for (b in altB) {
    Cis [[paste0("Cis_",   a, "_", b)]] <-
      (Hap1[,1] == a & Hap1[,2] == b) + (Hap2[,1] == a & Hap2[,2] == b)

    Trans[[paste0("Trans_", a, "_", b)]] <-
      (Hap1[,1] == a & Hap2[,2] == b) + (Hap2[,1] == a & Hap1[,2] == b)
    }
  }

  Xfull<- data.frame(GA, GB, DA, DB, Cis, Trans)

  #Testing
  fit_full   <- glm(Y ~ ., data = Xfull, family = fam)
  fit_noCis  <- glm(Y ~ ., data = Xfull[ , !(names(Xfull) %in% names(Cis)),   drop = FALSE],
                    family = fam)
  fit_noTrans<- glm(Y ~ ., data = Xfull[ , !(names(Xfull) %in% names(Trans)), drop = FALSE],
                    family = fam)

  p_cis<-lrtest(fit_noCis,fit_full)$`Pr(>Chisq)`[2]
  p_trans<-lrtest(fit_noTrans,fit_full)$`Pr(>Chisq)`[2]

  c(p_cis = p_cis, p_trans = p_trans)

}


