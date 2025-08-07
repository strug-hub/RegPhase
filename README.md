# RegPhase: Regression-on-Phase Test for Genetic Variants

The **RegPhase** package implements the Regerssion on Phpase **(RoP)** test for detecting phase-specific (*cis* or *trans*) effects between **two genetic variants**—either *biallelic* or *multiallelic*—on a continuous or binary outcome. By explicitly modelling *cis* and *trans* configurations, RoP helps determine whether the phase relationship between alleles contributes to trait variation beyond marginal genotype effects.

## Installation

```r
# Install devtools if you don't already have it
install.packages("devtools")

# Install the development version of RoP from GitHub
devtools::install_github("strug-hub/RegPhase")
```

## Function overview
```r
RoP(Y, Hap1, Hap2, family = c("gaussian", "binomial"))
```
| Argument | Description                                                                                                                                              |
| -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Y`      | Numeric vector of length *n* containing the outcome.<br> `family = "gaussian"`: treated as continuous.<br>`family = "binomial"`: treated as binary, coded as 0/1. |
| `Hap1`   | n by 2 integer matrix: haplotype on one homologous chromosome at the two variants.  Entries are integers indexing the observed alleles; within each column the smallest integer is taken as the reference allele.|
| `Hap2`   | n by 2 integer matrix for the other homologous chromosome, structured identically to `Hap1`.                                                            |
| `family` | `"gaussian"` or `"binomial"`, specifying the GLM family.                                                                                       |

**Returns**

A named vector with p-values:
| **Name**      | **Meaning**                                            |
| --------- | ---------------------------------------------------------- |
| `p_cis`   | P-value for the *cis* (same-chromosome) phase effect       |
| `p_trans` | P-value for the *trans* (opposite-chromosome) phase effect |

## Examples
**Example 1:** Continous Y, biallelic variants, Cis effect

```r
N=1000
Hap1<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
Hap2<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
Cis<-Hap1[,1]*Hap1[,2]+Hap2[,1]*Hap2[,2]
Y<-Cis+rnorm(n=N)
RoP(Y,Hap1,Hap2,family="gaussian")
```
**Example 2:** Binary Y,biallelic variants, Trans effect

```r
N=1000
Hap1<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
Hap2<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
Trans<-Hap1[,1]*Hap2[,2]+Hap2[,1]*Hap1[,2]
p<-exp(2*Trans)/(1+exp(2*Trans))
Y<-rbinom(n=1000,size = 1,prob=p)
RoP(Y,Hap1,Hap2,family="binomial")
```
**Example 3:** Continous Y, multi-allelic variants, Cis effect

```r
N=1000
Hap1<-cbind(matrix(sample(c(0,1,2),N,prob=c(0.3,0.4,0.3),replace=TRUE),ncol=1),
            matrix(sample(c(0,1,2,3),N,prob=c(0.25,0.25,0.25,0.25),replace=TRUE),ncol=1))
Hap2<-cbind(matrix(sample(c(0,1,2),N,prob=c(0.3,0.4,0.3),replace=TRUE),ncol=1),
            matrix(sample(c(0,1,2,3),N,prob=c(0.25,0.25,0.25,0.25),replace=TRUE),ncol=1))

Cis<-as.numeric((Hap1[,1]==1)&(Hap1[,2]==3))+as.numeric((Hap2[,1]==1)&(Hap2[,2]==3))
Y<-Cis+rnorm(n=N)
RoP(Y,Hap1,Hap2,family="gaussian")
```
