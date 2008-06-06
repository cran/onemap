#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: cr2pts.R                                                      #
# Contains: cr2pts                                                    #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

cr2pts <-
function (mrk1, mrk2, segr.type1, segr.type2) {
  # checking for correct types of segregation
  if (!any(!is.na(match(c("A.1","A.2","A.3","A.4","B1.5","B2.6","B3.7",
                          "C.8","D1.9","D1.10","D1.11","D1.12","D1.13",
                          "D2.14","D2.15","D2.16","D2.17","D2.18"),
                        segr.type1)))) 
    stop("unknown segregation type for 'mrk1'")
  
  if (!any(!is.na(match(c("A.1","A.2","A.3","A.4","B1.5","B2.6","B3.7",
                          "C.8","D1.9","D1.10","D1.11","D1.12","D1.13",
                          "D2.14","D2.15","D2.16","D2.17","D2.18"),
                        segr.type2)))) 
    stop("unknown segregation type for 'mrk2'")

  # the matrixes 'I' and the vectors 'm' are defined for both markers,
  # according to the segregation type
  # 'classes' are vectors indicating possible offspring molecular phenotypes
  switch(EXPR=segr.type1,
         'A.1' = {I1 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m1 <- diag(4);
                  classes1 <- c("ac","ad","bc","bd")},
         'A.2' = {I1 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m1 <- diag(4);
                  classes1 <- c("a","ac","ba","bc")},
         'A.3' = {I1 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m1 <- diag(4);
                  classes1 <- c("ac","a","bc","b")},
         'A.4' = {I1 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m1 <- diag(4);
                  classes1 <- c("ab","a","b","o")},
         'B1.5' = {I1 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m1 <- diag(3);
                   classes1 <- c("a","ab","b")},
         'B2.6' = {I1 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE);
                   m1 <- diag(3);
                   classes1 <- c("a","ab","b")},
         'B3.7' = {I1 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m1 <- diag(3);
                   classes1 <- c("a","ab","b")},
         'C.8' = {I1 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE);
                  m1 <- diag(2);
                  classes1 <- c("a","o")},
         'D1.9' = {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                   m1 <- diag(2);
                   classes1 <- c("ac","bc")},
         'D1.10' = {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("a","ab")},
         'D1.11' = {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("a","b")},
         'D1.12' = {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("ab","a")},
         'D1.13' = {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("a","o")},
         'D2.14' = {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("ac","bc")},
         'D2.15' = {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("a","ab")},
         'D2.16' = {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("a","b")},
         'D2.17' = {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("ab","a")},
         'D2.18' = {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2);
                    classes1 <- c("a","o")})
  
  switch(EXPR=segr.type2,
         'A.1' = {I2 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m2 <- diag(4);
                  classes2 <- c("ac","ad","bc","bd")},
         'A.2' = {I2 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m2 <- diag(4);
                  classes2 <- c("a","ac","ba","bc")},
         'A.3' = {I2 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m2 <- diag(4);
                  classes2 <- c("ac","a","bc","b")},
         'A.4' = {I2 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=TRUE);
                  m2 <- diag(4);
                  classes2 <- c("ab","a","b","o")},
         'B1.5' = {I2 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m2 <- diag(3);
                   classes2 <- c("a","ab","b")},
         'B2.6' = {I2 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE);
                   m2 <- diag(3);
                   classes2 <- c("a","ab","b")},
         'B3.7' = {I2 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m2 <- diag(3);
                   classes2 <- c("a","ab","b")},
         'C.8' = {I2 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE);
                  m2 <- diag(2);
                  classes2 <- c("a","o")},
         'D1.9' = {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                   m2 <- diag(2);
                   classes2 <- c("ac","bc")},
         'D1.10' = {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("a","ab")},
         'D1.11' = {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("a","b")},
         'D1.12' = {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("ab","a")},
         'D1.13' = {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("a","o")},
         'D2.14' = {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("ac","bc")},
         'D2.15' = {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("a","ab")},
         'D2.16' = {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("a","b")},
         'D2.17' = {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("ab","a")},
         'D2.18' = {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2);
                    classes2 <- c("a","o")})
  
  # procedure to count the number of individuals in each class
  n <- numeric(length(classes1)*length(classes2))
  k <- 1
  for(p2 in 1:length(classes2))
    for(p1 in 1:length(classes1)) {
      n[k] <- length(which(mrk1==classes1[p1] & mrk2==classes2[p2]))
      k <- k+1
    }
  ntot <- sum(n) # total number of individuals

  # calling C routine
  rcmb <- .C("r2pts",
             as.double(t(I1)),
             as.integer(length(classes1)),
             as.double(I2),
             as.integer(length(classes2)),
             as.integer(n),
             as.integer(ntot),
             r=as.double(numeric(4)),
             log_like=as.double(numeric(4)),
             posterior=as.double(numeric(4)),
             LOD=as.double(numeric(4)),
             PACKAGE="onemap")
  
  # results
  final <- cbind(rcmb$r,rcmb$log_like,rcmb$posterior,rcmb$LOD)
  dimnames(final) <- list(c("A1","A2","A3","A4"),list("Theta","log-Like","Posterior","LODs"))
  final
}

