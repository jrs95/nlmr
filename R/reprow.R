# Repeat rows
#
# This function creates a matrix of a repeated vector where each row is the same
# @param x vector to be repeated
# @param n number of repeats
# @export
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
