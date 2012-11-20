##

char.im <- function(x,r=0.2){
  z <- x[1] + 1i * x[2]
  f <- z+r*exp(-z)
  c(real=Re(f),imag=Im(f))
}

library("BB")



BBsolve(fn=char.im,par=c(0.1,0.1),r=0.8)

xs <- seq(0.05,1,by=0.05)

lag.log.eigen <- sapply(xs,function(x){
  BBsolve(fn=char.im,par=c(0.1,0.1),r=x,quiet=TRUE)$par
}
                      )

pdf("figure4-3.pdf")
matplot(x=xs,y=t(lag.log.eigen),type='l',ylim=c(-2,2),xlab="r",
        ylab=expression(lambda))
abline(v=0.4,lty=2)
dev.off()
