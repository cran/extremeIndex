#'@export
#'@import boot
.gp.pw.fit=function(data,init,CI=T,R=200,ncpus=1){
  x <- seq(0,max(data),by=0.05);
  dEGPs <- c();
  qEGPs <- c();

  fit.PWM <- .EGP.fitPWM(x=data,kappa0=init[1],sigma0=init[2],xi0=init[3]);
  if(CI){
    fit.PWM.boot <- boot::boot(data=data,statistic=.EGP.fitPWM.boot,R=R,kappa0=init[1],sigma0=init[2],xi0=init[3],parallel="multicore",ncpus=ncpus);
    CI.PWM.kappa <- boot::boot.ci(boot.out=fit.PWM.boot,index=1,type="perc")$perc[4:5];
    CI.PWM.sigma <- boot::boot.ci(boot.out=fit.PWM.boot,index=2,type="perc")$perc[4:5];
    CI.PWM.xi <- boot::boot.ci(boot.out=fit.PWM.boot,index=3,type="perc")$perc[4:5];
    CIs.PWM <- cbind(CI.PWM.kappa,CI.PWM.sigma,CI.PWM.xi);
  }





  fits <- list(PWM=fit.PWM)
  if(CI){
    CIs <- list(PWM=CIs.PWM)
  } else{
    CIs <- NULL
  }




  return(list(fit=fits,CI=CIs))
}
