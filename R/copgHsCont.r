copgHsCont <- function(p1, p2, teta, teta.st, VC){
     
########################################################################################   
# Transformations on theta parameter
########################################################################################   

   
if(VC$BivD %in% c("N") ) {

derteta.derteta.st <- 1/cosh(teta.st)^2
der2teta.derteta.stteta.st <- -(2 * (sinh(teta.st) * cosh(teta.st))/(cosh(teta.st)^2)^2)
       
}  

if(VC$BivD %in% c("F") ) {

derteta.derteta.st <- 1
der2teta.derteta.stteta.st <- 0
       
} 
   

if(VC$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) derteta.derteta.st <- der2teta.derteta.stteta.st <-  exp(teta.st) 
if(VC$BivD %in% c("C90","C270","J90","J270","G90","G270") ) derteta.derteta.st <- der2teta.derteta.stteta.st <- -exp(teta.st)  

  

########################################################################################
########################################################################################


########################################################################################
# Rotations
########################################################################################

if(VC$BivD %in% c("C90","J90","G90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(VC$BivD %in% c("C180","J180","G180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(VC$BivD %in% c("C270","J270","G270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################
     
     
     
     
     
     
if(VC$BivD == "N"){


der2h.derp2teta <- (((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - teta^2) * (dnorm((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - 
    teta^2)) * (qnorm(p2)/sqrt(1 - teta^2) - (qnorm(p1) - teta * qnorm(p2)) * (0.5 * 
    (2 * teta * (1 - teta^2)^-0.5))/sqrt(1 - teta^2)^2)) * sqrt(2 * 
    pi) * (-teta) - dnorm((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - teta^2)) * 
    sqrt(2 * pi))/sqrt(1 - teta^2) + dnorm((qnorm(p1) - teta * qnorm(p2))/sqrt(1 - 
    teta^2)) * sqrt(2 * pi) * (-teta) * (0.5 * (2 * teta * (1 - 
    teta^2)^-0.5))/sqrt(1 - teta^2)^2)/exp(-qnorm(p2)^2/2)
    
    
der2h.derp2p2 <-  -teta*(dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1-teta^2)))/
                  sqrt(1-teta^2)/dnorm(qnorm(p2))^2*(((qnorm(p1)-teta*qnorm(p2))/
                  sqrt(1-teta^2))/sqrt(1-teta^2)*teta + qnorm(p2))

der2h.derteta.teta.st <-(dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1-teta^2))*
                           (((qnorm(p1)-teta*qnorm(p2))*qnorm(p2)/(1-teta^2))-((qnorm(p1)-
                           teta*qnorm(p2))^2*teta/(1-teta^2)^2)))*((-qnorm(p2)*sqrt(1-teta^2)+
                           (qnorm(p1)-teta*qnorm(p2))*teta/sqrt(1-teta^2))/(1-teta^2)) + 
                           (dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1-teta^2)))*((2*qnorm(p1)*teta^2 - 
                           3*teta*qnorm(p2) + qnorm(p1))/((1-teta^2)^2.5)) 
  
der2h.derp1p2 <-   -((exp(-(teta^2*(qnorm(p2)^2+qnorm(p1)^2)-(2*teta*qnorm(p2)*qnorm(p1)))/
                     (1 - teta^2)/2))/sqrt(1 - teta^2))*(teta*(sqrt(2*pi)/exp(-qnorm(p2)^2/2))/
                    (1 - teta^2))*(teta*qnorm(p2)-qnorm(p1))

der2h.derp1teta <- (-2*(teta*(qnorm(p1)^2 + qnorm(p2)^2)-qnorm(p1)*qnorm(p2))*(1/(1 - teta^2)/2)-
                     (teta^2*(qnorm(p1)^2 + qnorm(p2)^2)-2*teta*qnorm(p1)*qnorm(p2))/(1 - teta^2)^2*
                     teta)*(exp(-(teta^2*(qnorm(p1)^2 + qnorm(p2)^2)-2*teta*qnorm(p1)*qnorm(p2))*
                     (1/(1 - teta^2)/2)))/sqrt(1 - teta^2)+(exp(-(teta^2*(qnorm(p1)^2 + qnorm(p2)^2)-
                     2*teta*qnorm(p1)*qnorm(p2))*(1/(1 - teta^2)/2)))/sqrt(1 - teta^2)/(1 - teta^2)*teta
                                                            

 
der2h.derp1p1 <-   -((exp(-(teta^2*(qnorm(p1)^2+qnorm(p2)^2)-(2*teta*qnorm(p1)*qnorm(p2)))/
	             (1 - teta^2)/2))/sqrt(1 - teta^2))*(teta*(sqrt(2*pi)/exp(-qnorm(p1)^2/2))/
	            (1 - teta^2))*(teta*qnorm(p1)-qnorm(p2))    


}




        


if(VC$BivD %in% c("C0","C90","C180","C270")){

der2h.derp2teta <- -(1/((-p2^teta + p1^teta* (-1 + p2^teta))^3 *teta^2))*
 p1^teta* p2^(-2 + teta)* (-1 + p1^-teta + p2^-teta)^(-1/
   teta)* (teta* (1 + 
       teta)* (p1^teta* teta + (-1 + p1^teta)* p2^teta* (1 + teta))* log(
      p1) + (-1 + p1^
       teta)* (teta* (1 + teta)* ((-1 + p1^teta)* p2^teta* teta + 
          p1^teta* (1 + teta))* log(
         p2) - (-p1^teta + (-1 + p1^teta)* p2^
           teta)* (teta^2 + (1 + teta) *log(-1 + p1^-teta + p2^-teta))))
    
    
t1 = -teta-1
t2 = p2^t1
t3 = t1^2
t5 = p2^2
t6 = 1/t5
t7 = p1^-teta
t8 = p2^-teta
t9 = t7+t8-1
t11 = -1-1/teta
t12 = t9^t11
t13 = t6*t12
t16 = t2*t1*t13
t18 = 1/t9
t23 = t2*t12
t24 = t11*t11
t26 = t8*t8
t27 = teta^2
t29 = t9*t9
t32 = t26*t27*t6/t29
t34 = t23*t11
t36 = t6*t18

    der2h.derp2p2 <- t2*t3*t13-t16-2*t16*t11*t8*teta*t18+t23*t24*t32+t34*t8*t27*t36+t34*t8*teta*t36-t34*t32


t2 = p2^(-1.0*teta-1.0);
t3 = log(p2);
t4 = t3*t3;
t6 = p1^(-1.0*teta);
t7 = p2^(-1.0*teta);
t8 = t6+t7-1.0;
t10 = -1.0-1/teta;
t11 = t8^(1.0*t10);
t14 = teta^2;
t15 = 1/t14;
t16 = log(t8);
t18 = log(p1);
t21 = -t6*t18-t7*t3;
t23 = 1/t8;
t25 = t15*t16+t10*t21*t23;
t29 = t2*t11;
t32 = t25*t25;
t12 = t18*t18;
t9 = t21*t21;
t5 = t8*t8;

der2h.derteta.teta.st <- t2*t4*t11-2.0*t2*t3*t11*t25+t29*t32+t29*(-2.0/t14/teta*t16+2.0*t15*t21*t23+t10*(t6*t12+t7*t4)*t23-t10*t9/t5);




t1 = 1.0+teta;
t3 = (p2*p1)^(-1.0*t1);
t4 = t1*t3;
t5 = 1/p2;
t7 = p2^(-1.0*teta);
t8 = p1^(-1.0*teta);
t9 = t7+t8-1.0;
t11 = -2.0-1/teta;
t12 = t9^(1.0*t11);

der2h.derp1p2 <-    -t4*t1*t5*t12-t4*t12*t11*t7*teta*t5/t9;


t1 = p1*p2;
t2 = -teta-1.0;
t3 = t1^(1.0*t2);
t4 = p1^(-1.0*teta);
t5 = p2^(-1.0*teta);
t6 = t4+t5-1.0;
t7 = -2.0-1/teta;
t8 = t6^(1.0*t7);
t9 = -t2*t3;
t10 = log(t1);
t11 = teta^2;
t12 = log(t6);
t13 = log(p1);
t14 = log(p2);

der2h.derp1teta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);






t1 = 1.0+teta;
t3 = (p2*p1)^(-1.0*t1);
t4 = t1*t3;
t5 = 1/p1;
t7 = p1^(-1.0*teta);
t8 = p2^(-1.0*teta);
t9 = t7+t8-1.0;
t11 = -2.0-1/teta;
t12 = t9^(1.0*t11);

der2h.derp1p1 <-    -t4*t1*t5*t12-t4*t12*t11*t7*teta*t5/t9;


}





if(VC$BivD == "F"){

der2h.derp2teta <-  -(1/((-exp((p1 + p2) *teta) + 
   exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^3))*
 exp(teta + p2* teta)* (-exp((1 + 2* p1 + p2)* teta)* (2 + (1 + p1 - 2* p2)* teta) + 
    exp((2 + p1)* teta)* (-2 + p1 *teta - 2* p2* teta) + 
    exp(teta + 3* p1* teta)* (-1 + teta - p2* teta) + 
    exp((3* p1 + p2)* teta)* (1 + teta - p2 *teta) + 
    exp((2 + p1 + p2)* teta)* (1 + p1* teta - p2* teta) + 
    exp((2 + p2)* teta)* (-1 + p2* teta) + exp(2* teta)* (1 + p2* teta) + 
    exp(2* (1 + p1)* teta)* (1 - p1* teta + p2* teta) + 
    exp((2* p1 + p2)* teta)* (-1 + (-1 + p1 + p2)* teta) - 
    exp(teta + p1* teta)* (1 + (-1 + p1 + p2)* teta) + 
    exp(teta + 2* p1* teta)* (2 + (-2 + p1 + 2 *p2)* teta) + 
    exp((1 + p1 + p2) *teta)* (2 - (-1 + p1 + 2* p2)* teta))
        
        
        
t1 = exp(teta)
t2 = teta*p1
t3 = exp(t2)
t5 = t1*(t3-1)
t6 = teta*p2
t8 = exp(t6+t2)
t10 = exp(t6+teta)
t12 = exp(t2+teta)
t13 = t8-t10-t12+t1
t14 = t13*t13
t20 = (teta*t8-teta*t10)^2
t24 = teta^2

       der2h.derp2p2 <- -2.0*t5/t14/t13*t20+t5/t14*(t24*t8-t24*t10)


t1 = exp(teta);
t2 = teta*p1;
t3 = exp(t2);
t5 = t1*(t3-1.0);
t6 = teta*p2;
t8 = exp(t6+t2);
t10 = exp(t6+teta);
t12 = exp(t2+teta);
t13 = t8-t10-t12+t1;
t14 = 1/t13;
t16 = t1*p1;
t18 = t3*t14;
t20 = t13*t13;
t21 = 1/t20;
t23 = p2+p1;
t25 = p2+1.0;
t26 = p1+1.0;
t28 = t23*t8-t25*t10-t26*t12+t1;
t32 = p1*p1;
t42 = t28*t28;
t44 = t23*t23;
t47 = t25*t25;
t49 = t26*t26;

der2h.derteta.teta.st <- -t5*t14-2.0*t16*t18+2.0*t5*t21*t28-t1*t32*t18+2.0*t16*t3*t21*t28-2.0*t5/t20/t13*t42+t5*t21*(t44*t8-t47*t10-t49*t12+t1);



			t1 = teta^2;
			t2 = exp(teta);
			t3 = t2 - 1;
			t5 = teta*p1;
			t6 = teta*p2;
			t8 = exp(t5+t6+teta);
			t10 = exp(t5+t6);
			t12 = exp(t5+teta);
			t14 = exp(t6+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			
der2h.derp1p2 <- t1*t3*t8/t16-2*teta*t3*t8/t16/t15*(teta*t10-teta*t14);




			t2 = exp(teta);
			t3 = t2-1.0;
			t4 = teta*p2;
			t5 = teta*p1;
			t7 = exp(t4+t5+teta);
			t10 = exp(t4+t5);
			t12 = exp(t4+teta);
			t14 = exp(t5+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			t17 = 1/t16;
			t21 = teta*t3;
der2h.derp1teta <- t3*t7*t17+teta*t2*t7*t17+t21*(p2+p1+1.0)*t7*t17-2.0*t21*t7/t15/t16*((p2+p1)*t10-(p2+1.0)*t12-(p1+1.0)*t14+t2);
		


			t1 = teta^2;
			t2 = exp(teta);
			t3 = t2 - 1;
			t5 = teta*p2;
			t6 = teta*p1;
			t8 = exp(t5+t6+teta);
			t10 = exp(t5+t6);
			t12 = exp(t5+teta);
			t14 = exp(t6+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			
der2h.derp1p1 <-  t1*t3*t8/t16-2*teta*t3*t8/t16/t15*(teta*t10-teta*t14);




}



if(VC$BivD %in% c("G0","G90","G180","G270")){

 der2h.derp2teta <-  -(1/(p2^2* teta^2))*
 exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
   teta))) *((-log(p1))^teta + (-log(p2))^teta)^(-3 + 1/
   teta)* (-log(p2))^(-2 + 
   teta)* (teta* (-log(p1))^
     teta* (-(-log(p1))^
         teta* (-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
          teta))* (-1 + teta - log(p2)) + (-log(p2))^
        teta* (((-log(p1))^teta + (-log(p2))^teta)^(
          2/teta) + (-1 + teta)* (teta + log(p2)) + ((-log(p1))^
             teta + (-log(p2))^teta)^(1/
           teta)* (-2 + 2 *teta + log(p2))))* log(-log(
        p1)) + ((-log(p1))^teta + (-log(p2))^
       teta)* ((-log(p1))^
        teta* (-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/teta))* (-1 +
           teta - log(p2)) + (-log(p2))^
        teta *(-((-log(p1))^teta + (-log(p2))^teta)^(
           2/teta) - ((-log(p1))^teta + (-log(p2))^teta)^(1/
           teta)* (-2 + log(p2)) + log(p2)))* log((-log(p1))^
       teta + (-log(p2))^teta) + 
    teta* ((((-log(p1))^teta + (-log(p2))^teta)^(
          2/teta) + ((-log(p1))^teta + (-log(p2))^teta)^(1/
           teta)* (-2 + log(p2)) - log(p2))* (-log(p2))^(2* teta)*
         log(-log(p2)) + 
       teta* (-log(p1))^(
        2* teta)* (1 + (-1 + teta - log(p2))* log(-log(p2))) + (-log(
           p1))^teta* (-log(p2))^
        teta* (teta - (1 - 2* teta + 
             teta^2 + ((-log(p1))^teta + (-log(p2))^teta)^(1/
              teta)* (-1 + 3 *teta - log(p2)) + log(p2) + 
             teta* log(p2))* log(-log(p2)))))
        
  
 t1 = log(p2)
t2 = (-t1)^teta 
t3 = log(p1)
t4 = (-t3)^teta  
t5 = t2+t4
t6 = 1/teta
t7 = t5^t6
t8 = exp(-t7)
t9 = t6-1
t10 = t5^t9
t11 = t8*t10
t12 = p2^2
t14 = 1/t12/p2
t15 = t2*t14
t20 = t1*t1
t21 = 1/t20
t26 = 1/t20/t1
t30 = t7*t7
t31 = t2*t2
t32 = t31*t2
t34 = t5*t5
t36 = 1/t34
t37 = t26*t36
t38 = t37*t11
t40 = t11*t2
t41 = teta*t14
t48 = t7*t31
t49 = t48*t14
t50 = 1/t5
t51 = t21*t50
t55 = t26*t50
t59 = t7*t32
t62 = t14*t26
t65 = t10*teta
t69 = t59*t62
t70 = t36*t8
t79 = t11*t9*t31
t86 = t9*t9
t18 = teta^2
t16 = t18*t14
t13 = t16*t37
 

 
der2h.derp2p2 <- -2*t11*t15/t1-3*t11*t15*t21-2*t11*t15*t26-t30*t32*t14*t38+3*t40*t41*t21+3*t40*t41*t26-3*t49*t51*t11-3*t49*t55*t11+t59*t14*t38+3*t48*t62*t50*t8*t65-t69*t70*t65+2*t69*t70*t10*t9*teta+3*t79*t41*t51+3*t79*t41*t55-t11*t86*t32*t13-3*t79*t16*t55+t11*t9*t32*t13-t40*t16*t26 
  
   
  
			  t1 = log(p2);
			  t2 = (-t1)^teta
			  t3 = log(p1);
			  t4 = (-t3)^teta
			  t5 = t2+t4;
			  t6 = 1/teta;
			  t7 = t5^t6
			  t8 = teta^2
			  t9 = 1/t8;
			  t10 = log(t5);
			  t11 = t9*t10;
			  t12 = log(-t1);
			  t13 = t2*t12;
			  t14 = log(-t3);
			  t16 = t13+t4*t14;
			  t18 = 1/t5;
			  t20 = -t11+t6*t16*t18;
			  t21 = t20^2
			  t23 = exp(-t7);
			  t25 = t6-1.0;
			  t26 = t5^t25
			  t28 = 1/p2
			  t29 = 1/t1;
			  t30 = t28*t29;
			  t31 = t26*t2*t30;
			  t36 = 2.0/t8/teta*t10;
			  t39 = 2.0*t9*t16*t18;
			  t40 = t12^2;
			  t42 = t14^2;
			  t44 = t2*t40+t4*t42;
			  t47 = t16^2
			  t49 = t5^2;
			  t50 = 1/t49;
			  t56 = t7*t7;
			  t61 = t23*t26;
			  t62 = t7*t20*t61;
			  t65 = -t11+t25*t16*t18;
			  t70 = t13*t30;
			  t73 = t65*t65;
			  t15 = t2*t28*t29;
			  
der2h.derteta.teta.st <- t7*t21*t23*t31+t7*(t36-t39+t6*t44*t18-t6*t47*t50)*t23*t31-t56*t21*t23*t31+2.0*t62*t65*t2*t30+2.0*t62*t70-t61*t73*t15-t61*(t36-t39+t25*t44*t18-t25*t47*t50)*t15-2.0*t61*t65*t70-t61*t2*t40*t28*t29;  
  
  
  
  
  
  
  			t3 = log(p2);
  			t4 = (-t3)^teta;
  			t5 = log(p1);
  			t6 = (-t5)^teta
  			t7 = t4+t6;
  			t8 = 1/teta;
  			t9 = t7^t8
  			t11 = p2^2;
  			t12 = 1/t11;
  			t13 = 1/t3;
  			t15 = 1/t7;
  			t18 = exp(-t9);
  			t19 = 1/p1;
  			t21 = -1 + t8;
  			t22 = (t7)^(2*t21)
  			t24 = teta - 1
  			t25 = (t3*t5)^t24
  			t27 = t7^-t8
  			t28 = t24*t27;
  			t29 = 1 + t28;
  			t30 = t22*t25*t29;
  			t33 = t18*t12;
  			t36 = t19*t22;
  			
der2h.derp1p2 <-  -t9*t4*t12*t13*t15*t18*t19*t30-t33*t19*t30+2.0*t33*t36*t21*t4*teta*t13*t15*t25*t29+t33*t36*t25*t24*t13*t29-t33*t36*t25*t28*t4*t13*t15;
  
  
  
  
  			t3 = log(p1);
  			  t4 = (-t3)^teta
  			  t5 = log(p2);
  			  t6 = (-t5)^teta
  			  t7 = t4+t6;
  			  t8 = 1/teta;
  			  t9 = t7^t8
  			  t10 = teta^2;
  			  t12 = log(t7);
  			  t13 = 1/t10*t12;
  			  t14 = log(-t3);
  			  t16 = log(-t5);
  			  t18 = t4*t14+t6*t16;
  			  t20 = 1/t7;
  			  t22 = -t13+t8*t18*t20;
  			  t24 = exp(-t9);
  			  t26 = t24/p1;
  			  t28 = 1/p2;
  			  t29 = -1.0+t8;
  			  t30 = t7^(2*t29)
  			  t32 = t3*t5;
  			  t33 = teta-1.0;
  			  t34 = t32^t33
  			  t35 = t7^(-1.0*t8)
  			  t36 = t33*t35;
  			  t17 = 1.0+t36;
  			  t15 = t34*t17;
  			  t11 = t26*t28;
  			  t2 = t30*t34;
  			  t1 = log(t32);
  			  
 der2h.derp1teta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
  

  			t3 = log(p1);
  			t4 = (-t3)^teta;
  			t5 = log(p2);
  			t6 = (-t5)^teta
  			t7 = t4+t6;
  			t8 = 1/teta;
  			t9 = t7^t8
  			t11 = p1^2;
  			t12 = 1/t11;
  			t13 = 1/t3;
  			t15 = 1/t7;
  			t18 = exp(-t9);
  			t19 = 1/p2;
  			t21 = -1 + t8;
  			t22 = (t7)^(2*t21)
  			t24 = teta - 1
  			t25 = (t3*t5)^t24
  			t27 = t7^-t8
  			t28 = t24*t27;
  			t29 = 1 + t28;
  			t30 = t22*t25*t29;
  			t33 = t18*t12;
  			t36 = t19*t22;
  			
der2h.derp1p1 <-  -t9*t4*t12*t13*t15*t18*t19*t30-t33*t19*t30+2.0*t33*t36*t21*t4*teta*t13*t15*t25*t29+t33*t36*t25*t24*t13*t29-t33*t36*t25*t28*t4*t13*t15;  




}






if(VC$BivD %in% c("J0","J90","J180","J270")){

 der2h.derp2teta <-  (1/(teta^2))*(1 - p1)^teta *((1 - p1)^
   teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-3 + 1/
  teta)* (1 - p2)^(-2 + 
  teta)* ((-1 + 
      teta)* teta* (-(-1 + (1 - p1)^teta)* (1 - p2)^
       teta* ((1 - p1)^teta - teta) + (1 - p1)^
       teta* (-1 + (1 - p1)^teta + teta))* log(
     1 - p1) + (-1 + (1 - p1)^
      teta)* ((-(1 - p1)^teta + (-1 + (1 - p1)^teta)* (1 - p2)^
          teta)* (-1 + teta)* log((1 - p1)^
         teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta) + 
      teta *((1 - p1)^
          teta* teta* (1 + (-1 + teta)* log(1 - p2)) + (-1 + (1 - p1)^
            teta)* (1 - p2)^teta *(-teta + (-1 + teta)^2 *log(1 - p2)))))
 
			t2 = (1 - p1)^teta
			t3 = 1 - p2
			t4 = t3^teta
			t5 = t2*t4
			t6 = t2+t4-t5
			t8 = 1/teta - 1
			t9 = t6^t8
			t10 = t8^2
			t12 = t4*teta
			t13 = 1/t3
			t14 = -t12*t13+t5*teta*t13
			t18 = t14^2
			t20 = t6^2
			t22 = teta - 1
			t23 = t3^t22
			t24 = 1 - t2
			t26 = 1/t20*t23*t24
			t27 = t9*t8
			t29 = teta^2
			t31 = t3*t3
			t32 = 1/t31
			t41 = 1/t6
			t51 = t9*t23
			t55 = t22*t22
			
der2h.derp2p2 <- t9*t10*t18*t26+t27*(t4*t29*t32-t12*t32-t5*t29*t32+t5*teta*t32)*t41*t23*t24-t27*t18*t26-2.0*t27*t14*t41*t23*t22*t13*t24+t51*t55*t32*t24-t51*t22*t32*t24 
 
 
 


			  t1 = 1.0-p1;
			  t2 = t1^teta
			  t3 = 1.0-p2;
			  t4 = t3^teta;
			  t5 = t2*t4;
			  t6 = t2+t4-t5;
			  t8 = 1/teta-1.0;
			  t9 = t6^t8;
			  t10 = teta^2
			  t11 = 1/t10;
			  t12 = log(t6);
			  t14 = log(t1);
			  t15 = t2*t14;
			  t16 = log(t3);
			  t18 = t4*t16;
			  t20 = t15+t18-t15*t4-t5*t16;
			  t21 = 1/t6;
			  t23 = -t11*t12+t8*t20*t21;
			  t25 = t23^2;
			  t28 = t3^(teta-1)
			  t29 = 1 - t2;
			  t30 = t28*t29;
			  t39 = t14^2;
			  t40 = t2*t39;
			  t42 = t16^2;
			  t50 = t20^2;
			  t56 = t6^2;
			  t13 = t9*t23;
			  t7 = t9*t28;
			  
			  
der2h.derteta.teta.st <- t9*t25*t30+t9*(2.0/t10/teta*t12-2.0*t11*t20*t21+t8*(t40+t4*t42-t40*t4-2.0*t15*t18-t5*t42)*t21-t8*t50/t56)*t30+2.0*t13*t28*t16*t29-2.0*t13*t28*t2*t14+t7*t42*t29-2.0*t7*t16*t2*t14-t7*t40;







t1 = 1.0-p2;
t2 = t1^(1.0*teta);
t3 = 1.0-p1;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t11 = t2*teta;
t12 = 1/t1;
t16 = -t11*t12+t11*t12*t4;
t19 = teta-1.0;
t20 = t1^(1.0*t19);
t22 = t3^(1.0*t19);
t23 = teta-1.0+t2+t4-t5;
t27 = t9*t20;

der2h.derp1p2 <- t9*t8*t16/t6*t20*t22*t23-t27*t19*t12*t22*t23+t27*t22*t16;




t1 = 1.0-p1;
t2 = t1^(1.0*teta);
t3 = 1.0-p2;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t10 = teta^2;
t11 = log(t6);
t12 = log(t1);
t13 = t2*t12;
t14 = log(t3);
t15 = t4*t14;
t16 = t13*t4;
t19 = t5*t14;
t21 = teta-1.0;
t27 = t1^(1.0*t21);
t28 = t3^(1.0*t21);
t30 = teta-1.0+t2+t4-t5;
t33 = t9*t27;

der2h.derp1teta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
            

t1 = 1.0-p1;
t2 = t1^(1.0*teta);
t3 = 1.0-p2;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t11 = t2*teta;
t12 = 1/t1;
t16 = -t11*t12+t11*t12*t4;
t19 = teta-1.0;
t20 = t1^(1.0*t19);
t22 = t3^(1.0*t19);
t23 = teta-1.0+t2+t4-t5;
t27 = t9*t20;

der2h.derp1p1 <- t9*t8*t16/t6*t20*t22*t23-t27*t19*t12*t22*t23+t27*t22*t16;


}




if( VC$BivD %in% c("C270","J270","G270") ) {

der2h.derp1p2   <- -der2h.derp1p2
der2h.derp1teta <- -der2h.derp1teta

}


if( VC$BivD %in% c("C90","J90","G90") ) {

der2h.derp1p1         <- -der2h.derp1p1
der2h.derp2p2         <- -der2h.derp2p2
der2h.derp1teta       <- -der2h.derp1teta
der2h.derteta.teta.st <- -der2h.derteta.teta.st 

}

if( VC$BivD %in% c("C180","J180","G180") ) {

der2h.derp2p2              = -der2h.derp2p2
der2h.derteta.teta.st      = -der2h.derteta.teta.st
derteta.derteta.st         = -derteta.derteta.st
der2teta.derteta.stteta.st = -der2teta.derteta.stteta.st 
der2h.derp1p2              = -der2h.derp1p2 
der2h.derp1teta            = -der2h.derp1teta                                   
der2h.derp2teta            = -der2h.derp2teta 
der2h.derp1p1              = -der2h.derp1p1

}




list(
der2h.derp2p2              = der2h.derp2p2, 
der2h.derteta.teta.st      = der2h.derteta.teta.st,  
derteta.derteta.st         = derteta.derteta.st, 
der2teta.derteta.stteta.st = der2teta.derteta.stteta.st,  
der2h.derp1p2              = der2h.derp1p2,  
der2h.derp1teta            = der2h.derp1teta,                                     
der2h.derp2teta            = der2h.derp2teta,  
der2h.derp1p1              = der2h.derp1p1)     




}




     























