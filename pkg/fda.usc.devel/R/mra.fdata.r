wt.filter.equivalent <- function(wt.filter, J){
  L <- wt.filter$L
  h <- wt.filter$h
  g <- wt.filter$g
  L.last <- L
  h.last <- h
  g.last <- g
  for(j in 2:J){
    L.new <- (2^j - 1)*(L-1) + 1
    hj <- NULL
    gj <- NULL
    for(l in 0:(L.new - 1)){
      u <- l
      ifelse(u >= L, g.mult <- 0, g.mult <- g[u+1])
      hjl <- g.mult*h.last[1]
      gjl <- g.mult*g.last[1]
      for(k in 1:(L.last-1)){
        u <- u-2
        if((u < 0) | (u >= L)) g.mult <- 0 else g.mult <- g[u+1]
        hjl <- hjl + g.mult*h.last[k+1]
        gjl <- gjl + g.mult*g.last[k+1]
      }
      hj <- c(hj,hjl)
      gj <- c(gj,gjl)
    }
    h.last <- hj
    g.last <- gj
    L.last <- L.new
  }
  wt.filter$L <- as.integer(L.last)
  wt.filter$h <- h.last
  wt.filter$g <- g.last
  wt.filter$level <- as.integer(J)
  return(wt.filter)
}

# el wf.filter.equivalent crea las bases en el linel deseado
# percival pag 95,96
################################################################################
qmf <- function(x) {
   l<-length(x)
   h <- (-1)^(0:(l - 1)) * x[l:1]
  return(h)
}
################################################################################
# The character strings currently supported are derived from one of four classes of wavelet transform
# filters: Daubechies, Least Asymetric, Best Localized and Coifiet. The pre?xes for filters of these
# classes are d, la, bl and c, respectively. Following the prefix, the filter name consists of an integer
# indicating length. Supported lengths are as follows:
# Daubechies 2,4,6,8,10,12,14,16,18,20.
# Least Asymetric 8,10,12,14,16,18,20.
# Best Localized 14,18,20.
# Coifiet 6,12,18,24,30
################################################################################
wave.filter <- function(name,level=1)    {
  filter.haar <- function() {
    L <- 2
    g <- c(0.7071067811865475, 0.7071067811865475)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.d4 <- function() {
    L <- 4
    g <- c(0.4829629131445341,0.8365163037378077,
           0.2241438680420134,-0.1294095225512603)/sqrt(2)
    h <- qmf(g)   
    return(list(L = L, h = h, g = g))
  }
  filter.mb4 <- function() {
    L <- 4
    g <- c(4.801755e-01, 8.372545e-01, 2.269312e-01, -1.301477e-01)/sqrt(2) 
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
# coifiet  
 filter.c6 <- function(){
      L <- as.integer(6)
      g <- c(-0.0156557285289848,-0.0727326213410511,
             0.3848648565381134,0.8525720416423900,
             0.3378976709511590,-0.0727322757411889)/sqrt(2)
     h <- qmf(g)   
     return(list(L = L, h = h, g = g))                
    }
    filter.c12 <- function(){
      L <- as.integer(12)
      g <- c(-0.0007205494453679,-0.0018232088707116,0.0056114348194211,
             0.0236801719464464,-0.0594344186467388,-0.0764885990786692,
             0.4170051844236707,0.8127236354493977,0.3861100668229939,
             -0.0673725547222826,-0.0414649367819558,0.0163873364635998)/sqrt(2)
     h <- qmf(g)   
     return(list(L = L, h = h, g = g)) 
    }

    filter.c18 <- function(){
      L <- as.integer(18)
      g <- c( -0.0000345997728362,-0.0000709833031381,0.0004662169601129,
             0.0011175187708906, -0.0025745176887502,-0.0090079761366615,
             0.0158805448636158,0.0345550275730615,-0.0823019271068856,
             -0.0717998216193117, 0.4284834763776168,0.7937772226256169,
             0.4051769024096150,-0.0611233900026726, -0.0657719112818552,
             0.0234526961418362, 0.0077825964273254,-0.0037935128644910)/sqrt(2)
     h <- qmf(g)   
     return(list(L = L, h = h, g = g)) 
    }                                                             
    filter.d20 <- function(){
      L <- as.integer(20)
      g <- c(0.0266700579005546,0.1881768000776863,0.5272011889317202,
             0.6884590394536250,0.2811723436606485,-0.2498464243272283,
             -0.1959462743773399,0.1273693403357890,0.0930573646035802,
             -0.0713941471663697,-0.0294575368218480,0.0332126740593703,
             0.0036065535669880-0.0107331754833036,0.0013953517470692,
             0.0019924052951930,-0.0006858566949566,-0.0001164668551285,
             0.0000935886703202,-0.0000132642028945)/sqrt(2)
     h <- qmf(g)   
     return(list(L = L, h = h, g = g)) 
    }
    
    filter.c24 <- function(){
      L <- as.integer(24)
      g <- c(-0.0000017849850031,-0.0000032596802369,
             0.0000312298758654,0.0000623390344610,
             -0.0002599745524878,-0.0005890207562444,
             0.0012665619292991,0.0037514361572790,
             -0.0056582866866115,-0.0152117315279485,
             0.0250822618448678,0.0393344271233433,
             -0.0962204420340021, -0.0666274742634348,
             0.4343860564915321,0.7822389309206135,
             0.4153084070304910, -0.0560773133167630,
             -0.0812666996808907, 0.0266823001560570,
             0.0160689439647787, -0.0073461663276432,
             -0.0016294920126020,0.0008923136685824)/sqrt(2)
     h <- qmf(g)   
     return(list(L = L, h = h, g = g)) 
    }
    filter.c30 <- function(){
      L <- as.integer(30)
      g <- c(-0.0000000951765727, -0.0000001674428858,
             0.0000020637618516, 0.0000037346551755,
             -0.0000213150268122, -0.0000413404322768,
             0.0001405411497166, 0.0003022595818445,
             -0.0006381313431115, -0.0016628637021860,
             0.0024333732129107, 0.0067641854487565,
             -0.0091642311634348, -0.0197617789446276,
             0.0326835742705106, 0.0412892087544753,
             -0.1055742087143175, 0.7742896037334738,
             0.4215662067346898,  -0.0520431631816557,
             -0.0919200105692549,  0.0281680289738655,
             0.0234081567882734,  -0.0101311175209033,
             -0.0041593587818186,  0.0021782363583355,
             0.0003585896879330, -0.0002120808398259)/sqrt(2)
     h <- qmf(g)   
     return(list(L = L, h = h, g = g)) 
    }    
  filter.fk4 <- function() {
    L <- 4
    g <- c(.6539275555697651, .7532724928394872, .5317922877905981e-1,
          -.4616571481521770e-1)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.d6 <- function() {
    L <- 6
    g <- c(0.3326705529500827, 0.8068915093110928, 0.4598775021184915,
          -0.1350110200102546, -0.0854412738820267, 0.0352262918857096)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.fk6 <- function() {
    L <- 6
    g <- c(.4279150324223103, .8129196431369074, .3563695110701871,
          -.1464386812725773, -.7717775740697006e-1, .4062581442323794e-1)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.d8 <- function() {
    L <- 8
    g <- c(0.2303778133074431, 0.7148465705484058, 0.6308807679358788,
          -0.0279837694166834, -0.1870348117179132, 0.0308413818353661,
          0.0328830116666778, -0.0105974017850021)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.fk8 <- function() {
    L <- 8
    g <- c(.3492381118637999, .7826836203840648, .4752651350794712,
          -.9968332845057319e-1, -.1599780974340301, .4310666810651625e-1,
          .4258163167758178e-1, -.1900017885373592e-1)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.la8 <- function() {
    L <- 8
    g <- c(-0.07576571478935668, -0.02963552764596039, 0.49761866763256290, 
	  0.80373875180538600, 0.29785779560560505, -0.09921954357695636, 
	  -0.01260396726226383, 0.03222310060407815)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }   
  filter.mb8 <- function() {
    L <- 8 
    aa<-c(-1.673619e-01, 1.847751e-02, 5.725771e-01, 7.351331e-01,           
	   2.947855e-01, -1.108673e-01, 7.106015e-03, 6.436345e-02)/sqrt(2)
    g <- aa[L:1]
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.bl14 <- function() {
    L <- 14
    g <- c( 0.0120154192834842, 0.0172133762994439, -0.0649080035533744,
	  -0.0641312898189170, 0.3602184608985549, 0.7819215932965554,
	   0.4836109156937821, -0.0568044768822707, -0.1010109208664125,
	   0.0447423494687405, 0.0204642075778225, -0.0181266051311065,
	  -0.0032832978473081, 0.0022918339541009)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.fk14 <- function() {
    L <- 14
    g <- c(.2603717692913964, .6868914772395985, .6115546539595115,
          .5142165414211914e-1, -.2456139281621916, -.4857533908585527e-1,
          .1242825609215128, .2222673962246313e-1, -.6399737303914167e-1,
          -.5074372549972850e-2, .2977971159037902e-1, -.3297479152708717e-2,
          -.9270613374448239e-2, .3514100970435962e-2)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.d16 <- function() {
    L <- 16
    g <- c(0.0544158422431049, 0.3128715909143031, 0.6756307362972904,
	  0.5853546836541907, -0.0158291052563816, -0.2840155429615702,
	  0.0004724845739124, 0.1287474266204837, -0.0173693010018083,
	 -0.0440882539307952, 0.0139810279173995, 0.0087460940474061,
	 -0.0048703529934518, -0.0003917403733770, 0.0006754494064506,
	 -0.0001174767841248)/sqrt(2) 	 
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.la10 <- function(){
      L <- as.integer(10)
      g <- c(0.0195388827353869,-0.0211018340249298,
             -0.1753280899081075,0.0166021057644243,
             0.6339789634569490,0.7234076904038076,
             0.1993975339769955,-0.0391342493025834,
             0.0295194909260734,0.0273330683451645)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
    }
    filter.la12 <- function(){
      L <- as.integer(12)
      g <- c(0.0154041093273377,0.0034907120843304,
             -0.1179901111484105,-0.0483117425859981,
             0.4910559419276396,0.7876411410287941,
             0.3379294217282401,-0.0726375227866000,
             -0.0210602925126954,0.0447249017707482,
             0.0017677118643983,-0.0078007083247650)/sqrt(2)
	   
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
    }   
    filter.la14 <- function(){
      L <- as.integer(14)
      g <- c(0.0102681767084968,0.0040102448717033,
             -0.1078082377036168,-0.1400472404427030,
             0.2886296317509833,0.7677643170045710,
             0.5361019170907720,0.0174412550871099,
             -0.0495528349370410,0.0678926935015971,
             0.0305155131659062,-0.0126363034031526,
             -0.0010473848889657,0.0026818145681164)/sqrt(2)    	   
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
    }     
  filter.la16 <- function() {
    L <- 16
    g <- c(-0.0033824159513594, -0.0005421323316355, 0.0316950878103452, 
	   0.0076074873252848, -0.1432942383510542, -0.0612733590679088, 
	   0.4813596512592012, 0.7771857516997478, 0.3644418948359564, 
	  -0.0519458381078751, -0.0272190299168137, 0.0491371796734768, 
	   0.0038087520140601, -0.0149522583367926, -0.0003029205145516, 
	   0.0018899503329007) /sqrt(2)
	   
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.mb16 <- function() {
    L <- 16
    aa<-c(-1.302770e-02, 2.173677e-02, 1.136116e-01, -5.776570e-02, 
	  -2.278359e-01,1.188725e-01, 6.349228e-01, 6.701646e-01, 
	   2.345342e-01,-5.656657e-02, -1.987986e-02, 5.474628e-02, 
	  -2.483876e-02,-4.984698e-02, 9.620427e-03, 5.765899e-03)/sqrt(2)
    g <- aa[L:1]
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.la20 <- function() {
    L <- 20
    g <- c(0.0007701598091030, 0.0000956326707837, -0.0086412992759401,
	 -0.0014653825833465, 0.0459272392237649, 0.0116098939129724,
	 -0.1594942788575307, -0.0708805358108615, 0.4716906668426588,
	  0.7695100370143388, 0.3838267612253823, -0.0355367403054689,
	 -0.0319900568281631, 0.0499949720791560, 0.0057649120455518,
	 -0.0203549398039460, -0.0008043589345370, 0.0045931735836703,
	  0.0000570360843390, -0.0004593294205481)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.bl20 <- function() {
    L <- 20
    g <- c(0.0008625782242896, 0.0007154205305517, -0.0070567640909701,
	  0.0005956827305406, 0.0496861265075979, 0.0262403647054251,
	 -0.1215521061578162, -0.0150192395413644, 0.5137098728334054,
	  0.7669548365010849, 0.3402160135110789, -0.0878787107378667,
	 -0.0670899071680668, 0.0338423550064691, -0.0008687519578684,
	 -0.0230054612862905, -0.0011404297773324, 0.0050716491945793,
	  0.0003401492622332, -0.0004101159165852)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.fk22 <- function() {
    L <- 22
    g <- c(.1938961077599566, .5894521909294277, .6700849629420265,
          .2156298491347700, -.2280288557715772, -.1644657152688429,
          .1115491437220700, .1101552649340661, -.6608451679377920e-1,
          -.7184168192312605e-1, .4354236762555708e-1, .4477521218440976e-1,
          -.2974288074927414e-1, -.2597087308902119e-1, .2028448606667798e-1,
          .1296424941108978e-1, -.1288599056244363e-1, -.4838432636440189e-2,
          .7173803165271690e-2, .3612855622194901e-3, -.2676991638581043e-2,
          .8805773686384639e-3)/sqrt(2)
    h <- qmf(g)
    return(list(L = L, h = h, g = g))
  }
  filter.mb24 <- function() {
    L <- 24
    aa<-c(-2.132706e-05, 4.745736e-04, 7.456041e-04, -4.879053e-03,                     
	  -1.482995e-03, 4.199576e-02, -2.658282e-03, -6.559513e-03,
	   1.019512e-01, 1.689456e-01, 1.243531e-01, 1.949147e-01,
	   4.581101e-01, 6.176385e-01, 2.556731e-01, -3.091111e-01,
	  -3.622424e-01, -4.575448e-03, 1.479342e-01, 1.027154e-02,
	  -1.644859e-02, -2.062335e-03, 1.193006e-03, 5.361301e-05)/sqrt(2)
    g <- aa[L:1]
    h <- qmf(g)
    return(list(L = L, h= h, g = g))
  }
 out<- switch(name,                        
    "haar" = filter.haar(),
    "d4" = filter.d4(),
    "d6" = filter.d6(),    
    "d8" = filter.d8(),    
    "d16" = filter.d16(),
    "d16" = filter.d20(),    
    "la8" = filter.la8(),
    "la14" = filter.la14(),    
    "la16" = filter.la16(),
    "la20" = filter.la20(),
    "fk4" = filter.fk4(),
    "fk6" = filter.fk6(),  
    "fk8" = filter.fk8(),
    "fk14" = filter.fk14(),
    "fk22" = filter.fk22(),
    "mb4" = filter.mb4(),
    "mb8" = filter.mb8(),
    "mb16" = filter.mb16(),
    "mb24" = filter.mb24(),
    "bl14" = filter.bl14(),
    "bl20" = filter.bl20(),
    "c6" = filter.c6(),    
    "c12" = filter.c12(),    
    "c18" = filter.c18(),            
    "c24" = filter.c24(),    
    "c30" = filter.c30(),    
    stop("Invalid selection for wave.filter"))
   if(level > 1) out <- wt.filter.equivalent(out, J=level)
   else out$level <- as.integer(1)
    return(out)
}    


dwt.forward <- function(V, filter, j){
  N <- length(V)
  h <- filter$h
  g <- filter$g
  L <- filter$L
  Wj <- rep(NA, length=N)
  Vj <- rep(NA, length=N)
  for(t in 0:(N-1)){
    k <- t
    Wjt <- h[1]*V[k+1]
    Vjt <- g[1]*V[k+1]
    for(n in 1:(L-1)){
      k <- k - 2^(j-1)
      if(k < 0) k <- k + ceiling(-k/N)*N
      Wjt <- Wjt + h[n+1]*V[k+1]
      Vjt <- Vjt + g[n+1]*V[k+1]
    }
    Wj[t+1] <- Wjt
    Vj[t+1] <- Vjt    
  }
  results <- list(W = Wj, V = Vj)
  return(results)
}               
################################################################################
dwt.backward <- function(W, V, filter, j){
  N <- length(V)
  h <- filter$h
  g <- filter$g
  L <- length(h)#filter$L
  Vj <- rep(NA, length=N)
  for(t in 0:(N-1)){
    k <- t
    Vjt <- h[1]*W[k+1] + g[1]*V[k+1]
    for(n in 1:(L-1)){
      k <- k + 2^(j-1)
      if(k >= N) k <- k - floor(k/N)*N
      Vjt <- Vjt + h[n+1]*W[k+1] + g[n+1]*V[k+1]
    }
    Vj[t+1] <- Vjt
  }
  return(Vj)
}


################################################################################
dwt <- function(X, filter="la8", n.levels){
  if (is.character(filter))   filter <- wave.filter(filter) 
  else  stop("Invalid argument: 'filter' must be of class 'character'")
  L <- filter$L
  X <- as.matrix(X)
  dim.X <- dim(X)
  N <- dim(X)[1]
  n.series <- dim(X)[2]   
  # determine the level of decomposition
  if(missing(n.levels)){
    J <- as.integer(floor(log(((N-1)/(L-1))+1)/log(2)))
  } else if(!is.numeric(n.levels) | (round(n.levels) != n.levels)){
    stop("Invalid argument value: 'n.levels' must be an integer value")
  } else J <- as.integer(n.levels)   
  # initialize variables for pyramid algorithm
  Vj <- X
  W.coefs <- as.list(rep(NA, length=J))
  names(W.coefs) <- lapply(1:J, function(j, x){names(x)[j] <- paste("W",j,sep="")}, x=W.coefs)
  V.coefs <- as.list(rep(NA, length=J))
  names(V.coefs) <- lapply(1:J, function(j, x){names(x)[j] <- paste("V",j,sep="")}, x=V.coefs)

  # implement the pyramid algorithm
  for(j in 1:J){  
      analysis <- sapply(1:n.series,
                         function(i,v,f,j){
                           out <- dwt.forward(v[,i],f,j)
                          return(c(out$W,out$V))
                         }, v=Vj, f=filter, j=j)
    Wj <- matrix(analysis[1:N,], ncol=n.series)
    Vj <- matrix(analysis[(N+1):(2*N),], ncol=n.series)
    W.coefs[[j]] <- Wj
    V.coefs[[j]] <- Vj
    Lj <- (2^j-1)*(L-1)+1
  }
   # creat dwt object for ouput
  dwt <- list(
               W = W.coefs,
               V = V.coefs,
               filter = filter,
               level = J,
               series = X)     
  return(dwt)
}
################################################################################
mra.fdata <- function(fdataobj, filter="la8", n.levels){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  X<-t(fdataobj$data)
  tt<-fdataobj$argvals
  rtt<-fdataobj$rangeval  
  J <- n.levels
  details <- as.list(rep(NA, length=J))
  names(details) <- lapply(1:J, function(j, x){names(x)[j] <- paste("D",j,sep="")}, x=details)
  smooths <- as.list(rep(NA, length=J))
  names(smooths) <- lapply(1:J, function(j, x){names(x)[j] <- paste("S",j,sep="")}, x=smooths)
  wt <- dwt(X, filter,n.levels)
  filter <- wt$filter
  n.series <- dim(wt$W[[1]])[2]
  N <- dim(wt$series)[1]    
  for(j in 1:J){
    Wj <- wt$W[[j]]
    Vj <- wt$V[[j]]
    DWj <- Wj
    SWj <- matrix(rep(0, length=length(Wj)), ncol=n.series)
    DVj <- matrix(rep(0, length=length(Wj)), ncol=n.series)
    SVj <- Vj    
    for(k in j:1){
          DVj <- sapply(1:n.series,
                        function(i,w,v,f,k){
                          return(out <- dwt.backward(w[,i],v[,i],f,k))
                        }, w=DWj, v=DVj, f=filter, k=k)
          SVj <- sapply(1:n.series,
                        function(i,w,v,f,k){
                          return(out <- dwt.backward(w[,i],v[,i],f,k))
                        }, w=SWj, v=SVj, f=filter, k=k)
      DWj <- matrix(rep(0, length=length(DVj)), ncol=n.series)
      SWj <- DWj
    }
    details[[j]] <- fdata(t(DVj),tt,rtt)
    smooths[[j]] <- fdata(t(SVj),tt,rtt)
  }
  mra <-list(fdata.est=smooths[[J]]+details[[J]], D =details,S = smooths,filter = filter,
            level = as.integer(J),fdataobj= fdataobj)            
 return(mra)
}
#crear basis con wt.filter.equivalent, <fx,basis>=scores =dij?

