#' cross point poly(A) sites using expectation maximization and laplace mixture model
#' @name LMMcutoff
#' @usage LMMcutoff(data,GenomeSize,p.adjusted)
#' @param data A data.frame containing seqname, strand, poly(A) sites, and gene id.
#' @param GenomeSize A genome size of a species recommended.
#' @param p.adjusted Adjusted pvalue (BH).
#' @return The polyA data.frame will be return cross point between cleavage sites within or different PolyA site.

#'@export
LMMcutoff<-function(data,GenomeSize,p.adjusted){
  ##### pre-process
  TMR_M_1<-sum(data$score)  # Total number of uniquely mapped reads
  genome_size <-GenomeSize      # Genome size in base pairs
  # Calculate lambda
  lambda_MM <- TMR_M_1/(2*genome_size)
  #####
  #################################count per million transformation
  cpm<-function(x,y){
    c<-x*(1000000/y)
    return(c)
  }
  ####
  ###########logit2 transformation
  logit<-function(x){
    y = (log2(x+(median(x))))
    return(y)
  }
  #########################
  result11<-cbind(data,pCount=round(cpm(data$score,sum(data$score)),0))
  #
  sdata<-result11%>%
    mutate(pvalue=dpois(pCount,lambda_MM))%>%
    ungroup()

  sdata<-cbind(sdata,adPvalue=p.adjust(sdata$pvalue,method="BH"))
  significant_p_data<-sdata%>%
    filter(adPvalue<p.adjusted)
  ###
  result1<-significant_p_data%>%
    group_by(seqname,strand,gene_id)%>%
    arrange(coord) %>%
    mutate(distance = c(0, diff(coord))) %>%
    ungroup()
  ##
  result_m<- result1 %>%
    filter(distance > 0)
  #####
  q1<-quantile(result1$distance)[2]
  q2<-quantile(result_m$distance)[2]
  #####################
  if((q2-q1)==0){
    distance1<-result$distance
    m<-median(distance1)
  }else{
    distance1<-result_m$distance
    m<-median(distance1)
  }
  #########################
  distance2<-result1$distance
  ######################
  if(m<=4){
    dis<-log2(distance1+0.5)
    logit_val<-0.5
  }else{
    if(m>4 && m<12){
      dis<-logit(distance2)
      logit_val<-median(distance2)
    }else{
      dis<-log2(distance1)
      logit_val<-0
    }
  }
  ######
  dis.kmeans <- kmeans(dis, 2)
  dis.kmeans.cluster <-dis.kmeans$cluster

  dis.df <- data.frame(x =dis, cluster =dis.kmeans.cluster)

  dis.summary.df <- dis.df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x),size = n())

  dis.summary.df <-dis.summary.df %>%
    mutate(alpha = size / sum(size))

  #######
  scale<-function(x){
    s<-abs(sqrt(x/2))
  }
  #######
  dis.summary.df <- dis.summary.df %>%
    mutate(sl=scale(variance))

  # Expectation Step of the EM Algorithm
  e_step <- function(x, mu.vector,sl.vector,alpha.vector) {
    comp1.prod <- dlaplace(x, location =mu.vector[1], scale = sl.vector[1])* alpha.vector[1]
    comp2.prod <- dlaplace(x, location =mu.vector[2], scale = sl.vector[2])* alpha.vector[2]
    sum.of.comps <- comp1.prod + comp2.prod
    comp1.post <- comp1.prod / sum.of.comps
    comp2.post <- comp2.prod / sum.of.comps

    sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
    sum.of.comps.ln.sum <- sum(sum.of.comps.ln)

    list("loglik" = sum.of.comps.ln.sum,
         "posterior.df" = cbind(comp1.post, comp2.post))
  }

  # Maximization Step of the EM Algorithm

  m_step <- function(x, posterior.df) {
    comp1.n <- sum(posterior.df[, 1])
    comp2.n <- sum(posterior.df[, 2])

    comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
    comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)

    comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
    comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n
    comp1.sl <- sqrt(comp1.var/2)
    comp2.sl <- sqrt(comp2.var/2)
    comp1.alpha <- comp1.n / length(x)
    comp2.alpha <- comp2.n / length(x)

    list("mu" = c(comp1.mu, comp2.mu),
         "sl"=c(comp1.sl, comp2.sl),
         "alpha" = c(comp1.alpha, comp2.alpha))
  }
  #Now we just need to write a loop to go between the functions for each EM step. Each iteration will consist of us first calling the e_step function and then calling the m_step function (if needed). We will run this for 50 iterations or when the log likelihood difference between two iteration is less than 1e-6 (whichever comes first):

  for (i in 1:50) {
    if (i == 1) {
      # Initialization
      e.step <- e_step(dis, dis.summary.df[["mu"]],dis.summary.df[["sl"]],
                       dis.summary.df[["alpha"]])
      m.step <- m_step(dis, e.step[["posterior.df"]])
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step <- e_step(dis, m.step[["mu"]], m.step[["sl"]],
                       m.step[["alpha"]])
      m.step <- m_step(dis, e.step[["posterior.df"]])
      loglik.vector <- c(loglik.vector, e.step[["loglik"]])

      loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
      if(loglik.diff < 1e-6) {
        break
      } else {
        cur.loglik <- e.step[["loglik"]]
      }
    }
  }
  loglik.vector

  # Plot a Mixture Component

  plot_mix_comps <- function(x, mu, scale, lam) {
    lam * dlaplace(x, location =mu, scale =scale)
  }

  the.plot<-data.frame(x = dis) %>%
    ggplot() +
    geom_histogram(aes(x, after_stat(density)), binwidth = 1, colour = "white",
                   fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(m.step$mu[1], m.step$sl[1],
                              lam = m.step$alpha[1]),
                  colour = "blue", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(m.step$mu[2],m.step$sl[2],
                              lam = m.step$alpha[2]),
                  colour = "red", lwd = 1.5) +
    ylab("Density") +
    xlab(expression("Distance between adjacent cleavage sites (logit"["2"]*"nt)"))+
    ggtitle("")

  #extract coordinates of both lines from plot
  line.df = data.frame(x = ggplot_build(the.plot)$data[[2]]$x,
                       red = ggplot_build(the.plot)$data[[2]]$y,
                       blue = ggplot_build(the.plot)$data[[3]]$y)


  #find the minimal distance between lines along y axis
  line.df$delta = line.df$red - line.df$blue

  #find x coordinate for the minimal delta y
  x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
  x_coord

  plot_cutoff<-the.plot+
    ylab("Density") +
    theme_classic()

  ########
  cutoffPoint=((2^x_coord)-logit_val)
  return(list(cutoffFig =plot_cutoff, cutoffPoint = cutoffPoint))
}

