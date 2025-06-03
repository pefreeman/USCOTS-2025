library(shiny)
library(shinythemes)
library(ggplot2)

ui <- navbarPage("Die Tosses",
  
  tabPanel("Analyze Data",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Input: your experimental data."),
        p("Output: p-values for the chi-square GoF test and for multinomial simulations."),
        numericInput("num.1","Number of 1's",5),
        numericInput("num.2","Number of 2's",5),
        numericInput("num.3","Number of 3's",5),
        numericInput("num.4","Number of 4's",5),
        numericInput("num.5","Number of 5's",5),
        numericInput("num.6","Number of 6's",5),
        numericInput("num.sim0","Number of Simulations",100000),
        numericInput("seed0","Random Number Seed",101),
        submitButton("Submit")
      ),
      mainPanel(
        width = 8,
        verbatimTextOutput("oid0")
      )
    )
  ),

  tabPanel("Single Experiment",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Input: simulation settings."),
        p("Output: p-values for the chi-square GoF test and for multinomial simulations."),
        
        numericInput("trials1","Number of Tosses",30),
        numericInput("faces1","Number of Die Faces",6),
        numericInput("num.sim1","Number of Simulations",100000),
        numericInput("seed1","Random Number Seed",101),
        submitButton("Submit")
      ),
      mainPanel(
        width = 8,
        verbatimTextOutput("oid1")
      )
    )
  ),

  tabPanel("Repeated Experiments",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Input: simulation settings."),
        p("Output: information about the difference in p-values between the chi-square GoF test and multinomial simulations."),
        numericInput("trials2","Number of Tosses",30),
        numericInput("faces2","Number of Die Faces",6),
        numericInput("num.sim2","Number of Simulations",1000),
        numericInput("num.rep2","Number of Repetitions",1000),
        numericInput("seed2","Random Number Seed",101),
        submitButton("Submit")
      ),
      mainPanel(
        width = 8,
        verbatimTextOutput("oid2"),
        plotOutput("plot2")
      )
    )
  ),


  # change font
  tags$head(tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=Manrope:wght@200..800&display=swap');
    * { font-family: 'Manrope'; }
  ")))
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  output$oid0 <- renderPrint({
    d <- c(input$num.1,input$num.2,input$num.3,input$num.4,input$num.5,input$num.6)
    k <- sum(d)
    m <- length(d)
    p <- rep(1/m,m)
    set.seed(input$seed0)
    cat("The total number of tosses         ",k,"\n\n")
    W <- sum((d-k*p)^2/(k*p))
    p.val <- 1-pchisq(W,length(d)-1)
    cat("The chi-square value:              ",round(W,3),"\n")
    cat("The p-value (chi-square):          ",signif(p.val,digits=4),"\n\n")
    num.sim <- input$num.sim0
    obs     <- dmultinom(d,prob=p)
    X       <- rmultinom(num.sim,k,p)
    a       <- apply(X,2,function(x,p){dmultinom(x,prob=p)},p=p)
    cat("The p-value (multinomial):         ",signif(sum(a<=obs)/num.sim,digits=4),"\n")
    f <- function(p,k,y.obs,q)
    {
      pbinom(y.obs,size=k,prob=p) - q
    }
    p.lo <- uniroot(f,interval=c(0,1),k=num.sim,y.obs=sum(a<=obs)-1,q=0.975)$root
    if ( sum(a<=obs) == num.sim ) {
      p.hi = 1
    } else {
      p.hi <- uniroot(f,interval=c(0,1),k=num.sim,y.obs=sum(a<=obs),q=0.025)$root
    }
    cat("    95% CI:                       [",signif(p.lo,digits=4),",",signif(p.hi,digits=4),"]\n")
  })
  
  output$oid1 <- renderPrint({
    k <- input$trials1
    m <- input$faces1
    p <- rep(1/m,m)
    set.seed(input$seed1)
    d <- rmultinom(1,k,p)
    cat("The total number of tosses         ",k,"\n")
    cat("The observed number for each face: ",d,"\n\n")
    W <- sum((d-k*p)^2/(k*p))
    p.val <- 1-pchisq(W,length(d)-1)
    cat("The chi-square value:              ",round(W,3),"\n")
    cat("The p-value (chi-square):          ",signif(p.val,digits=4),"\n\n")
    num.sim <- input$num.sim1
    obs     <- dmultinom(d,prob=p)
    X       <- rmultinom(num.sim,k,p)
    a       <- apply(X,2,function(x,p){dmultinom(x,prob=p)},p=p)
    cat("The p-value (multinomial):         ",signif(sum(a<=obs)/num.sim,digits=4),"\n")
    f <- function(p,k,y.obs,q)
    {
      pbinom(y.obs,size=k,prob=p) - q
    }
    p.lo <- uniroot(f,interval=c(0,1),k=num.sim,y.obs=sum(a<=obs)-1,q=0.975)$root
    if ( sum(a<=obs) == num.sim ) {
      p.hi = 1
    } else {
      p.hi <- uniroot(f,interval=c(0,1),k=num.sim,y.obs=sum(a<=obs),q=0.025)$root
    }
    cat("    95% CI:                       [",signif(p.lo,digits=4),",",signif(p.hi,digits=4),"]\n")
  })

  output$oid2 <- renderPrint({
    k <- input$trials2
    m <- input$faces2
    p <- rep(1/m,m)
    set.seed(input$seed2)

    num.rep <- input$num.rep2
    p.chi   <- rep(NA,num.rep)
    p.mult  <- rep(NA,num.rep)
    delta.p <- rep(NA,num.rep)

    num.sim <- input$num.sim2
    for ( ii in 1:num.rep ) {
      X.obs       <- rmultinom(1,k,p)
      W           <- sum((X.obs-k*p)^2/(k*p))
      p.chi[ii]   <- 1 - pchisq(W,m-1)
      obs         <- dmultinom(X.obs,prob=p)
      X           <- rmultinom(num.sim,k,p)
      a           <- apply(X,2,function(x,p){dmultinom(x,prob=p)},p=p)
      p.mult[ii]  <- sum(a<=obs)/num.sim
      delta.p[ii] <- p.mult[ii] - p.chi[ii]
    }

    cat("Mean Difference:",round(mean(delta.p),5),"\n")
    cat("SE Difference:  ",round(sd(delta.p)/sqrt(num.rep),5),"\n\n")
    
    f <- function(mu,s,y.obs,q)
    {
      pnorm(y.obs,mean=mu,sd=s) - q
    }
    mu.lo <- uniroot(f,interval=c(-1,1),y.obs=mean(delta.p),
                     s=sd(delta.p)/sqrt(num.rep),q=0.975)$root
    mu.hi <- uniroot(f,interval=c(-1,1),y.obs=mean(delta.p),
                     s=sd(delta.p)/sqrt(num.rep),q=0.025)$root
    if ( mu.lo > 0 ) {
      cat("Multinomial simulations are more conservative.\n\n")
    } else if ( mu.hi < 0 ) {
      cat("Chi-square is more conservative.\n\n")
    } else {
      cat("Multinomial simulations and chi-square\ngive consistent results.\n\n")
    }
    
  })

  output$plot2 <- renderPlot({
    k <- input$trials2
    m <- input$faces2
    p <- rep(1/m,m)
    set.seed(input$seed2)

    num.rep <- input$num.rep2
    p.chi   <- rep(NA,num.rep)
    p.mult  <- rep(NA,num.rep)
    delta.p <- rep(NA,num.rep)

    num.sim <- input$num.sim2
    for ( ii in 1:num.rep ) {
      X.obs       <- rmultinom(1,k,p)
      W           <- sum((X.obs-k*p)^2/(k*p))
      p.chi[ii]   <- 1 - pchisq(W,m-1)
      obs         <- dmultinom(X.obs,prob=p)
      X           <- rmultinom(num.sim,k,p)
      a           <- apply(X,2,function(x,p){dmultinom(x,prob=p)},p=p)
      p.mult[ii]  <- sum(a<=obs)/num.sim
      delta.p[ii] <- p.mult[ii] - p.chi[ii]
    }

    ggplot(data=data.frame(delta.p),mapping=aes(x=delta.p,y=after_stat(density))) +
      geom_histogram(fill="dodgerblue",bins=25) +
      geom_vline(xintercept=0,col="seagreen",lty=2) +
      xlab(expression(Delta*"p"))
  })
}

shinyApp(ui,server)
