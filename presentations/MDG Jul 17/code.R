library(ggplot2)

# Plot of power against ICC

pow_ICC <- function(rho)
{
  # n_i <- 175
  # k=26 for rho=0.05, k=37 for rho=0.1222314 
  m <- 10
  k <- 26
  s <- 1
  p <- pnorm( sqrt((m*k*0.3^2)/(2*s^2*(1+(m-1)*rho))) - 1.96)
  return(p)
}

df <- data.frame(rho=seq(0,1,0.005))
df$p <- sapply(df$rho, pow_ICC)

# dataframe for necessary n
df_n  <- data.frame(rho=seq(0,1,0.005))
df_n$n <- 175*(1+(10-1)*df_n$rho)
df_n$p <- (df_n$n-min(df_n$n))/(max(df_n$n) - min(df_n$n))

ggplot(df, aes(rho, p)) + geom_line() + geom_vline(xintercept = 0.05, linetype=2) + ylab("Power")
 # geom_line(data=df_n)
ggsave("U:/Projects/MRC SDF/WP2/Notes/MDG pres Jul 17/power.pdf", height=4)

# Power vs n for two different values of rho
pow_n <- function(x)
{
  n <- x[1]
  rho <- x[2]
  m <- 10
  s <- 1
  p <- pnorm( sqrt((n*0.3^2)/(2*s^2*(1+(m-1)*rho))) - 1.96)
  return(p)
}

df3 <- data.frame(n = rep(seq(1, 500), 2), rho=c(rep(0.01, 500), rep(0.05, 500)))
df3$p <- apply(df3, 1, pow_n)
df3$rho <- as.factor(df3$rho)
ggplot(df3, aes(n, p, colour=rho)) + geom_line() + ylab("Power") + scale_y_continuous(breaks = seq(0,1,0.2)) +
  geom_vline(xintercept = 219, linetype=2)
ggsave("U:/Projects/MRC SDF/WP2/Notes/MDG pres Jul 17/power_2ICCs.pdf", height=4)

# Prior distrbution for the ICC
m <- 0.05
n <- 60
df2 <- data.frame(rho = rbeta(1000000, n*m+1, n*(1-m)+1))
# a = n*m+1 = 4
# b = n*(1-m)+1 = 58
# rho_u = qbeta(0.95, 4, 58) = 0.1222314 
ggplot(df2, aes(rho)) + geom_histogram(aes(y=..ncount..), bins=100, alpha=0.2, colour="black") +
  geom_line(data=df[df$rho <= max(df2$rho), ], aes(rho, p), colour="red") +
  scale_y_continuous(breaks = seq(0,1,0.2)) + ylab("Power")
ggsave("U:/Projects/MRC SDF/WP2/Notes/MDG pres Jul 17/prior_ICC.pdf", height=4)

df2$p <- sapply(df2$rho, pow_ICC)
df2$t <- df2$p < 0.8
quantile(df2$p, c(0.05))
ggplot(df2, aes(p)) + geom_histogram(aes(y=..ncount..), bins=100, alpha=0.2, colour="black") +
  scale_y_continuous(breaks = seq(0,1,0.2)) + ylab("")
ggsave("U:/Projects/MRC SDF/WP2/Notes/MDG pres Jul 17/prior_power.pdf", height=4)

# Prior for a treatment effect
power_known_s <- function(mu, alpha, n, s)
{
  prob <- 1-pnorm(qnorm(1-alpha)-mu/sqrt(2*(s^2)/n))
  return(prob)
}

df4 <- data.frame(mu = rnorm(10000, 0.2, 0.25))
df4$p <- apply(df4, 1, power_known_s, alpha=0.025, n=500, s=1)
mean(df4$p)
ggplot(df4, aes(p)) + geom_histogram(aes(y=..ncount..), bins=100, alpha=0.2, colour="black") +
  scale_y_continuous(breaks = seq(0,1,0.2))

# Expected power for different n
assur <- function(x)
{
  alpha <- 0.025
  n <- x[1]
  m <- 0.2
  nu <- 0.25^2
  tau <- sqrt((2*(1)^2)/n)
  z <- qnorm(1-alpha)
  return(pnorm((-tau*z+m)/(sqrt(tau^2 + nu))))
}

df5 <- data.frame(n=seq(1, 700))
df5$p <- apply(df5, 1, assur)
df5$y <- df5$p/df5$n
ggplot(df5, aes(n, p)) + geom_line() + ylab("Unconditional power")
ggsave("U:/Projects/MRC SDF/WP2/Notes/MDG pres Jul 17/stallard.pdf", height=4)

# MIDSHIIPS feasibility aspects example
# Pilot trial sample of 30 in each arm is designed to give a pooled SD estimate of
# sd +/- 0.39%, follow up +/- 7%, adherance +/- 8 to 18%
# Note that it is assumed recruitment rates will be precisely estimated
# Show what the assurance metrics would be
