library(ggplot2)
library(wesanderson)
library(latex2exp)

df <- read.csv("results.csv")

df$K <- as.factor(df$K)
df$sampler <- as.factor(df$sampler)
levels(df$sampler)
df$ess_per_time <- df$load_min_ess / df$time
df2 <- aggregate(.~ N + K + P + sampler, data=df, mean)

N_labs <- c("N = 50", "N = 100", "N = 500", "N = 1000")
# N_labs <- c(TeX("$N$ = 50"), TeX("$N$ = 100"), TeX("$N$ = 500"), TeX("$N$ = 1000"))
names(N_labs) <- c("50", "100", "500", "1000")
P_labs <- c("P = 10", "P = 100", "P = 1000")
names(P_labs) <- c("10", "100", "1000")


plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

df2$log_ess_per_time <- sapply(df2$ess_per_time, log10)

plt <- ggplot(df2, aes(x=K, y=ess_per_time, color=sampler, group=sampler)) + geom_point(size=0.5) + geom_line() +
  scale_y_continuous(trans="log10", labels=plain, breaks=c(0.0001, 0.01, 1, 100)) +
  facet_grid(col=vars(N), row=vars(P), scales="fixed",
             labeller = labeller(N=N_labs, P=P_labs)) +
  labs(y=expression("ESS/second")) +
  scale_color_manual(name="sampler", breaks=c("old", "gibbs", "hmc", "ortho"), labels=c("conditional Gibbs", "joint Gibbs", "HMC", "orthogonal"), values=wes_palette("Darjeeling2")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA))
# plt <- ggplot(df) + geom_histogram(aes(x=ess_per_time))
plt

ggsave("timing.pdf", width=5.5, height=4, units="in")