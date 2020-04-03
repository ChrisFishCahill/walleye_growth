library(purrr)
library(ggplot2)
library(tidyr)
library(dplyr)

paths <- dir(pattern = ".rds")



# which failed?
buggered <- out %>% dplyr::filter(convergence == 1)
print(buggered, n = Inf)
whichSims <- buggered$iter
message(paste0(
  paste0(length(whichSims), " out of "),
  length(unique(totest$iter)), " iterations did not converge"
))

table(buggered$sig_varies)

out %>%
  dplyr::filter(!(iter %in% whichSims)) %>%
  dplyr::mutate(sig_varies = paste0("Sim = ", sig_varies)) %>%
  dplyr::mutate(sig_varies_fitted = paste0("Fitted = ", sig_varies_fitted)) %>%
  ggplot(aes(ln_global_omega)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = out[["true_ln_global_omega"]][1]) +
  facet_grid(sig_varies_fitted ~ sig_varies) +
  xlab(expression(omega)) +
  ylab("Count")
ggsave("sim2/hist-sim.pdf", width = 7, height = 5)

out$sig_varies_fitted <- factor(out$sig_varies_fitted, levels = c("by lake", "by time", "both", "ar1 st"))
out$sig_varies <- factor(out$sig_varies, levels = c("by lake", "by time", "both", "ar1 st"))

out %>%
  dplyr::filter(!(iter %in% whichSims)) %>%
  dplyr::mutate(Matched = ifelse(sig_varies_fitted == sig_varies, TRUE, FALSE)) %>%
  ggplot(aes(sig_varies, exp(ln_global_omega), colour = sig_varies_fitted, fill = Matched)) +
  geom_boxplot() +
  geom_hline(yintercept = exp(out[["true_ln_global_omega"]][1])) +
  xlab(expression(Simulated ~ omega ~ variation)) +
  labs(colour = expression(Fitted ~ omega ~ variation)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_manual(values = c("white", "grey60")) +
  ylab(expression(omega[0]))
ggsave("sim2/boxplot-sim.pdf", width = 7, height = 5)
