library(tidyverse)
library(cowplot)
set.seed(12345)
s1 <- factor(c(rbinom(100, 20, .3), rbinom(13,20,.8), rep(19,7)), levels = 1:20)
s2 <- factor(c(rbinom(100, 20, .6), rbinom(20,20,.8)), levels = 1:20)

raw <- tibble(
  count = c(table(s1), table(s2)),
  x = c(1:20 -.25, 1:20 +.25),
  Sample = rep(c("One","Two"), each=20))
A <- ggplot(raw, aes(x, count, fill = Sample)) +
  geom_col() +
  scale_fill_viridis_d(end = .5, drop=FALSE) +
  theme_bw() +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = c(1,1),
        legend.justification = c("right","top"),
        legend.background = element_blank())

mover <- raw %>%
  mutate(x = round(x)) %>%
  pivot_wider(names_from = Sample, values_from = count) %>%
  mutate(removed = pmax(One - Two, 0), added = pmax(Two - One, 0)) %>%
  mutate(unchanged = pmax(One - removed, 0)) %>%
  select(-Two,-One) %>%
  pivot_longer(-x) %>%
  mutate(name = fct_relevel(name, "removed", "added", "unchanged")) %>%
  filter(x < 20)

our_green <- viridis::viridis_pal(begin = .5, end = .5)(1)
our_purple <- viridis::viridis_pal(begin = 0, end = 0)(1)
fills <- c("removed" = NA, "added" = our_purple, "unchanged" = our_purple)
alphas <- c("removed" = 0, "added" = .3, "unchanged" = 1)
B <- ggplot(mover, aes(x, value, fill=name, alpha=name)) +
  geom_col(color = our_purple, width = .5) +
  scale_fill_manual(values = fills) +
  scale_alpha_manual(values = alphas) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  xlab("") +
  ylab("count") +
  theme(legend.position = c(1,1),
        legend.background = element_blank(),
        legend.justification = c("right","top"),
        legend.title = element_blank())

ecdfs <- raw %>%
  mutate(x = round(x)) %>%
  group_by(Sample) %>%
  mutate(ecdf = cumsum(count),
         ecdf = ecdf / max(ecdf),
         count = NULL)
ecdf_area <- ecdfs %>%
  pivot_wider(names_from = Sample, values_from = ecdf)
ecdf_area <- bind_rows(
  old = ecdf_area,
  new = ecdf_area %>% mutate(across(-x, lag)),
  .id = "source") %>%
  arrange(x, source) %>%
  mutate(top = ifelse(Two >= One, "Two","One")) %>%
  filter(!is.na(top))
C <- ggplot(ecdfs, aes(x)) +
  geom_step(aes(y = ecdf, color = Sample)) +
  scale_color_viridis_d(end = .5, drop = FALSE) +
  theme_bw() +
  xlab("") +
  theme(legend.position = c(1, 0),
        legend.justification = c("right","bottom"),
        legend.background = element_blank()) +
  geom_ribbon(data = ecdf_area, aes(ymin = Two, ymax = One, fill = top),
              alpha = .3) +
  scale_fill_viridis_d(end = .5) +
  annotate("text", x = c(8.25,16.5), y = c(.5, .7), label = c("+","-"),
           color = c(our_purple, our_green), size=10) +
  annotate("curve", x = 17, y = .715, xend = 18.75, yend = .9, curvature = .3,
           color = our_green, size = 1, lineend = "round",
           arrow = arrow()) +
  guides(color = guide_legend("Sample"), fill = "none")

all_plots <- plot_grid(A, B, C, align = "v", labels = LETTERS[1:3], nrow=1)

Rcpp::sourceCpp("minimal-shiftscore.cpp")

ems <- ShiftScore(Matrix::drop0(as.vector(s1)), Matrix::drop0(as.vector(s2)),
                  nresamp = 1000)

ggsave2("emd-display.pdf", width = 12, height = 4)
