library("tidyverse")
library("patchwork")
library("geomtextpath")

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = 0.27, 
                      label = as.character(0.27)),
                  colour = "#0072B2", size = 3, linewidth = 1) +
  labs(y = "Predation", x = "Individual fecundity", 
       title = "Jones & Comita 2010") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        plot.title = element_text(face = "italic")) -> predi_jc

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = 0.26, 
                      label = as.character(0.26)),
                  colour = "#0072B2", size = 3, linewidth = 1) +
  labs(y = "Predation", x = "Individual fecundity", title = "Jackson et al. 2023") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        plot.title = element_text(face = "italic")) -> predi_us

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = 0.21, 
                      label = as.character(0.21)),
                  colour = "#0072B2", size = 3, linewidth = 1) +
  labs(y = "Predation", x = "Neighbourhood fecundity") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm")))) -> predn_jc

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = -0.13, 
                      label = as.character(-0.13)),
                  colour = "#D55E00", size = 3, linewidth = 1) +
  labs(y = "Predation", x = "Connectivity") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm")))) -> predn_us

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = 0.11, 
                      label = as.character(0.11)),
                  colour = "#0072B2", size = 3, linewidth = 1) +
  labs(y = "Individual fecundity", x = "Neighbourhood fecundity") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm")))) -> set_jc

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = -0.12, 
                      label = as.character(-0.12)),
                  size = 3, linewidth = 1) +
  labs(y = "Individual fecundity", x = "Connectivity") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm")))) -> set_us

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = 0.10, 
                      label = as.character("0.10")),
                  colour = "#0072B2", size = 3, linewidth = 1) +
  labs(y = "Realised fecundity", x = "Neighbourhood fecundity") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm")))) -> real_jc

ggplot() +
  geom_textabline(aes(intercept = 0.3, slope = -0.06, 
                      label = as.character(-0.06)),
                  size = 3, linewidth = 1) +
  labs(y = "Realised fecundity", x = "Connectivity") +
  theme_classic(base_size = 7) +
  theme(axis.line.y = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm"))),
        axis.line.x = element_line(arrow = grid::arrow(ends = "last", length = unit(2, "mm")))) -> real_us

(predi_jc + predi_us) /
  (predn_jc + predn_us) /
  (set_jc + set_us) /
  (real_jc + real_us) + plot_annotation(tag_levels = 'a')

ggsave(here::here("output", "figures", "comparison_schematic.png"),
                  dpi = 600, width = 80, height = 150, units = "mm")
