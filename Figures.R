library(raster)
library(tidyverse)
library(tmap)
library(lubridate)
library(ggplot2)
library(cowplot)
library(grid)
library(treemapify)
library(ggsci)
library(corrplot)
library(ggrepel)
library(mgcv)

source("Functions.R")

########################## Input data #################################
#=========================== Sichuan map ==============================
SC.pref <- shapefile(".\\Data\\Sichuan_City.shp")
SC.pref$name <- c("Chengdu", "Zigong","Panzhihua","Luzhou","Deyang","Mianyang","Guangyuan","Suining","Neijiang","Leshan","Nanchong","Meishan","Yibin","Guang'an","Dazhou","Yaan","Bazhong","Ziyang","Aba","Ganzi","Liangshan")   # add prefecture names


#=========================== HFMD data ==============================
dat <- read_csv(".\\Data\\Sichuan_pref_yearly.csv") %>% as.data.frame() 

dat.6 <- dat %>%                                       # filter to data before 2015
  filter(Year < 2015)







#####################################################################
#
#                           Main text
#
#####################################################################

########################## Figure 2 #################################
#======== Figure 2A ===========
# aggregate data by year and severity
dat.agg.year <- dat %>%
  group_by(Year) %>% 
  summarize(Y_severe = sum(Y_severe), Y_mild = sum(Y_mild), population = sum(pop)) %>%
  mutate(Severe = Y_severe/population*100000, Mild = Y_mild/population*100000) %>%
  gather(Type, inc, Severe:Mild) %>%
  mutate(Type = factor(Type, levels = c("Severe", "Mild")))


F2.A <- ggplot(dat.agg.year, aes(x = Year, y = inc, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(x = "Year", y = "Incidence rate (1/100,000)", fill = "") +
  theme(legend.position = c(0.0, 0.95), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))


#======== Figure 2B ===========
# incidence rate map
dat.inc.mean <- dat %>% 
  mutate(inc = Y/pop*100000) %>%
  group_by(Pref) %>%
  summarize(inc = mean(inc))

SC.pref@data$inc <- dat.inc.mean$inc

tm_shape(SC.pref) + 
  tm_polygons("inc", n = 6, style = "pretty", palette = "OrRd", title = "Annual mean\nincidence rate\n(per 100,000)") +
  tm_text("name", size = 0.7) +
  tm_layout(legend.position = c("left", "bottom"), legend.text.size = 0.7, main.title.size = 0.7, inner.margins = c(0,0.1,0,0), main.title.fontface = "bold") +
  tm_layout(legend.title.size = 0.8, frame = FALSE)
F2.B <- grid.grab()

#======== Figure 2C ===========
# test count ts
test.agg.year <- dat %>%
  mutate(Mild = Z_mild_1 + Z_mild_2 + Z_mild_3, Severe =  Z_severe_1 + Z_severe_2 + Z_severe_3) %>%
  gather(Type, Test, Mild:Severe) %>%
  mutate(Type = factor(Type, levels = c("Severe", "Mild")))

F2.C <- ggplot(test.agg.year, aes(x = Year, y = Test, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(x = "Year", y = "Number of lab-typed samples", fill = "") +
  theme(legend.position = c(0.0, 0.95), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))

#======== Figure 2D ===========
# test count map
dat.test.count <- dat %>%
  mutate(test.total = Z_mild_1 + Z_mild_2 + Z_mild_3 + Z_severe_1 + Z_severe_2 + Z_severe_3) %>%
  group_by(Pref) %>%
  summarize(test.total = sum(test.total)) %>%
  mutate(test.prop = test.total/sum(test.total))
SC.pref$test.prop <- dat.test.count$test.prop

tm_shape(SC.pref) + 
  tm_polygons("test.prop", n = 7, style = "pretty", palette = "Greens", title = "Proportion of\nall typed\ncases") +
  tm_text("name", size = 0.7) +
  tm_layout(legend.position = c("left", "bottom"), legend.text.size = 0.7, main.title.size = 0.8, inner.margins = c(0,0.05,0,0), main.title.fontface = "bold") +
  tm_layout(legend.title.size = 0.8, frame = FALSE)
F2.D <- grid.grab()

# save_plot("Result\\Figure\\Fig 2.pdf", plot_grid(F2.A, F2.B, F2.C, F2.D, ncol = 2, labels = c("(A)", "(B)", "(C)", "(D)"), rel_widths = c(1, 1.3), hjust = 0, label_x = 0.02), base_width = 8, base_height = 8)








########################## Figure 3 #################################
# get the archetypal designs existing and IncRate
load("Result\\Disease model\\Disease_model_realizations_6.RData")    # output of Disease system model realizations.R

# process data to get the data for alternative designs
dat.agg.alter <- dat.6 %>%
  mutate(Z_mild = Z_mild_1 + Z_mild_2 + Z_mild_3, Z_severe = Z_severe_1 + Z_severe_2 + Z_severe_3, inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Year, Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop)) %>%
  mutate(inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop), inc = mean(inc), inc_severe = mean(inc_severe)) %>%
  mutate(Z = Z_mild + Z_severe)

# Append the data to the shapefile
SC.pref$existing <- dat.agg.alter$Z/sum(dat.agg.alter$Z)
SC.pref$IncRate <- dat.agg.alter$inc/sum(dat.agg.alter$inc)
SC.pref$inc <- dat.agg.alter$inc

# get the optmization result
load("Result\\Optimization\\GA_all_TestX1.Rdata")
SC.pref$optimal.all <- convert.var(n = 21, x = GA.All.X1@solution[1:20])

load("Result\\Optimization\\GA_severe_TestX1.Rdata")
SC.pref$optimal.severe <- convert.var(n = 21, x = GA.severe.X1@solution[1:20])

Fig3.dat <- SC.pref@data %>%
  arrange(-inc) %>%
  dplyr::select(name, inc, existing, optimal.all, optimal.severe, IncRate) %>%
  gather("Type", "design", existing:IncRate) %>%
  mutate(Type = factor(Type, levels = c("existing", "IncRate","optimal.all", "optimal.severe"))) %>%
  mutate(Type = factor(Type, labels = c("(A) Existing","(B) IncRate", "(C) Optimal for all", "(D) Optimal for severe")))

# treemaps
Fig3.abcd <- Fig3.dat %>%
  ggplot(aes(area = design, fill = inc, label = name)) +
  geom_treemap(layout = "fixed", start = "topleft", col = "white", size = 1.5) +
  geom_treemap_text(layout = "fixed", start = "topleft", grow = FALSE, reflow = TRUE, place = "center", size = 10, min.size = 1) +
  facet_wrap(~Type) +
  scale_fill_material("deep-orange") +
  labs(fill = "Annual mean incidence rate (1/100,000)") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = "bold", margin = margin(0,0,5,0)), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.key.width = unit(1.5, "cm"), 
        legend.title = element_text(size = 12, vjust = 1), 
        legend.text = element_text(size = 10), 
        legend.key.height = unit(0.3, "cm"), 
        panel.border = element_rect(size = 3, color = "red"),
        panel.spacing.y = unit(1, "lines"),
        legend.justification = "center", 
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0)
  )


# scatterplots
Fig3e <- SC.pref@data %>%
  filter(!name %in% c("Chengdu", "Meishan")) %>%
  ggplot(aes(x = inc, y = optimal.all, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = IncRate, yend = optimal.all, col = as.factor(IncRate - optimal.all > 0)), size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "royalblue1", shape = 17) +
  geom_text_repel(size = 3, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = IncRate), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(E) IncRate and Optimal for all")+
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0, unit = "cm")) 

Fig3e.inset <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.all, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = IncRate, yend = optimal.all, col = as.factor(IncRate - optimal.all > 0)), size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 1, col = "royalblue1", shape = 17) +
  geom_text_repel(data = SC.pref@data[SC.pref$name %in% c("Chengdu", "Meishan"), ], aes(x = inc, y = optimal.all, label = name), size = 3, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = IncRate), size = 1) +
  guides(col = FALSE) +
  xlab("") +
  ylab("") +
  annotate(geom = "rect", ymax = 0.1, ymin = 0.01, xmin = 15, xmax = 100, fill = NA, col = "red", linetype = 2) +
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 8, hjust = 0.5), axis.title = element_text(size = 12), axis.text.x = element_text(size = 8), plot.margin = margin(t = 0, unit = "cm"))

Fig3e.all <- ggdraw() +
  draw_plot(Fig3e) +
  draw_plot(Fig3e.inset, x = 0.07, y = .55, width = .5, height = .4)

Fig3f <- SC.pref@data %>%
  filter(!name %in% c("Chengdu", "Meishan")) %>%
  ggplot(aes(x = inc, y = optimal.severe, label = name)) +
  geom_smooth(method = "lm") + 
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_segment(aes(x = inc, xend = inc, y = IncRate, yend = optimal.severe, col = as.factor(IncRate - optimal.severe > 0)), size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  geom_point(size = 2, col = "blue3", shape = 15) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = IncRate), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(F) IncRate and Optimal for severe")+
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0, unit = "cm")) +
  ylim(0.015, 0.1)

Fig3f.inset <- SC.pref@data %>%
  #  ggplot(aes(x = inc, y = GA.inc.optimal.pest, label = name)) +
  ggplot(aes(x = inc, y = optimal.severe, label = name)) +
  geom_smooth(method = "lm") + 
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_segment(aes(x = inc, xend = inc, y = IncRate, yend = optimal.severe, col = as.factor(IncRate - optimal.severe > 0)), size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  geom_point(size = 1, col = "blue3", shape = 15) +
  geom_text_repel(data = SC.pref@data[SC.pref$name %in% c("Chengdu", "Meishan"), ], size = 3, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = IncRate), size = 1) +
  guides(col = FALSE) +
  annotate(geom = "rect", ymax = 0.1, ymin = 0.015, xmin = 15, xmax = 100, fill = NA, col = "red", linetype = 2) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 8, hjust = 0.5), axis.title = element_text(size = 12), axis.text.x = element_text(size = 8), plot.margin = margin(t = 0, unit = "cm"))

Fig3f.all <- ggdraw() +
  draw_plot(Fig3f) +
  draw_plot(Fig3f.inset, x = 0.07, y = .55, width = .5, height = .4)

Fig3.row2 <- plot_grid(Fig3e.all, Fig3f.all)

# save_plot("Result\\Figure\\Fig 3.pdf", plot_grid(Fig3.abcd, Fig3.row2, ncol = 1, rel_heights = c(8, 4)), base_width = 8, base_height = 11)









########################## Figure 4 #################################
# for both all and severe, with finer scale
load("Result\\Optimal_psevere_optimal_designs.Rdata")    # output of Optimal psevere optimal designs.R

result.all$ofv <- -result.all$ofv

result.all.smooth <- result.all %>%
  group_by(goal) %>%
  do(data.frame(., pred = predict(gam(ofv ~ s(p_severe, bs = "gp"), data = .))))

result.all.min.smooth <- result.all %>%
  group_by(goal) %>%
  do(data.frame(., pred = predict(gam(ofv ~ s(p_severe, bs = "gp"), data = .)))) %>%
  summarize(ofv.min = min(pred), p_severe = p_severe[which.min(pred)])

result.all.GA.optimal <- result.all.min.smooth %>%
  mutate(p_severe = c(0.17, 0.70), ofv.min = c(1.978, 0.02733), goal = c("(A) For all", "(B) For severe"))

levels(result.all.smooth$goal) <- levels(result.all.min.smooth$goal) <- c("(A) For all", "(B) For severe")

result.all.smooth %>%
  ggplot(aes(x = p_severe, y = ofv, group = goal, col = goal)) +
  geom_point(aes(shape = goal)) +
  geom_line(aes(x = p_severe, y = pred, group = goal, col = goal)) +
  facet_wrap(~goal, scales = "free") +
  guides(col = FALSE, shape = FALSE) +
  xlab("Proportion of severe cases subjected to serotyping") +
  ylab("MAE (1/100,000)") +
  geom_point(data = result.all.min.smooth, aes(x = p_severe, y = ofv.min, shape = goal), col = "black", size = 2.5) +
  geom_point(data = result.all.GA.optimal, aes(x = p_severe, y = ofv.min, shape = goal), col = "blue", size = 2.5) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold", size = 14), axis.text.y = element_text(angle = 90, hjust = 0.5))

# ggsave("Result\\Figure\\Fig 4.pdf", width = 8, height = 4)







########################## Figure 5 #################################
load("Result\\Compare_OFV_of_benchmark_all_09to14.Rdata")   # output of Compare OFV of archetypal designs.R
result.all.inc <- result.all %>%
  mutate(Type = factor(Type, levels = unique(Type)), ofv = -ofv) 
levels(result.all.inc$Type) <- c("Optimal", "Existing", "Equal", "PopSize", "Case", "IncRate", "SevereCase", "SevereIncRate")

Fig5A <- result.all.inc %>%
  ggplot(aes(x = Type, y = ofv, fill = Type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, notch = TRUE) +
  xlab("") + 
  ylab("MAE (1 per 100,000)") +
  labs(fill = "")+
  ggtitle("(A) For all, 2009-2014") +
  guides(fill = FALSE) +
  theme(legend.position = c(0.75, 0.9), plot.title = element_text(hjust = 0, size = 12), axis.text.x = element_text(size = 10, face = "italic", angle = 30, hjust = 1), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 10), legend.text = element_text(size = 10, face = "italic"), axis.title.x  = element_blank()) +
  scale_fill_d3()



load("Result\\Compare_OFV_of_benchmark_severe_09to14.Rdata")    # output of Compare OFV of archetypal designs.R
result.all.inc.severe <- result.all %>%
  mutate(Type = factor(Type, levels = unique(Type)), ofv = -ofv) 
levels(result.all.inc.severe$Type) <- c("Optimal", "Existing", "Equal", "PopSize", "Case", "IncRate", "SevereCase", "SevereIncRate")

Fig5B <- result.all.inc.severe %>%
  ggplot(aes(x = Type, y = ofv, fill = Type)) +
  geom_violin()  +
  geom_boxplot(width = 0.3, notch = TRUE) +
  xlab("") + 
  ylab("MAE (1 per 100,000)") +
  labs(fill = "")+
  ggtitle("(B) For severe, 2009-2014") +
  guides(fill = FALSE) +
  theme(legend.position = c(0.75, 0.9), plot.title = element_text(hjust = 0, size = 12), axis.text.x = element_text(size = 10, face = "italic", angle = 30, hjust = 1), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 10), legend.text = element_text(size = 10, face = "italic"), axis.title.x  = element_blank()) +
  scale_fill_d3()



load("Result\\Compare_OFV_of_benchmark_all_15.Rdata")     # output of Compare OFV of archetypal designs.R
result.all.inc.2015 <- result.all %>%
  mutate(Type = factor(Type, levels = unique(Type))) 

levels(result.all.inc.2015$Type) <- c("Optimal", "Existing", "Equal", "PopSize", "Case", "IncRate", "SevereCase", "SevereIncRate")

Fig5C <- result.all.inc.2015 %>%
  ggplot(aes(x = Type, y = ofv, fill = Type)) +
  geom_violin() +
  geom_boxplot(width = 0.3, notch = TRUE) +
  xlab("") + 
  ylab("MAE (1 per 100,000)") +
  labs(fill = "")+
  ggtitle("(C) For all, 2015") +
  guides(fill = FALSE) +
  theme(legend.position = c(0.75, 0.9), plot.title = element_text(hjust = 0, size = 12), axis.text.x = element_text(size = 10, face = "italic", angle = 30, hjust = 1), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 10), legend.text = element_text(size = 10, face = "italic"), axis.title.x  = element_blank()) +
  scale_fill_d3()



load("Result\\Compare_OFV_of_benchmark_severe_15.Rdata")     # output of Compare OFV of archetypal designs.R
result.all.inc.severe.2015 <- result.all %>%
  mutate(Type = factor(Type, levels = unique(Type))) 

levels(result.all.inc.severe.2015$Type) <- c("Optimal", "Existing", "Equal", "PopSize", "Case", "IncRate", "SevereCase", "SevereIncRate")

Fig5D <- result.all.inc.severe.2015 %>%
  ggplot(aes(x = Type, y = ofv, fill = Type)) +
  geom_violin()  +
  geom_boxplot(width = 0.3, notch = TRUE) +
  xlab("") + 
  ylab("MAE (1 per 100,000)") +
  labs(fill = "")+
  ggtitle("(D) For severe, 2015") +
  guides(fill = FALSE) +
  theme(legend.position = c(0.75, 0.9), plot.title = element_text(hjust = 0, size = 12), axis.text.x = element_text(size = 10, face = "italic", angle = 30, hjust = 1), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 10), legend.text = element_text(size = 10, face = "italic"), axis.title.x  = element_blank()) +
  scale_fill_d3()

# save_plot("Result\\Figure\\Fig 5.pdf", plot_grid(Fig5A, Fig5B, Fig5C, Fig5D), base_width = 8, base_height = 8)







########################## Figure 6 #################################
load("Result\\Optimization\\GA_all_TestX0.5.Rdata")
load("Result\\Optimization\\GA_all_TestX2.Rdata")
load("Result\\Optimization\\GA_all_TestX5.Rdata")
load("Result\\Optimization\\GA_severe_TestX0.5.Rdata")
load("Result\\Optimization\\GA_severe_TestX2.Rdata")
load("Result\\Optimization\\GA_severe_TestX5.Rdata")

SC.pref$optimal.all.half <- convert.var(n = 21, x = GA.All.X0.5@solution[1:20])
SC.pref$optimal.all.2 <- convert.var(n = 21, x = GA.All.X2@solution[1:20])
SC.pref$optimal.all.5 <- convert.var(n = 21, x = GA.All.X5@solution[1:20])

SC.pref$optimal.severe.half <- convert.var(n = 21, x = GA.severe.X0.5@solution[1:20])
SC.pref$optimal.severe.2 <- convert.var(n = 21, x = GA.severe.X2@solution[1:20])
SC.pref$optimal.severe.5 <- convert.var(n = 21, x = GA.severe.X5@solution[1:20])



SC.copy.resource <- SC.pref@data %>%
  arrange(-inc) %>%
  dplyr::select(name, inc, optimal.all.half, optimal.all.2, optimal.all.5, optimal.severe.half, optimal.severe.2, optimal.severe.5) %>%
  gather("Type", "design", optimal.all.half:optimal.severe.5) %>%
  mutate(Type = factor(Type, levels = c("optimal.all.half", "optimal.all.2", "optimal.all.5", "optimal.severe.half", "optimal.severe.2", "optimal.severe.5")))

SC.copy.resource$Type <- factor(SC.copy.resource$Type, labels = c("(A) For all, resource*0.5", "(B) For all, resource*2", "(C) For all, resource*5", "(D) For severe, resource*0.5", "(E) For severe, resource*2", "(F) For severe, resource*5"))

SC.copy.resource %>%
  ggplot(aes(area = design, fill = inc, label = name)) +
  geom_treemap(layout = "fixed", start = "topleft", col = "white", size = 1.5) +
  geom_treemap_text(layout = "fixed", start = "topleft", grow = FALSE, reflow = TRUE, place = "center", size = 10, min.size = 1) +
  facet_wrap(~Type) +
  scale_fill_material("deep-orange") +
  labs(fill = "Annual mean incidence rate (1/100,000)") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = "bold", margin = margin(0,0,5,0)), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.key.width = unit(1.5, "cm"), 
        legend.title = element_text(size = 12, vjust = 1), 
        legend.text = element_text(size = 10), 
        legend.key.height = unit(0.4, "cm"), 
        panel.border = element_rect(size = 3, color = "red"),
        panel.spacing.y = unit(1, "lines"),
        legend.justification = "center", 
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0)
  )
# ggsave("Result\\Figure\\Fig 6.pdf", width = 8, height = 6.5)







########################## Figure 7 #################################
p_severe_pattern <- data.frame(obj = rep(c("(A) For all", "(B) For severe"), each = 4), resource = c("*0.5", "*1", "*2", "*5"), p.severe = c(GA.All.X0.5@solution[21], 0.11, GA.All.X2@solution[21], GA.All.X5@solution[21], GA.severe.X0.5@solution[21], 0.65,  GA.severe.X2@solution[21], GA.severe.X5@solution[21]))

ggplot(p_severe_pattern, aes(x = as.factor(resource), y = p.severe, group = obj)) +
  geom_point() +
  geom_line() +
  facet_wrap(~obj, scales = "free") +
  xlab("Typing resources") +
  ylab("Proportion of severe cases to be serotyped") +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 14, face = "bold", margin = margin(0,0,5,0)), axis.title = element_text(size = 13)) +
  scale_x_discrete(labels = c("*0.5", "Historical", "*2", "*5"))
# ggsave("Result\\Figure\\Fig 7.pdf", width = 8, height = 4)







#####################################################################
#
#                    Supplementary Information
#
#####################################################################

######################## Figure ST1.1 ###############################
load("Result\\Disease model\\Disease_model_realizations_7.RData")  # output of Disease system model realizations.R

dat.agg <- dat %>%
  group_by(Year) %>% 
  summarize(Y = sum(Y), population = sum(pop)) %>%
  mutate(inc.rate = Y/population*100000) 

# province-level serotype-specific incidence rate 4000*7*3
inc.fit <- dm.realization$lambda_ikt
for(i in 1:4000)
{
  for(k in 1:3)
  {
    inc.fit[i,,k,] <- dm.realization$lambda_ikt[i,,k,] * dm.realization$N_it *10000/100000
  }
}

count.fit <- apply(inc.fit, c(1, 3, 4), sum)

result.3 <- NULL
for(i in 1:4000)
{
  current <- data.frame(rep = i, t = 2009:2015, k = rep(1:3, each = 7), inc = as.vector(apply(count.fit[i,,], 1, function(x) x/dat.agg$population*1000000)))
  result.3 <- rbind(result.3, current)
}

random.ID <- sample(1:4000, 100)
ggplot(result.3[result.3$rep %in% random.ID, ], aes(x = t, y = inc, col = as.factor(k), group = interaction(rep, k))) +
  geom_line(alpha = 0.1) +
  labs(col = "") +
  ylab("Incidence rate by serotype (per 100,000)") +
  xlab("Year") +
  scale_color_discrete(name = "", labels = c("CV-A16", "EV-A71", "Other")) +
  theme(legend.position = c(0.05, 0.9), axis.text = element_text(size = 10), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

# ggsave("Result\\Figure\\Fig ST1.1.pdf", width = 7, height = 3.5)





######################## Figure ST1.2 ###############################
inf.fit.df <- data.frame(repi = 1:4000, i = rep(SC.pref$name, each = 4000), k = rep(c("CV-A16", "EV-A71", "Other"), each = 4000*21), t = rep(2009:2015, each = 4000*21*3), inc = as.vector(inc.fit)/10)

inf.fit.df$i <- factor(inf.fit.df$i, levels = SC.pref$name)

ggplot(inf.fit.df[inf.fit.df$repi %in% random.ID, ], aes(x = t, y = inc, group = interaction(repi, k), col = as.factor(k))) +
  geom_line(alpha = 0.1) +
  facet_wrap(~i, scales = "free_y") +
  xlab("Year") +
  ylab("Incidence rate by serotype (1/100,000)") +
  labs(col = "")+
  theme(legend.position = c(0.8, 0.1), axis.text = element_text(size = 8), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.title = element_text(size = 10), legend.text = element_text(size = 10), strip.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

# ggsave("Result\\Figure\\Fig ST1.2.pdf", width = 8, height = 7)





######################## Figure ST1.3 ###############################
p_severe.dat <- data.frame(rep = 1:4000, k = rep(c("CV-A16", "EV-A71", "Other"), each = 4000), p.severe = as.vector(dm.realization$p_severe_k))

ggplot(p_severe.dat, aes(x = p.severe, fill = as.factor(k), group = k)) +
  geom_histogram(bins = 200) +
  #  facet_wrap(~k, scales = "free") +
  #  guides(fill = FALSE) +
  xlab("Probability of severe case") +
  labs(fill = "") +
  theme(legend.position = c(0.8, 0.9), axis.text = element_text(size = 10), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Frequency")

# ggsave("Result\\Figure\\Fig ST1.3.pdf", width = 8, height = 4.5)






######################## Figure S1 ###############################
dat.test <- dat %>%
  mutate(Z_mild = Z_mild_1 + Z_mild_2 + Z_mild_3, Z_severe = Z_severe_1 + Z_severe_2 + Z_severe_3, Z = Z_mild + Z_severe) %>%
  group_by(Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Y_mild = sum(Y_mild), Z = sum(Z), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe)) %>%
  ungroup() %>%
  mutate(test.rate.all = Z/Y*100, test.rate.mild = Z_mild/Y_mild*100, test.rate.severe = Z_severe/Y_severe*100)

SC.pref$test.rate.all <- dat.test$test.rate.all
SC.pref$test.rate.mild <- dat.test$test.rate.mild
SC.pref$test.rate.severe <- dat.test$test.rate.severe

FS1.A <- tm_shape(SC.pref) +
  tm_polygons("test.rate.all", palette = "Greens", title = "Percentage\ncases being\ntested")+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(A) All clinical cases", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FS1.B <- tm_shape(SC.pref) +
  tm_polygons("test.rate.mild", palette = "Greens", title = "Percentage\ncases being\ntested")+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(B) Mild cases", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FS1.C <- tm_shape(SC.pref) +
  tm_polygons("test.rate.severe", palette = "Greens", title = "Percentage\ncases being\ntested")+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(C) Severe cases", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

# tmap_save(tmap_arrange(FS1.A, FS1.B, FS1.C, ncol = 3), "Result\\Figure\\Fig S1.pdf", width = 12, height = 4)






######################## Figure S3 ###############################
# existing result as candidate solution for GA
# for all
load("Disease_model_realizations_6.RData")    # output of Disease system model realizations.R

# process data to get the data for alternative designs
dat.agg.alter <- dat.6 %>%
  mutate(Z_mild = Z_mild_1 + Z_mild_2 + Z_mild_3, Z_severe = Z_severe_1 + Z_severe_2 + Z_severe_3, inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Year, Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop)) %>%
  mutate(inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop), inc = mean(inc), inc_severe = mean(inc_severe)) %>%
  mutate(Z = Z_mild + Z_severe)

# Append the data to the shapefile for plotting maps
SC.pref$existing <- dat.agg.alter$Z/sum(dat.agg.alter$Z)
SC.pref$equal <- 1/21
SC.pref$popsize <- dat.agg.alter$pop/sum(dat.agg.alter$pop)
SC.pref$case <- dat.agg.alter$Y/sum(dat.agg.alter$Y)
SC.pref$IncRate <- dat.agg.alter$inc/sum(dat.agg.alter$inc)
SC.pref$CaseSevere <- dat.agg.alter$Y_severe/sum(dat.agg.alter$Y_severe)
SC.pref$IncRateSevere <- dat.agg.alter$inc_severe/sum(dat.agg.alter$inc_severe)

FigS3.A <- tm_shape(SC.pref) +
  tm_polygons("existing", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(A)")~bolditalic(Existing)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FigS3.B <- tm_shape(SC.pref) +
  tm_polygons("equal", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(B)")~bolditalic(Equal)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FigS3.C <- tm_shape(SC.pref) +
  tm_polygons("popsize", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(C)")~bolditalic(PopSize)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FigS3.D <- tm_shape(SC.pref) +
  tm_polygons("case", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(D)")~bolditalic(Case)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FigS3.E <- tm_shape(SC.pref) +
  tm_polygons("IncRate", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(E)")~bolditalic(IncRate)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FigS3.F <- tm_shape(SC.pref) +
  tm_polygons("CaseSevere", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(F)")~bolditalic(SevereCase)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

FigS3.G <- tm_shape(SC.pref) +
  tm_polygons("IncRateSevere", palette = "Blues", title = "Proportion\nof subtyping")+
  tm_text("name", size = 0.7) +
  tm_layout(title = expression(bold("(G)")~bolditalic(SevereIncRate)), legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

# tmap_save(tmap_arrange(FigS3.A, FigS3.B, FigS3.C, FigS3.D, FigS3.E, FigS3.F, FigS3.G, nrow = 3), "Result\\Figure\\Fig S3.pdf", width = 12, height = 12)







######################## Figure S4 ###############################
SC.pref$inc <- dat.agg.alter$inc
SC.copy <- SC.pref@data

SC.copy.long <- SC.copy %>%
  dplyr::select(name, inc, existing, equal, popsize, case, IncRate, CaseSevere, IncRateSevere) %>%
  gather("Type", "design", existing:IncRateSevere) %>%
  mutate(Type = factor(Type, levels = c("existing", "equal", "popsize", "case", "IncRate", "CaseSevere", "IncRateSevere")))

levels(SC.copy.long$Type) <- c("(A)*italic(' Existing')", "(B)*italic(' Equal')", "(C)*italic(' PopSize')", "(D)*italic(' Case')", "(E)*italic(' IncRate')", "(F)*italic(' SevereCase')", "(G)*italic(' SevereIncRate')")


SC.copy.long %>%
  arrange(-inc) %>%
  ggplot(aes(area = design, fill = inc, label = name)) +
  geom_treemap(layout = "fixed", start = "topleft", col = "white", size = 1.5) +
  geom_treemap_text(layout = "fixed", start = "topleft", grow = FALSE, reflow = TRUE, place = "center", size = 10, min.size = 1) +
  facet_wrap(~Type, ncol = 3, labeller = label_parsed) +
  scale_fill_material("deep-orange") +
  labs(fill = "Annual mean incidence rate (1/100,000)") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = "bold", margin = margin(0,0,5,0)), 
        legend.position = c(0.7, 0.15), 
        legend.key.width = unit(1.5, "cm"), 
        legend.direction = "horizontal",
        legend.title = element_text(size = 12, vjust = 1), 
        legend.text = element_text(size = 10), 
        legend.key.height = unit(0.3, "cm"), 
        panel.border = element_rect(size = 3, color = "red"),
        panel.spacing.y = unit(1, "lines"),
        legend.justification = "center", 
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0)
  ) +
  guides(fill = guide_colorbar(title.position = "top"))
# ggsave("Result\\Figure\\Fig S4.pdf", width = 8, height = 9)








######################## Figure S5 ###############################
load("Result\\Optimal_psevere_benchmark_scenario.Rdata")  # output of Optimal psevere archetypal designs.R

levels(result.all$goal) <- c("(A) For all", "(B) For severe")
levels(result.all$design) <- c("Existing", "Equal", "PopSize", "Case", "IncRate", "SevereCase", "SevereIncRate")

# smooth
result.all.smooth <- result.all %>%
  mutate(ofv = -ofv) %>%
  group_by(goal, design) %>%
  do(data.frame(., pred = predict(gam(ofv ~ s(p_severe, bs = "gp"), data = .)))) 

result.all.min.smooth <- result.all %>%
  mutate(ofv = -ofv) %>%
  group_by(goal, design) %>%
  do(data.frame(., pred = predict(gam(ofv ~ s(p_severe, bs = "gp"), data = .)))) %>%
  summarize(ofv.min = min(pred), p_severe = p_severe[which.min(pred)])

write.csv(result.all.min.smooth, "Result\\Optimal_psevere_for_benchmarking.csv", row.names = FALSE)

result.all.smooth %>%
  ggplot() +
  geom_point(aes(x = p_severe, y = ofv, col = design, group = design)) +
  geom_line(aes(x = p_severe, y = pred, col = design)) +
  facet_wrap(~goal, scales = "free") +
  labs(col = "") +
  xlab("Probability of severe cases being serotyped") +
  ylab("MAE (1/100,000)") +
  theme(legend.position = c(0.8, 0.85), legend.text = element_text(face = "italic", size = 10), strip.background = element_blank(), strip.text = element_text(size = 14, hjust = 0, face = "bold"), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  #  geom_point(data = result.all.min, aes(x = p_severe, y = ofv.min)) +
  geom_point(data = result.all.min.smooth, aes(x = p_severe, y = ofv.min, shape = goal), size = 2.5) +
  guides(shape = FALSE)

# ggsave("Result\\Figure\\Fig S5.pdf", width = 8, height = 4)







######################## Figure S6 ###############################
# original constraint, incidence rate, including suggestions
tm.optimal.all <- tm_shape(SC.pref) +
  tm_polygons("optimal.all", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(A) For all", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))



# original constraint, pathogen composition, including suggestions
tm.optimal.severe <- tm_shape(SC.pref) +
  tm_polygons("optimal.severe", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(B) For severe", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

tmap_arrange(tm.optimal.all, tm.optimal.severe)

# tmap_save(tmap_arrange(tm.optimal.all, tm.optimal.severe), "Result\\Figure\\Fig S6.pdf", width = 8, height = 3.75)






######################## Figure S7 ###############################
SC.pref$Y <- dat.agg.alter$Y
SC.pref$Z.optimal.inc <- sum(dat.agg.alter$Z) * SC.pref$optimal.all
SC.pref$Z.optimal.severe.inc <- sum(dat.agg.alter$Z) * SC.pref$optimal.severe
SC.pref$p.test.inc <- SC.pref$Z.optimal.inc/SC.pref$Y
SC.pref$p.test.inc.severe <- SC.pref$Z.optimal.severe.inc/SC.pref$Y

# original constraint, incidence rate, including suggestions
tm.optimal.all.prop <- tm_shape(SC.pref) +
  tm_polygons("p.test.inc", palette = "Blues", title = "Proportion of\ncases being\nsubtyped", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(A) For all", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))



# original constraint, pathogen composition, including suggestions
tm.optimal.severe.prop <- tm_shape(SC.pref) +
  tm_polygons("p.test.inc.severe", palette = "Blues", title = "Proportion of\ncases being\nsubtyped", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(B) For severe", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))


tmap_arrange(tm.optimal.all.prop, tm.optimal.severe.prop)

# tmap_save(tmap_arrange(tm.optimal.all.prop, tm.optimal.severe.prop), "Result\\Figure\\Fig S7.pdf", width = 8, height = 3.75)






######################## Figure S8 ###############################
obj.ofv <- data.frame(constraints = rep(c("*0.5", "*1", "*2", "*5"), 2), obj = rep(c("(A) For all", "(B) For severe"), each = 4), ofv = c(-GA.All.X0.5@fitnessValue, -GA.All.X1@fitnessValue, -GA.All.X2@fitnessValue, -GA.All.X5@fitnessValue, -GA.severe.X0.5@fitnessValue, -GA.severe.X1@fitnessValue, -GA.severe.X2@fitnessValue, -GA.severe.X5@fitnessValue))

obj.ofv %>%
  ggplot(aes(x = constraints, y = ofv, group = obj)) +
  geom_point() +
  geom_line() +
  facet_wrap(~obj, scales = "free") +
  xlab("Typing resources") +
  ylab("MAE (1/100,000)") +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 14, face = "bold", margin = margin(0,0,5,0))) +
  scale_x_discrete(labels = c("*0.5", "Observed", "*2", "*5"))
# ggsave("Result\\Figure\\Fig S8.pdf", width = 8, height = 4)






######################## Figure S9 ###############################
# For all, reousrce * 0.5
FigS9.A <- tm_shape(SC.pref) +
  tm_polygons("optimal.all.half", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(A) For all, resource*0.5", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))


# For all, reousrce * 2
FigS9.B <- tm_shape(SC.pref) +
  tm_polygons("optimal.all.2", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(B) For all, resource**2", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

# For all, reousrce * 5
FigS9.C <- tm_shape(SC.pref) +
  tm_polygons("optimal.all.5", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(C) For all, resource*5", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))



# severe inc rate, reousrce * 0.5
FigS9.D <- tm_shape(SC.pref) +
  tm_polygons("optimal.severe.half", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(D) For severe, resource*0.5", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))

# severe inc rate, reousrce * 2
FigS9.E <- tm_shape(SC.pref) +
  tm_polygons("optimal.severe.2", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(E) For severe, resource*2", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))


# severe inc rate, reousrce * 5
FigS9.F <- tm_shape(SC.pref) +
  tm_polygons("optimal.severe.5", palette = "Blues", title = "Proportion\nof subtyping", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12))+
  tm_text("name", size = 0.7) +
  tm_layout(title = "(F) For severe, resource*5", legend.position = c("left", "bottom"), legend.text.size = 0.7, inner.margins = c(0,0.1,0,0), outer.margins = c(0,0,0,0), legend.title.size = 0.9, frame = FALSE, title.size = 1.1, title.fontface = "bold", title.position = c(0, 0.98))


# tmap_save(tmap_arrange(FigS9.A, FigS9.D, FigS9.B, FigS9.E, FigS9.C, FigS9.F, nrow = 3), "Result\\Figure\\Fig S9.pdf", width = 8, height = 12)






######################## Figure S10 ###############################
SC.pref$inc.design <- SC.pref$inc/sum(SC.pref$inc)

FigS10.A <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.all.half, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = inc.design, yend = optimal.all.half , col = as.factor(inc.design - optimal.all.half  > 0)), size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "royalblue1", shape = 17) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = inc.design), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(A) For all, resource*0.5") +
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 10), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0.2, unit = "cm"), axis.title.x = element_text(margin = margin(-1,0,0,0)))

FigS10.B <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.all.2, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = inc.design, yend = optimal.all.2, col = as.factor(inc.design - optimal.all.2 > 0)),  size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "royalblue1", shape = 17) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = inc.design), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(B) For all, resource*2")+
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 10), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0.2, unit = "cm"), axis.title.x = element_text(margin = margin(-1,0,0,0)))

FigS10.C <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.all.5, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = inc.design, yend = optimal.all.5, col = as.factor(inc.design - optimal.all.5 > 0)),  size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "royalblue1", shape = 17) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = inc.design), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(C) For all, resource*5")+
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 10), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0.2, unit = "cm"), axis.title.x = element_text(margin = margin(-1,0,0,0)))

FigS10.D <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.severe.half, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = inc.design, yend =optimal.severe.half, col = as.factor(inc.design - optimal.severe.half  > 0)),  size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "blue3", shape = 15) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = inc.design), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(D) For severe, resource*0.5") +
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 10), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0.2, unit = "cm"), axis.title.x = element_text(margin = margin(-1,0,0,0)))

FigS10.E <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.severe.2, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = inc.design, yend = optimal.severe.2, col = as.factor(inc.design - optimal.severe.2 > 0)),  size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "blue3", shape = 15) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = inc.design), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(E) For severe, resource*2")+
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 10), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0.2, unit = "cm"), axis.title.x = element_text(margin = margin(-1,0,0,0)))

FigS10.F <- SC.pref@data %>%
  ggplot(aes(x = inc, y = optimal.severe.5, label = name)) +
  geom_smooth(method = "lm") + 
  geom_segment(aes(x = inc, xend = inc, y = inc.design, yend = optimal.severe.5, col = as.factor(inc.design - optimal.severe.5 > 0)),  size = 0.5, arrow = arrow(angle = 10, type = "closed", length = unit(0.03, "npc"))) +
  scale_color_manual(values = c("#d73027", "#1a9850")) +
  geom_point(size = 2, col = "blue3", shape = 15) +
  geom_text_repel(size = 2.8, angle = 90) +
  geom_abline(slope = 1/sum(SC.pref$inc), intercept = 0) + 
  geom_point(aes(x = inc, y = inc.design), size = 2) +
  guides(col = FALSE) +
  xlab("Annual mean incidence rate (1/100,000)") +
  ylab("Proportion of typing") +
  ggtitle("(F) For severe, resource*5")+
  theme(plot.title = element_text(hjust = -0.3, size = 12), axis.text.y = element_text(angle = 90, size = 10, hjust = 0.5), axis.title = element_text(size = 10), axis.text.x = element_text(size = 10), plot.margin = margin(t = 0.2, unit = "cm"), axis.title.x = element_text(margin = margin(-1,0,0,0)))


# save_plot("Result\\Figure\\Fig S10.pdf", plot_grid(FigS10.A, FigS10.B, FigS10.C, FigS10.D, FigS10.E, FigS10.F, ncol = 3), base_width = 10, base_height =6.7)



######################## Figure S11 ###############################
load("Result\\OBJ_testlimit_cor.Rdata")     # output of "Obj value of random designs under different resource scenarios.R"

ofv.result[ofv.result == -10] <- NA
M1 <- cor(ofv.result[, 1:4], use = "pairwise.complete.obs")
colnames(M1) <- rownames(M1) <- c("*0.5", "Observed", "*2", "*5")


M2 <- cor(ofv.result[, 5:8], use = "pairwise.complete.obs")
colnames(M2) <- rownames(M2) <- c("*0.5", "Observed", "*2", "*5")

pdf("Result\\Figure\\Fig S11.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
corrplot(M1, method = "color", type = "upper", addCoef.col = "white", tl.col = "black", tl.cex = 0.8, number.cex = 0.7, addgrid.col = "white", tl.srt = 30, cl.ratio = .2, cl.align = "l")
corrplot(M2, method = "color", type = "upper", addCoef.col = "white", tl.col = "black", tl.cex = 0.8, number.cex = 0.7, addgrid.col = "white", tl.srt = 30, cl.ratio = .2, cl.align = "l")
dev.off()
