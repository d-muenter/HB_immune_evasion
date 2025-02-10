library(tidyverse)
library(ggpubr)
library(RColorBrewer)


dir.create("Output/cell_quantification/", r=T)


### lymphoid cells
tumor_lymphoid <- c(0.008419779, 0.001076368, 0.0072331, 0.028723964, 0.001832268, 
                    0.025223084, 0.002102141) # proportion of lymphoid cells in tumor cell areas
non_tumor_lymphoid <- c(0.135850779, 0.233406556, 0.020482809, 0.082828637, 
                        0.048758155, 0.161944751, 0.021841705) # proportion of lymphoid cells in non-tumor cell areas

test_ly <- wilcox.test(tumor_lymphoid, non_tumor_lymphoid, paired=T)

lymphoid_df <- data.frame(proportion = c(tumor_lymphoid, non_tumor_lymphoid), group = rep(c("Tumor", "Non-tumor"), each=7))

lymphoid_df_mean <- lymphoid_df %>% 
  group_by(group) %>% 
  summarise(
    n = n(), 
    mean = mean(proportion), 
    sd=sd(proportion)
  ) %>% 
  mutate(se = sd/sqrt(n))

colors <- brewer.pal(7, "Spectral")[c(7,1)]

ly <- ggplot(lymphoid_df, aes(x = group, y = proportion)) + 
  geom_col(data=lymphoid_df_mean, aes(x=group, y=mean, fill=group), width=0.7) + 
  geom_errorbar(data=lymphoid_df_mean, aes(x=group, y=mean, ymin=mean-se, ymax=mean+se), width=0.3) + 
  geom_point(data = lymphoid_df, aes(x = group, y = proportion), color="black", size=2.5) + 
  geom_bracket(xmin=1, xmax=2, y.position=0.26, 
               label = paste0("p = ", round(test_ly$p.value, 3)), 
               tip.length=c(0.02, 0.2), vjust=-.5) + 
  theme_minimal() + 
  scale_fill_manual(lymphoid_df_mean$group, values = colors) + 
  theme(axis.title.x = element_blank(), panel.border = element_rect(fill=NA), 
        legend.position = "none", text = element_text(size=20), 
        axis.title.y = element_text(vjust=2), panel.grid=element_blank()) + 
  labs(y = "Lymphoid Cell Proportion") + 
  ylim(0, 0.27)

tiff("Output/cell_quantification/lymphoid_cell_proportions.tiff", width=300, height=4800)
print(ly)
dev.off()



### myeloid cells
tumor_myeloid <- c(0.068832953, 0.018782627, 0.01910007, 0.011655012, 
                   0.014582114, 0.006300025, 0.01828374) # proportion of myeloid cells in tumor cell areas
non_tumor_myeloid <- c(0.094665101, 0.256607945, 0.011704462, 0.042287174, 
                                0.001049178, 0.067064954, 0.050431613) # proportion of myeloid cells in non-tumor cell areas

test_my <- wilcox.test(tumor_myeloid, non_tumor_myeloid, paired=T)

myeloid_df <- data.frame(proportion = c(tumor_myeloid, non_tumor_myeloid), group = rep(c("Tumor", "Non-tumor"), each=7))
myeloid_df_mean <- myeloid_df %>% 
  group_by(group) %>% 
  summarise(
    n = n(), 
    mean = mean(proportion), 
    sd=sd(proportion)
  ) %>% 
  mutate(se = sd/sqrt(n))

colors <- brewer.pal(7, "Spectral")[c(7,1)]

my1 <- ggplot(myeloid_df, aes(x = group, y = proportion)) + 
  geom_col(data=myeloid_df_mean, aes(x=group, y=mean, fill=group), width=0.7) + 
  geom_errorbar(data=myeloid_df_mean, aes(x=group, y=mean, ymin=mean-se, ymax=mean+se), width=0.3) + 
  geom_point(data = myeloid_df, aes(x = group, y = proportion), color="black", size=2.5) + 
  geom_bracket(xmin=1, xmax=2, y.position=0.28, 
               label = "ns", 
               tip.length=c(0.02, 0.2), vjust=-.5) + 
  theme_minimal() + 
  scale_fill_manual(myeloid_df_mean$group, values = colors) + 
  theme(axis.title.x = element_blank(), panel.border = element_rect(fill=NA), 
        legend.position = "none", text = element_text(size=20), 
        axis.title.y = element_text(vjust=2), panel.grid=element_blank()) + 
  labs(y = "Myeloid Cell Proportion") +
  ylim(0, 0.29)


tiff("Output/cell_quantification/myeloid_cell_proportions.tiff", width=300, height=480)
print(my1)
dev.off()





### CD204+ cells
# number of CD204+ cells in tumor/non-tumor
tumor_CD204_n <- c(924, 352, 189, 10, 1414, 94, 35)
non_tumor_CD204_n <- c(35, 153, 0, 144, 1, 260, 189)

# number of myeloid cells in tumor/non-tumor
tumor_myeloid_n <- c(5461, 698, 1793, 155, 1918, 281, 1122)
non_tumor_myeloid_n <- c(2903, 17873, 16, 4336, 55, 1493, 2629)

# CD204+ proportion of myeloid cells
tumor_CD204_by_myeloid <- tumor_CD204_n / tumor_myeloid_n
non_tumor_CD204_by_myeloid <- non_tumor_CD204_n / non_tumor_myeloid_n

CD204_by_mye_df <- data.frame(proportion = c(tumor_CD204_by_myeloid, non_tumor_CD204_by_myeloid), group = rep(c("Tumor", "Non-tumor"), each=7))

test_CD204 <- wilcox.test(tumor_CD204_by_myeloid, non_tumor_CD204_by_myeloid, paired=T)

CD204_by_mye_df <- data.frame(proportion = c(tumor_CD204_by_myeloid, non_tumor_CD204_by_myeloid), group = rep(c("Tumor", "Non-tumor"), each=7))
CD204_by_mye_df_mean <- CD204_by_mye_df %>% 
  group_by(group) %>% 
  summarise(
    n = n(), 
    mean = mean(proportion), 
    sd=sd(proportion)
  ) %>% 
  mutate(se = sd/sqrt(n))

colors <- brewer.pal(7, "Spectral")[c(7,1)]

CD204 <- ggplot(CD204_by_mye_df, aes(x = group, y = proportion)) + 
  geom_col(data=CD204_by_mye_df_mean, aes(x=group, y=mean, fill=group), width=0.7) + 
  geom_errorbar(data=CD204_by_mye_df_mean, aes(x=group, y=mean, ymin=mean-se, ymax=mean+se), width=0.3) + 
  geom_point(data = CD204_by_mye_df, aes(x = group, y = proportion), color="black", size=2.5) + 
  geom_bracket(xmin=1, xmax=2, y.position=0.8, 
               label = paste0("p = ", round(test_CD204$p.value, 3)), 
               tip.length=c(0.2, 0.02), vjust=-.5) + 
  theme_minimal() + 
  scale_fill_manual(CD204_by_mye_df_mean$group, values = colors) + 
  theme(axis.title.x = element_blank(), panel.border = element_rect(fill=NA), 
        legend.position = "none", text = element_text(size=20), 
        axis.title.y = element_text(vjust=2), panel.grid=element_blank()) + 
  labs(y = "CD204+ Cell Proportion") + 
  ylim(0, 0.82)

tiff("Output/cell_quantification/CD204_cell_proportions.tiff", width=300, height=480)
print(CD204)
dev.off()

