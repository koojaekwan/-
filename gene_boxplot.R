
# library -----------------------------------------------------------------
library(tidyverse)
library(data.table)

library(glmnet)
library(caret)
library(Metrics)

library(doParallel)


doParallel::registerDoParallel(3)


load("D:/Jae Kwan/4학년여름/연구생/연구주제/glm.RData")


# one_hot_encoding --------------------------------------------------------

one_hot_meta_dat <- model.matrix(~.-1, meta_dat)
colnames(one_hot_meta_dat) <- one_hot_meta_dat %>% colnames %>% 
  str_remove(.,"condition")

# dummyVars(" ~ .", data = meta_dat)->kkkk
# predict(kkkk,meta_dat) %>% data.frame 



# cancer ------------------------------------------------------------------

normal_t <- normal %>% t %>% as.matrix
normal_t[1:3,1:2]
dim(normal_t)


cancer_fit <- glmnet(x = normal_t, 
                     y = one_hot_meta_dat,
                     family = "multinomial", 
                     type.multinomial = "grouped", 
                     alpha = 1) # alpha=1 : lasso

par(mfrow=c(2,2))
plot(cancer_fit)

cancer_fit$npasses
cancer_fit$lambda

# cv fit ------------------------------------------------------------------
set.seed(100)

cv_fit <- 
cv.glmnet(x = normal_t,
          y = one_hot_meta_dat,
          family = "multinomial",
          type.multinomial = "grouped",
          alpha = 1,
          parallel = T)

par(mfrow=c(1,1))
plot(cv_fit)

par(mfrow=c(2,2))
plot(cv_fit$glmnet.fit, xvar="lambda")
plot(cv_fit$glmnet.fit, xvar="norm")
par(mfrow=c(1,1))


cv_fit$lambda.1se
cv_fit$lambda.min

# predict(cv_fit, newx = normal_t, s = "lambda.1se", type = "class") -> a
# 
# table(a, meta_dat$condition)
# mean(a==meta_dat)


# cancer_fit  coef --------------------------------------------------------


# temp <- predict(cancer_fit, newx = normal_t, 
#                 s = cv_fit$lambda.min, 
#                 type = "class")   # s : lambda
# 
# temp %>% head
# temp %>% tail
# 
# mean(temp==meta_dat)
# table(prediction = temp, actual = meta_dat$condition)



tmp_coeffs <- coef(cancer_fit, s = cv_fit$lambda.1se)

data.frame(name = tmp_coeffs[[1]]@Dimnames[[1]][tmp_coeffs[[1]]@i + 1], coefficient = tmp_coeffs[[1]]@x) -> k1

data.frame(name = tmp_coeffs[[2]]@Dimnames[[1]][tmp_coeffs[[2]]@i + 1], coefficient = tmp_coeffs[[2]]@x) -> k2

data.frame(name = tmp_coeffs[[3]]@Dimnames[[1]][tmp_coeffs[[3]]@i + 1], coefficient = tmp_coeffs[[3]]@x) -> k3

data.frame(name = tmp_coeffs[[4]]@Dimnames[[1]][tmp_coeffs[[4]]@i + 1], coefficient = tmp_coeffs[[4]]@x) -> k4


all(k1$name==k2$name)
all(k2$name==k3$name)
all(k3$name==k4$name)

k1 %>% 
  left_join(k2, by = "name") %>% 
  left_join(k3, by = "name") %>% 
  left_join(k3, by = "name")



# normalized count & gene condition ---------------------------------------
normal_temp <- normal %>% t %>% data.frame
normal_temp[1:2,1:2]

normal_z <- scale(normal_temp, center = T, scale = T)
normal_z <- normal_z %>% t %>% data.frame
normal_z[1:2,1:2]




sig_gene <- k1 %>% 
  select(name) %>% 
  left_join(normal_z %>% rownames_to_column(var="name"), by = c("name"))

sig_gene <- sig_gene[-1,]
sig_gene[1:2,1:3]


rownames(sig_gene) <- sig_gene[,1]
sig_gene <- sig_gene[,-1]

sig_gene <- sig_gene %>% t %>% data.frame
sig_gene <- sig_gene %>% rownames_to_column(var = "name")

sig_gene$name <- gsub(pattern = "\\.",replacement = "-",x = sig_gene$name)

sig_gene <- sig_gene %>% left_join(meta_dat %>% rownames_to_column("name"), by = c("name"))

colnames(sig_gene)

work <- sig_gene %>% select(1:11,"condition")
work %>% colnames
sig_gene %>% dim

work[1:2,1:4]

work_temp <- 
work %>%
  pivot_longer(cols = -c(condition, name),
               names_to = "gene", values_to = "count")

work_temp %>% head

work_temp$count <- as.numeric(work_temp$count)


library(hrbrthemes)
library(viridis)

work_temp %>%
  ggplot(aes(x=gene,y=count,fill=condition))+
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.06, alpha=0.2)+
  theme_ipsum() +
  theme(
    legend.position=c(.9, .8),
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot by subtype") +
  xlab("") +
  ylab("Normalized Count") +
  labs(fill="Subtype") + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

