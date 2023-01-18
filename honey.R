# load libraries

library(tidyverse)
library(brms)
library(patchwork)
library(cmdstanr)
library(here)
library(GGally)
library(tidybayes)
library(kableExtra)
library(bayesplot)
library(emmeans)
library(broom.mixed)
library(ggdist)
library(posterior)
library(easystats)

options(mc.cores = 4, brms.backend = "cmdstanr")

options(warn = -1)

theme_set(theme_bw())

# import data

all_honey <- read.csv(here("data","all_honeys.csv"), header = TRUE, sep=",")
all_honey$origin <- as.factor(all_honey$origin)
# to produce centered moisture values:
all_honey <- all_honey %>% mutate(mc=m-mean(m))
# dataset to validate from Beckh:
beckh <- read.csv(here("data","beckh.csv"), header = TRUE, sep=",")
beckh$origin <- as.factor(beckh$origin)
beckh <- beckh %>% mutate(mc=m-mean(m))

model_honey <- read.csv(here("data","model_honey.csv"), header = TRUE, sep=",")
model_honey$source <- as.factor(model_honey$source)
model_honey <- model_honey %>% mutate(mc=m-mean(all_honey$m))

all_honey_means <- all_honey %>% group_by(origin) %>% summarize(mean_aw=mean(aw))
all_honey_means$origin <- as.factor(all_honey_means$origin)

glucose <- read.csv(here("data","glucose.csv"), header=TRUE, sep=",")

# code for Figure 2

glucose %>%  ggplot(aes(x=m, y=aw))+
  geom_line()+
  geom_segment(x = 12, xend = 28, y=0.43, yend=0.75,lty=2, size=1.1, color="red")+
  labs(x="moisture content m (%)", y=expression("a"[w]))+
  geom_segment(x=12, xend=12, y=0, yend=0.43, lty=2, size=1.1,color="red")+
  geom_segment(x=28, xend=28, y=0, yend=0.75, lty=2, size=1.1, color="red")

# code for Figure 3:

(all_honey_plot <- ggplot(all_honey, aes(x=m, y=aw), )+
  geom_point(aes(color=origin), shape=21)+
   labs(x="moisture content m (%)", y=expression("a"[w]))
)

# Complete pooling with the null model with brms:

honeys_pooled <- 
  brm(data=all_honey, family=gaussian, aw ~ 1, 
      prior = c(set_prior("normal(0.55, 0.1)", class = "Intercept"),
                set_prior("exponential(1)", class = "sigma")),
                iter=4000, warmup=2000, chains=4, refresh=0,
                file=here("fits","honeys_pooled"))
  
print(honeys_pooled, digits=4)
honeys_pooled_post <- as_draws_df(honeys_pooled)


# brms code for no-pooling with the null model

honeys_unpooled <- brm(data=all_honey, family=gaussian,
                       aw ~ 0 + factor(origin),
                         prior(normal(0.55,0.1), class=b),
                         prior(exponential(1), class=sigma),
                         iter=4000, warmup=2000, chains=4, refresh=0,
                         file=here("fits","honeys_unpooled")
                         )
print(honeys_unpooled)
honeys_unpooled_post <- as_draws_df(honeys_unpooled)

# brms code for partial-pooling with the null model

honeys_partial <- 
  brm(data = all_honey, family = gaussian,
      aw ~ 1+(1|origin),
      prior = c(set_prior("normal(0.55, 0.1)", class = "Intercept"),
                set_prior("exponential(1)", class = "sd"),
                set_prior("exponential(1)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, cores = 4, refresh=0, file=here("fits","honeys_partial"))

print(honeys_partial)
ICC_null <- performance::icc(honeys_partial)
bayestestR:: plot(honeys_partial)
honeys_partial_post <- as_draws_df(honeys_partial)

# code for Figure 4 with bayesplot:

predictions_multi_null <- posterior_predict(honeys_partial, newdata=all_honey_means)

ppc_intervals(all_honey_means$mean_aw, yrep=predictions_multi_null, prob_outer = 0.95)+
  ggplot2::scale_x_continuous(labels = all_honey_means$origin, breaks = 1:nrow(all_honey_means))+
  labs(x="origin", y="mean water activity")+
    geom_hline(yintercept = fixef(honeys_partial)[1], lty=2)+
  theme(legend.position = "none")


# code for Figure 5

mean_plot <- ggplot(data=honeys_pooled_post, aes(x=b_Intercept))+
  geom_density(fill="blue", alpha=0.5)+
  geom_density(data=honeys_partial_post, aes(x=b_Intercept), fill="red")+
  annotate("label", x=0.58, y=275, label="partially pooled", fill="red", alpha=0.5, size=3)+
  annotate("label", x=0.58, y=300, label="completely pooled", fill="blue", alpha=0.5, size=3)+
  labs(x=expression("mean a"[w]), y="density", subtitle = "A")

sd_plot <- ggplot(data=honeys_pooled_post, aes(x=sigma))+
  geom_density(fill="blue", alpha=0.5)+
  geom_density(data=honeys_partial_post, aes(x=sigma), fill="red")+
  geom_density(data=honeys_partial_post,aes(x=sd_origin__Intercept), fill="cyan", alpha=0.5)+
  annotate("label", x=0.05, y=800, label=expression(paste("partially pooled  ",sigma[r])), fill="red", alpha=0.5, size=3)+
  annotate("label", x=0.05, y=600, label=expression(paste("completely pooled  ",sigma[r])), fill="blue", alpha=0.5, size=3)+
    annotate("label", x=0.05, y=700, label=expression(paste("partially pooled  ",sigma[b])), fill="cyan", alpha=0.5, size=3)+
  labs(x=expression("standard deviation a"[w]), y="",subtitle = "B")
mean_plot+sd_plot


# code for Table 2 with kable

honeys_var1 <- summary(honeys_partial)
honeys_var2 <- rbind(data.frame(honeys_var1$fixed), data.frame(honeys_var1$random$origin),
                     data.frame(honeys_var1$spec_pars))
rownames(honeys_var2) <- c("$\\mu$", "$\\sigma_b$", "$\\sigma_r$")
colnames(honeys_var2) <- c("mean","SE", "lower bound", "upper bound")
honeys_var2[1:3,1:4] %>% 
    rownames_to_column(var = "parameter") %>% kbl(booktabs=T, escape=F, digits = c(2,3,3,3,3)) %>% kable_styling(position="center", full_width = F)

# brms code for pooled regression with predictor

honeys_pooled2 <- brm(data=all_honey, family = gaussian,
                      aw ~ 1 + mc,
                      prior(normal(0.58,0.1), class=Intercept),
                      prior(normal(0.05,0.01), class=b),
                      prior(exponential(1), class=sigma),
                      iter=4000, warmup=2000, chains=4, refresh=0,
                      file=here("fits","honeys_pooled2"))
print(honeys_pooled2, digits=3)

honeys_pooled2_post <- as_draws_df(honeys_pooled2)

# subset of data (origin 1-8) to show Simpson's paradox
all_honey2 <- all_honey[1:823,]

honeys_pooled2a <- brm(data=all_honey2, family = gaussian,
                      aw ~ 1 + mc,
                      prior(normal(0.58,0.1), class=Intercept),
                      prior(normal(0.05,0.01), class=b),
                      prior(exponential(1), class=sigma),
                      iter=4000, warmup=2000, chains=4, refresh=0,
                      file=here("fits","honeys_pooled2a"))
print(honeys_pooled2a, digits=3)
honeys_pooled2a_post <- as_draws_df(honeys_pooled2a)

# code for Figure 6

all_honey_plot1 <- all_honey %>% ggplot(aes(x=mc, y=aw, color=origin))+
  geom_point(shape=21)+
 geom_abline(data=honeys_pooled2_post,aes(intercept=fixef(honeys_pooled2)[1,1], slope = fixef(honeys_pooled2)[2,1]))+
   labs(x=expression("m"[c]), y=expression("a"[w]), subtitle = "A")+
  theme(legend.position = "none")


all_honey_plot2 <- all_honey2 %>% ggplot(aes(x=mc, y=aw, color=origin))+
  geom_point(shape=21)+
 geom_abline(data=honeys_pooled2a_post,aes(intercept=fixef(honeys_pooled2a)[1,1], slope = fixef(honeys_pooled2a)[2,1]))+
   labs(x=expression("m"[c]), y=expression("a"[w]), subtitle = "B")

all_honey_plot1 + all_honey_plot2

#  code for Table 3 using kable

c_honey1 <- summary(honeys_pooled2)
summary_honey1 <- rbind(data.frame(c_honey1$fixed), data.frame(c_honey1$spec_pars)) %>% select(-Rhat,-Bulk_ESS,-Tail_ESS)

rownames(summary_honey1) <- c("intercept $\\beta_0$", "slope $\\beta_1$", "$\\sigma_r$")
colnames(summary_honey1) <- c("mean","SE", "lower bound", "upper bound")

kbl(summary_honey1, digits = c(3,3,2,2,3,0,0), escape=F, booktabs=TRUE) %>% kable_styling(position = "center", full_width = F)

#  brms code for regression of the unpooled data using the map function from the purrr package

nested_honeys <- all_honey %>% nest(-origin)

model_nested_honeys <- nested_honeys %>% mutate(model= map(data, ~brm(data=., family=gaussian,
             aw ~ 1 + mc,
             prior(normal(0.55,0.1), class=Intercept),
             prior(normal(0.03,0.01), class=b),
             prior(exponential(1), class=sigma),
             iter=4000, warmup=2000, chains=4, refresh=0)))
# saving the result:
saveRDS(model_nested_honeys, here("fits","model_nested_honeys.rds"))

# loading the unpooled regression results:

model_nested_honeys <- readRDS(here("fits","model_nested_honeys.rds"))

#make prediction for every cultivar:
mc.seq <- data.frame(mc = seq(from = -5, to = 11, by = 0.1))
model_nested_honeys <- model_nested_honeys %>% mutate(model_pred = map(model, ~predict(., newdata=mc.seq))) 
# predict new data
model_nested_honeys <- model_nested_honeys %>% mutate(newdata=map(data,~bind_cols(mc.seq)))   
# bind new data to predictions and add posterior draws
model_nested_honeys <- model_nested_honeys %>% mutate(post=map(model,~as_draws_df(.)))

# code for Figure 7

pred1 <- model_nested_honeys$model_pred[[1]] %>% as_tibble() %>% bind_cols(mc.seq)
pred1$origin<- as.factor(rep("1", length(nrow(mc.seq))))

pred2 <- model_nested_honeys$model_pred[[2]] %>% as_tibble() %>% bind_cols(mc.seq)
pred2$origin<- as.factor(rep("2", length(nrow(mc.seq))))

pred3 <- model_nested_honeys$model_pred[[3]] %>% as_tibble() %>% bind_cols(mc.seq)
pred3$origin<- as.factor(rep("3", length(nrow(mc.seq))))

pred4 <- model_nested_honeys$model_pred[[4]] %>% as_tibble() %>% bind_cols(mc.seq)
pred4$origin<- as.factor(rep("4", length(nrow(mc.seq))))

pred5 <- model_nested_honeys$model_pred[[5]] %>% as_tibble() %>% bind_cols(mc.seq)
pred5$origin<- as.factor(rep("5", length(nrow(mc.seq))))

pred6 <- model_nested_honeys$model_pred[[6]] %>% as_tibble() %>% bind_cols(mc.seq)
pred6$origin<- as.factor(rep("6", length(nrow(mc.seq))))

pred7 <- model_nested_honeys$model_pred[[7]] %>% as_tibble() %>% bind_cols(mc.seq)
pred7$origin<- as.factor(rep("7", length(nrow(mc.seq))))

pred8 <- model_nested_honeys$model_pred[[8]] %>% as_tibble() %>% bind_cols(mc.seq)
pred8$origin<- as.factor(rep("8", length(nrow(mc.seq))))

pred9 <- model_nested_honeys$model_pred[[9]] %>% as_tibble() %>% bind_cols(mc.seq)
pred9$origin<- as.factor(rep("9", length(nrow(mc.seq))))

pred10 <- model_nested_honeys$model_pred[[10]] %>% as_tibble() %>% bind_cols(mc.seq)
pred10$origin<- as.factor(rep("10", length(nrow(mc.seq))))

pred11 <- model_nested_honeys$model_pred[[11]] %>% as_tibble() %>% bind_cols(mc.seq)
pred11$origin<- as.factor(rep("11", length(nrow(mc.seq))))

pred12 <- model_nested_honeys$model_pred[[12]] %>% as_tibble() %>% bind_cols(mc.seq)
pred12$origin<- as.factor(rep("12", length(nrow(mc.seq))))

pred13 <- model_nested_honeys$model_pred[[13]] %>% as_tibble() %>% bind_cols(mc.seq)
pred13$origin<- as.factor(rep("13", length(nrow(mc.seq))))

pred14 <- model_nested_honeys$model_pred[[14]] %>% as_tibble() %>% bind_cols(mc.seq)
pred14$origin<- as.factor(rep("14", length(nrow(mc.seq))))

pred15 <- model_nested_honeys$model_pred[[15]] %>% as_tibble() %>% bind_cols(mc.seq)
pred15$origin<- as.factor(rep("15", length(nrow(mc.seq))))

pred16 <- model_nested_honeys$model_pred[[16]] %>% as_tibble() %>% bind_cols(mc.seq)
pred16$origin<- as.factor(rep("16", length(nrow(mc.seq))))

pred17 <- model_nested_honeys$model_pred[[17]] %>% as_tibble() %>% bind_cols(mc.seq)
pred17$origin<- as.factor(rep("17", length(nrow(mc.seq))))

pred18 <- model_nested_honeys$model_pred[[18]] %>% as_tibble() %>% bind_cols(mc.seq)
pred18$origin<- as.factor(rep("18", length(nrow(mc.seq))))

pred19 <- model_nested_honeys$model_pred[[19]] %>% as_tibble() %>% bind_cols(mc.seq)
pred19$origin<- as.factor(rep("19", length(nrow(mc.seq))))

pred20 <- model_nested_honeys$model_pred[[20]] %>% as_tibble() %>% bind_cols(mc.seq)
pred20$origin<- as.factor(rep("20", length(nrow(mc.seq))))

pred21 <- model_nested_honeys$model_pred[[21]] %>% as_tibble() %>% bind_cols(mc.seq)
pred21$origin<- as.factor(rep("21", length(nrow(mc.seq))))

pred22 <- model_nested_honeys$model_pred[[22]] %>% as_tibble() %>% bind_cols(mc.seq)
pred22$origin<- as.factor(rep("22", length(nrow(mc.seq))))

pred23 <- model_nested_honeys$model_pred[[23]] %>% as_tibble() %>% bind_cols(mc.seq)
pred23$origin<- as.factor(rep("23", length(nrow(mc.seq))))

pred24 <- model_nested_honeys$model_pred[[24]] %>% as_tibble() %>% bind_cols(mc.seq)
pred24$origin<- as.factor(rep("24", length(nrow(mc.seq))))

pred25 <- model_nested_honeys$model_pred[[25]] %>% as_tibble() %>% bind_cols(mc.seq)
pred25$origin<- as.factor(rep("25", length(nrow(mc.seq))))

pred26 <- model_nested_honeys$model_pred[[26]] %>% as_tibble() %>% bind_cols(mc.seq)
pred26$origin<- as.factor(rep("26", length(nrow(mc.seq))))

pred27 <- model_nested_honeys$model_pred[[27]] %>% as_tibble() %>% bind_cols(mc.seq)
pred27$origin<- as.factor(rep("27", length(nrow(mc.seq))))

pred28 <- model_nested_honeys$model_pred[[28]] %>% as_tibble() %>% bind_cols(mc.seq)
pred28$origin<- as.factor(rep("28", length(nrow(mc.seq))))

pred29 <- model_nested_honeys$model_pred[[29]] %>% as_tibble() %>% bind_cols(mc.seq)
pred29$origin<- as.factor(rep("29", length(nrow(mc.seq))))

pred_unpooled <- rbind(pred1,pred2,pred3,pred4,pred5,pred6, pred7, pred8, pred9, pred10, pred11, pred12, pred13, pred14, pred15, pred16,pred17, pred18,pred19,pred20, pred21,pred22,pred23, pred24, pred25, pred26, pred27,pred28,pred29)

pred_unpooled_plot <- all_honey %>% ggplot(aes(x=mc, y=aw))+
  geom_point()+
  geom_ribbon(data=pred_unpooled, aes(x=mc, y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),fill = "lightblue", alpha=0.5)+
    geom_line(data=pred_unpooled, aes(x=mc, y=Estimate), color="blue")+
  facet_wrap(~origin, ncol=6)+
  labs(x=expression("centered m"[c]), y=expression("a"[w]))

pred_unpooled_plot  

# brms code for partially pooled regresion with predictor, varying intercept

priors_pooled <- c(prior(normal(0.55, 0.1), class = Intercept),
                   prior(normal(0.05,0.01), class = b),
                   prior(exponential(1), class = sd),
                   prior(exponential(1), class = sigma))
honeys_partial2 <- 
  brm(data = all_honey, family = gaussian,
      aw ~ 1+(1|origin) + mc,
      prior = priors_pooled,
      chains = 4, iter = 4000, warmup = 2000, cores = 4, refresh=0, file=here("fits","honeys_partial2"))

honeys_partial2_post <- as_draws_df(honeys_partial2)
print(honeys_partial2)
honeys_partial2$fit
performance::icc(honeys_partial2)


# Code for Figure 8

mean_plot2a <- ggplot(data=honeys_pooled2_post, aes(x=b_Intercept))+
  geom_density(fill="blue", alpha=0.5)+
  geom_density(data=honeys_partial2_post, aes(x=b_Intercept), fill="red")+
  annotate("label", x=0.57, y=275, label="partially pooled", fill="red", alpha=0.5, size=1.8)+
  annotate("label", x=0.57, y=300, label="completely pooled", fill="blue", alpha=0.5, size=1.8)+
  labs(x=expression(beta[0]), y="density", subtitle = "A")+
           theme(axis.text.x = element_text(angle=70, hjust=1))

mean_plot2b <- ggplot(data=honeys_pooled2_post, aes(x=b_mc))+
  geom_density(fill="blue", alpha=0.5)+
  geom_density(data=honeys_partial2_post, aes(x=b_mc), fill="red")+
  annotate("label", x=0.011, y=900, label="partially pooled", fill="red", alpha=0.5, size=1.8)+
  annotate("label", x=0.011, y=1000, label="completely pooled", fill="blue", alpha=0.5, size=1.8)+
  labs(x=expression(beta[1]), y="", subtitle = "B")+
           theme(axis.text.x = element_text(angle=70, hjust=1))

sd_plot2a <- ggplot(data=honeys_pooled2_post, aes(x=sigma))+
  geom_density(fill="blue", alpha=0.5)+
  geom_density(data=honeys_partial2_post, aes(x=sigma), fill="red")+
  geom_density(data=honeys_partial2_post,aes(x=sd_origin__Intercept), fill="cyan", alpha=0.5)+
  annotate("label", x=0.0355, y=1100, label=expression(paste("partially pooled  ",sigma[r])), fill="red", alpha=0.5, size=1.8)+
  annotate("label", x=0.0355, y=900, label=expression(paste("completely pooled  ",sigma[r])), fill="blue", alpha=0.5, size=1.8)+
    annotate("label", x=0.0355, y=1000, label=expression(paste("partially pooled  ",sigma[b])), fill="cyan", alpha=0.5, size=1.8)+
  labs(x=expression(sigma),y="", subtitle = "C")+
           theme(axis.text.x = element_text(angle=70, hjust=1))

mean_plot2a+mean_plot2b+sd_plot2a


# code for Figure 9

all_honey %>% ggplot(aes(x=mc, y=aw,color=origin))+
  geom_point(shape=21)+
  geom_abline(aes(intercept=fixef(honeys_partial2)[1,1], slope = fixef(honeys_partial2)[2,1]))+
   geom_abline(data=honeys_pooled2_post,aes(intercept=fixef(honeys_pooled2)[1,1], slope = fixef(honeys_pooled2)[2,1]), lty=2)+
   theme(legend.position = "none")

#| code for Figure 10
# coefficient plot via spread_rvars from the tidybayes package in combination with posterior package
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-posterior.html

varying_ic <- honeys_partial2 %>% 
  spread_rvars(r_origin[origin,coef], b_Intercept)%>%
  mutate(origin_intercept = b_Intercept + r_origin) %>% 
  ggplot(aes(y = factor(origin), dist = origin_intercept)) +
  stat_dist_halfeye() +
  labs(x=expression("intercept value a"[w]), y="origin")
varying_ic

# code for Figure 11

mc.seq <- expand_grid(origin=1:29, mc=seq(from = -5, to = 11, by = 0.1))
mc.seq$origin <- as.factor(mc.seq$origin)
pred_model_honeys1 <- cbind(mc.seq, predict(honeys_partial2, newdata=mc.seq, re_formula = NULL)[,-2])
pred_model_honeys2 <- cbind(mc.seq, predict(honeys_partial2, newdata=mc.seq, re_formula = NA)[,-2])

ggplot(pred_model_honeys1, aes(x=mc, y=Estimate))+
  geom_line(color="blue")+ # including group effects
  geom_line(data=pred_model_honeys2, aes(x=mc, y=Estimate), lty=2, color="red")+ # without group effects
  geom_point(data=all_honey, aes(x=mc, y=aw))+
      geom_line(data=pred_unpooled, aes(x=mc, y=Estimate), color="black", lty=3)+
  facet_wrap(~origin, ncol=6)+
  labs(x=expression("centered m"[c]), y=expression("a"[w]))

# code for Table 4 using kable

honeys_var2 <- summary(honeys_partial2)
honeys_var3 <- rbind(data.frame(honeys_var2$fixed), data.frame(honeys_var2$random$origin),
                     data.frame(honeys_var2$spec_pars))
rownames(honeys_var3) <- c("intercept $\\beta_0$", "slope $\\beta_1$","$\\sigma_b$", "$\\sigma_r$")
colnames(honeys_var3) <- c("mean","SE", "lower bound", "upper bound")
honeys_var3[1:4,1:4] %>% 
    rownames_to_column(var = "parameter") %>% kbl(booktabs=T, escape=F, digits = c(2,3,3,3,3)) %>% kable_styling(position="center", full_width = F)


# code for Figure 12

newvary=expand.grid(origin=1:29, mc=seq(from = -5, to = 12, by = 0.1))

pred_honey <- cbind(newvary, predict(honeys_partial2, newdata=newvary,re_formula = NA)[,-2])
names(pred_honey) <- c("origin", "mc", "aw", "lower", "upper")
pred_honey$origin <- as.factor(pred_honey$origin)

fit_honey <- cbind(newvary, fitted(honeys_partial2, newdata=newvary,re_formula = NA)[,-2])
names(fit_honey) <- c("origin", "mc", "aw", "lower", "upper")
fit_honey$origin <- as.factor(fit_honey$origin)

modelplot1 <- all_honey_plot <- ggplot(all_honey, aes(x=m, y=aw), )+
  geom_point(aes(color=origin), shape=21)+
   labs(x="moisture content m (%)", y=expression("a"[w]), subtitle = "B")+
  geom_point(data=model_honey, aes(x=m, y=aw))+
    theme(legend.position = "none")


modelplot2 <- model_honey %>% ggplot(aes(x=mc, y=aw))+
  geom_point()+
  geom_line(data=pred_honey, aes(x=mc, y=aw), color="black", size=1.2)+
  geom_ribbon(data = pred_honey, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = .6)+
  geom_ribbon(data = fit_honey, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .3)+
  labs(x=expression("centered m"[c]), y=expression("a"[w]), subtitle = "A")

modelplot2 + modelplot1

# code for Figure 13


# add new origin
newvary2 <- expand_grid(mc=seq(from = -5, to = 12, by = 0.1), origin=c("29","new origin"))

#new level, partial pooling, based on data uncertainty (alternative is to use sample_new_levels="gaussian" based on model's parameters)
predavg_new <- cbind(newvary2, predict(honeys_partial2, newdata=newvary2, re_formula=NULL, allow_new_levels=TRUE, sample_new_levels="gaussian")[,-2]) 

new_origin_plot <- ggplot(data=predavg_new[which(predavg_new$origin %in% c("new origin")),], aes(x=mc, y=Estimate))+
  geom_ribbon( aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`), fill="lightblue", alpha=0.5)+
  geom_line(aes(y=Estimate))+
  geom_ribbon(data=predavg_new[which(predavg_new$origin %in% c("29")),], 
              aes(y=Estimate,ymin = `Q2.5`, ymax = `Q97.5`), fill="red", alpha=0.5)+
  geom_point(data=beckh, aes(x=mc, y=aw))+
  labs(x=expression("centered m"[c]), y=expression("a"[w]))

new_origin_plot


# Code for the Supplemental material


# code for Figure S1:

all_honey %>% ggplot(aes(x=m, y=aw))+
  geom_point(aes(color=origin))+
  facet_wrap(~origin, ncol=6)+
  labs(x="moisture content (%)", y= "water activity")+
  theme(legend.position = "none")

# code for Figure S2

plot(honeys_pooled, widths=c(1,2))

# code for Figure S3, produced with bayesplot

mcmc_rank_overlay(honeys_pooled_post,pars=vars(b_Intercept:sigma) )

# Code for Figure S4, using bayesplot 

mcmc_acf(honeys_pooled_post, pars=vars(b_Intercept:sigma))

# raw brms output showing Rhat and n_eff:
print(honeys_pooled)

# code for Figure S5:

plot_prior_mean <- ggplot(data.frame((x=c(0.3,0.8))), aes(x=x))+
  stat_function(fun=dnorm, args=list(0.55,0.1),color="black", lty=2)+
  stat_function(fun=dnorm, args=list(0.552,0.001))+
  labs(x=expression(mu), y="density", subtitle = "A")

plot_prior_sigma <- ggplot(data.frame((xxxx=c(0,0.1))), aes(x=xxxx))+
  stat_function(fun=dexp, args=list(1),color="black", lty=2)+
  stat_function(fun=dnorm, args=list(0.05,0.001))+
  labs(x=expression(sigma), y="density", subtitle = "B")

plot_prior_mean + plot_prior_sigma

# Code for Table S1:

honeys_offsets <- ranef(honeys_partial)$origin %>% as_tibble(rownames="origin")
honeys_offsets <- honeys_offsets %>% select(-(c(Q2.5.Intercept,Q97.5.Intercept))) %>% rename(Estimate=Estimate.Intercept, SE=Est.Error.Intercept)
knitr::kable(honeys_offsets, booktabs=T, escape=F, digits = c(2,2,3)) %>% kable_styling(position="center", full_width = F)

# Code for Figure S6:
ppc1_plot <- pp_check(honeys_partial, type ="dens_overlay")+ggplot2::labs(subtitle = "A")
ppc2_plot <- pp_check(honeys_partial, type="stat")+ggplot2::labs(subtitle = "B")+ theme(axis.text.x = element_text(angle=70, hjust=1))
ppc3_plot <- pp_check(honeys_partial, type="ecdf_overlay")+ggplot2::labs(subtitle = "C")
ppc1_plot+ppc2_plot+ppc3_plot

# Code for Figure S7:
honeys_prior <- brm(data=all_honey, 
                    family = gaussian,
                    aw ~ 1 + mc,
                    prior=c(prior(normal(0.55,0.1), class=Intercept),
                            prior(normal(0.05,0.01), class=b),
                            prior(exponential(1), class=sigma)),
                    iter=4000, warmup=2000, chains=4, refresh=0,
                    sample_prior = T,
                    file=here("fits","honeys_prior"))
prior_honeys <- prior_draws(honeys_prior)

prior_honeys %>% slice_sample(n=100) %>% rownames_to_column("draw") %>% 
  expand(nesting(draw, Intercept, b),
         a = c(-5, 5)) %>% 
  mutate(d = Intercept + b * a) %>% 
  
  ggplot(aes(x = a, y = d)) +
  geom_line(aes(group = draw),
            color = "firebrick", alpha = .4) +
  labs(x = expression("m"[c]), y = expression("a"[w]))

# Code for Figure S8, first regression without and with centering, then the plot:

# completely pooled with centered predictor
honeys_pooled2 <- brm(data=all_honey, family = gaussian,
                      aw ~ 1 + mc,
                      prior(normal(0.55,0.1), class=Intercept),
                      prior(normal(0.5,0.1), class=b),
                      prior(exponential(1), class=sigma),
                      iter=4000, warmup=2000, chains=4, refresh=0,
                      file=here("fits","honeys_pooled2"))
print(honeys_pooled2, digits=3)
honeys_pooled2_post <- as_draws_df(honeys_pooled2)

# completely pooled predictor, no centering
honeys_pooled2b <- brm(data=all_honey, family = gaussian,
                       aw ~ 1 + m,
                       prior(normal(0.2,0.1), class=Intercept),
                       prior(normal(0.5,0.1), class=b),
                       prior(exponential(1), class=sigma),
                       iter=4000, warmup=2000, chains=4, refresh=0,
                       file=here("fits","honeys_pooled2b"))
print(honeys_pooled2b, digits=3)
honeys_pooled2b_post <- as_draws_df(honeys_pooled2b)
# producing the actual plots with GGally:
cor_db_2b <- dplyr::select(honeys_pooled2b_post,b_Intercept:sigma)
cor_db_2b <- setNames(cor_db_2b, c(expression(beta[0]), expression(beta[1]),expression(sigma[e])))

cor_db_2b_plot <-cor_db_2b  %>% 
  ggpairs(diag=list(continuous="densityDiag"),
          mapping=aes(fill="red"),
          upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
          labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 10, color = "red"),
        strip.text.y = element_text(size = 10, color = "red"))+
  theme(axis.text.x = element_text(angle=70, hjust=1))

cor_db_2 <- dplyr::select(honeys_pooled2_post,b_Intercept:sigma)
cor_db_2 <- setNames(cor_db_2, c(expression(beta[0]), expression(beta[1]),expression(sigma[e])))

cor_db_2_plot <-cor_db_2  %>% 
  ggpairs(diag=list(continuous="densityDiag"),
          mapping=aes(fill="red"),
          upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
          labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 10, color = "red"),
        strip.text.y = element_text(size = 10, color = "red"))+
  theme(axis.text.x = element_text(angle=70, hjust=1))

cor_db_2b_plot
cor_db_2_plot

# Code for Figure S9

p_c <- bind_rows(
  model_nested_honeys$post[[1]],
  model_nested_honeys$post[[2]],      
  model_nested_honeys$post[[3]],
  model_nested_honeys$post[[4]],
  model_nested_honeys$post[[5]],
  model_nested_honeys$post[[6]],
  model_nested_honeys$post[[7]],
  model_nested_honeys$post[[8]],
  model_nested_honeys$post[[9]],
  model_nested_honeys$post[[10]],
  model_nested_honeys$post[[11]],
  model_nested_honeys$post[[12]],
  model_nested_honeys$post[[13]],
  model_nested_honeys$post[[14]],
  model_nested_honeys$post[[15]],
  model_nested_honeys$post[[16]],
  model_nested_honeys$post[[17]],
  model_nested_honeys$post[[18]],
  model_nested_honeys$post[[19]],
  model_nested_honeys$post[[20]],
  model_nested_honeys$post[[21]],
  model_nested_honeys$post[[22]],
  model_nested_honeys$post[[23]],
  model_nested_honeys$post[[24]],
  model_nested_honeys$post[[25]],
  model_nested_honeys$post[[26]],
  model_nested_honeys$post[[27]],
  model_nested_honeys$post[[28]],
  model_nested_honeys$post[[29]]
)
iter <- 8000 

p_c <- 
  p_c %>% 
  mutate(origin = rep(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14", "15","16","17","18","19","20","21","22","23", "24","25","26","27","28","29"),each = iter))%>%  mutate(origin=fct_relevel(origin, "1","2","3","4", "5", "6", "7", "8", "9", "10",  "11", "12",  "13", "14",  "15", "16", "17", "18","19","20","21","22","23", "24","25","26","27","28","29"))

p_intercept_plot <- p_c %>% 
  ggplot(aes(x = b_Intercept, y = origin)) +
  stat_halfeye(fill = "green4", alpha=0.5,point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("intercept ",beta[0])), y="origin")+
  coord_cartesian(xlim=c(0.25,1.0))

p_slope_plot <- p_c %>% 
  ggplot(aes(x = b_mc, y = origin)) +
  stat_halfeye(fill = "green4", alpha=0.5,point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("slope ",beta[1])), y="origin")+
  coord_cartesian(xlim=c(-0.1,0.1))

p_intercept_plot+p_slope_plot

# Code for Figure S9:

ppc_plot1a <- pp_check(honeys_partial2, type = "dens_overlay")+labs(subtitle = "A")
ppc_plot1b <- pp_check(honeys_partial2, type="ecdf_overlay") + labs(subtitle = "B")
ppc_plot1a+ppc_plot1b

# Code for Table S2:

honeys_offsets2 <- ranef(honeys_partial2)$origin %>% as_tibble(rownames="origin")
honeys_offsets2 <- honeys_offsets2 %>% select(-(c(Q2.5.Intercept,Q97.5.Intercept))) %>% rename(Estimate=Estimate.Intercept, SE=Est.Error.Intercept)
knitr::kable(honeys_offsets2, booktabs=T, escape=F, digits = c(2,3,3)) %>% kable_styling(position="center", full_width = F)

# Overview of all R packages used

sessionInfo()

