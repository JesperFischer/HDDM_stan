#Predictive checks?

 #RLDMM 
library(posterior)

n_check = 20

rts = as_draws_df(fit1$draws()) %>% select(matches("out\\[\\d+,1\\]")) %>% slice(1:n_check) %>% mutate(draw = 1:n_check) %>% pivot_longer(-draw, values_to = "predicted_rt",names_to = "trial")%>%
  mutate(trial = gsub(".*\\[(\\d+).*\\]", "\\1", trial))

expect = as_draws_df(fit1$draws()) %>% select(matches("expect\\[\\d+\\]")) %>% slice(1:n_check) %>% mutate(draw = 1:n_check) %>% dplyr::select(-"expect[201]") %>% pivot_longer(-draw, values_to = "predicted_expectation",names_to = "trial")%>%
  mutate(trial = gsub(".*\\[(\\d+).*\\]", "\\1", trial))

expect %>% ggplot()+geom_line(aes(x = as.numeric(trial), y = predicted_expectation, group = as.factor(draw)))+theme_classic()

rts %>% ggplot()+geom_density(aes(x = predicted_rt, group = draw))+geom_density(data = resp, aes(x = q), col = "red")+theme_classic()
```


# ```{r}
# rts = as_draws_df(fit1$draws()) %>% select(matches("out\\[\\d+,1\\]")) %>% mutate(draw = 1:nrow(.)) %>% pivot_longer(-draw, values_to = "predicted_rt",names_to = "trial")%>%
#   mutate(trial = gsub(".*\\[(\\d+).*\\]", "\\1", trial))
# 
# expect = as_draws_df(fit1$draws()) %>% select(matches("expect\\[\\d+\\]")) %>% 
#   mutate(draw = 1:nrow(.)) %>% 
#   dplyr::select(-"expect[201]") %>% 
#   pivot_longer(-draw, values_to = "predicted_expectation",names_to = "trial")%>%
#   
#   mutate(trial = gsub(".*\\[(\\d+).*\\]", "\\1", trial))
# 
# deltat = as_draws_df(fit1$draws()) %>% select(matches("deltat\\[\\d+\\]")) %>% 
#   mutate(draw = 1:nrow(.)) %>% 
#   pivot_longer(-draw, values_to = "deltat",names_to = "trial")%>%
#   mutate(trial = gsub(".*\\[(\\d+).*\\]", "\\1", trial))
# 
# 
# 
# inner_join(deltat,rts) %>% group_by(trial) %>% 
#   summarize(sd_expect = sd(deltat),
#             deltat = median(deltat),
#             sd_rt = sd(predicted_rt),
#             predicted_rt = median(predicted_rt)) %>% 
#   ggplot(aes(x = deltat, y = predicted_rt)) +geom_point()+geom_smooth()+theme_classic()+
#   scale_x_continuous(breaks = seq(-2,2,by = 0.5))
# 
# 
# inner_join(expect,rts) %>% group_by(trial) %>% 
#   summarize(sd_expect = sd(predicted_expectation),
#             predicted_expectation = median(predicted_expectation),
#             sd_rt = sd(predicted_rt),
#             predicted_rt = median(predicted_rt)) %>% 
#   ggplot(aes(x = predicted_expectation, y = predicted_rt)) +geom_point()+geom_smooth()+theme_classic()+
#   scale_x_continuous(breaks = seq(0,1,by = 0.1))
# 
# 
# 
# 
# 
# 
# resps = as_draws_df(fit1$draws()) %>% select(matches("out\\[\\d+,2\\]")) %>% mutate(draw = 1:nrow(.)) %>% pivot_longer(-draw, values_to = "predicted_resp",names_to = "trial")%>%
#   mutate(trial = gsub(".*\\[(\\d+).*\\]", "\\1", trial))
# 
# 
# inner_join(resps,rts) %>% group_by(trial) %>% 
#   summarize(predicted_resp = mean(predicted_resp),
#             sd_rt = sd(predicted_rt),
#             predicted_rt = median(predicted_rt)) %>% 
#   ggplot(aes(x = predicted_resp, y = predicted_rt)) +geom_point()+geom_smooth()+theme_classic()+
#   scale_x_continuous(breaks = seq(0,1,by = 0.1))
# 
# 
# inner_join(rts,expect) %>% group_by(trial) %>% 
#   summarize(sd_expect = sd(predicted_expectation),
#             predicted_expectation = median(predicted_expectation),
#             sd_rt = sd(predicted_rt),
#             predicted_rt = median(predicted_rt)) %>% 
#   ggplot(aes(x = predicted_expectation, y = predicted_rt)) +geom_pointrange(aes(ymin = predicted_rt-sd_rt,ymax = predicted_rt+sd_rt))+geom_smooth()+theme_classic()+
#   scale_x_continuous(breaks = seq(0,1,by = 0.1))
# 
# resp %>% ggplot(aes(x = expectation, y = q))+
#   geom_point()+
#   theme_classic()+
#   geom_smooth()
# ```
