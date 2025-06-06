rm(list = ls())
require("ggplot2")
load("data/dfm_filtered_resolved_unigram_ManualvocFilt_sentiment_accronym.rda")

sentiment = dat@docvars$sentiment
bw = 2
n_obs = sum(!is.na(dat@docvars$sentiment))

g <- ggplot(dat@docvars, aes(sentiment))  + 
  geom_histogram( colour = "gray50") + geom_vline(xintercept = quantile(sentiment,0.33333),col = "red",size=2)+
  theme_bw() + xlab("Sentiment") + ylab("Count")

pdf("figures/sentiment.pdf",
    height = 6,
    width = 8,
    paper = "special")
g
dev.off()
