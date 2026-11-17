obs <- list(S1=c(.018,.23), S2=c(.080,.12), S3=c(.025,.052), S4=c(.011,.07), S5=c(.071,.073))
a <- commandArgs(trailingOnly=TRUE)
for (st in a) {
 fs <- list.files(pattern=sprintf("^fallsize_%s_seed", st)); 
 d <- do.call(rbind, lapply(fs, function(f){x<-read.csv(f); x$study<-NULL; x}))
 o <- obs[[st]]; cat(sprintf("\n== %s  observed top %.3f bottom %.3f\n", st, o[1], o[2]))
 for (nm in unique(d$truth)) { s <- d[d$truth==nm,]; tp <- s$top[s$top>0]
  cat(sprintf("%-32s n %d | top fall %3.0f%% | size if any %.3f [%.3f, %.3f] | >= obs top %3.0f%% | bottom %.3f [%.3f, %.3f], >= obs bottom %3.0f%% | bot>=top given top %3.0f%%\n",
   nm, nrow(s), 100*mean(s$top>0), median(tp), quantile(tp,.1), quantile(tp,.9), 100*mean(s$top>=o[1]),
   median(s$bot), quantile(s$bot,.1), quantile(s$bot,.9), 100*mean(s$bot>=o[2]), 100*mean(s$bot[s$top>0]>=s$top[s$top>0]))) }
}
