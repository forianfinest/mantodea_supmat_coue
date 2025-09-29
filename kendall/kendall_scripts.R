##FOR CALCULATING KENDALL'S TAU

library(readr)
for_kendall_table <- read_delim("for_kendall_table.csv", 
    delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
        grouping_mark = "."), trim_ws = TRUE)
d=as.numeric(for_kendall_table$D)
pres_nat=as.numeric(for_kendall_table$Native)
pres_ali=as.numeric(for_kendall_table$Alien)
res_nat<-cor.test(pres_nat,d, method="kendall", exact = FALSE, conf.level = 0.95)
res_ali<-cor.test(pres_ali,d, method="kendall", exact = FALSE, conf.level = 0.95)

##FOR VISUALIZING THE RESULTS

res_nat
res_ali
