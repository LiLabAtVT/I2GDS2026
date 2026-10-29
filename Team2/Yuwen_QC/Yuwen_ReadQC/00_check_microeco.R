# Inspect how microeco stores FAPROTAX function groups
library(microeco)
src <- deparse(trans_func$public_methods$cal_spe_func)
cat(grep("func_group_list|func_groups", src, value = TRUE), sep = "\n")
data(prok_func_FAPROTAX)
str(prok_func_FAPROTAX, max.level = 1)
if (!is.null(prok_func_FAPROTAX$func_groups)) str(prok_func_FAPROTAX$func_groups, max.level = 2)
