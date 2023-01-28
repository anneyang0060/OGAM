setwd('.../OGAM/datasets/credit')

file1 <- 'P2P_Macro_Data_processed_part1.Rdata'
file2 <- 'P2P_Macro_Data_processed_part2.Rdata'
load(file1); df1 <- df
load(file2); df2 <- df
file.remove(file1)
file.remove(file2)

df <- rbind(df1, df2)
save(df,file = 'P2P_Macro_Data_processed.Rdata')
