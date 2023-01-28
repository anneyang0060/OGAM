setwd('.../OGAM/datasets/flight')

for(year in 1988:2005){
  
  file1 <-paste0(year,'_gam_part1.Rdata')
  file2 <-paste0(year,'_gam_part2.Rdata')
  load(file1); df1 <- df
  load(file2); df2 <- df
  file.remove(file1)
  file.remove(file12
  
  if(year<=1996){
    df <- rbind(df1, df2)
  }else{
    file3 <-paste0(year,'_gam_part3.Rdata')
    load(file3); df3 <- df
    df <- rbind(df1, df2)
    df <- rbind(df, df3)
    file.remove(file3)
  }
  
  save(df,file = paste0(year,'_gam.Rdata'))
}
