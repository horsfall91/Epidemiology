library(R0)
library(tidyverse)
library(inflection)

cov_url <- "https://covid.ourworldindata.org/data/owid-covid-data.csv"
covid_data <- read.csv(cov_url)

countries <- unique(covid_data$location)
GT <- generation.time("gamma", c(5.2, 2.8))

country_vec <- c()
country_MLr0_vec <- c()
country_MLr0_vec_loCI <- c()
country_MLr0_vec_hiCI <- c()

country_EGr0_vec <- c()
country_EGr0_vec_loCI <- c()
country_EGr0_vec_hiCI <- c()
time_start <- Sys.time()

for (c in countries) {
  
  country <- covid_data[which(covid_data$location == c), ]
  country <- as.data.frame(country)
  
  tot <- country %>% select(total_cases, date)
  tot[is.na(tot)] <- 0
  tot <- tot[complete.cases(tot),]
  
  x <- c(1:length(tot$total_cases))
  x <- x[1:round(length(x)/2)]
  y <- tot$total_cases
  y <- y[1:round(length(y)/2)]
  inflec <- findipiterplot(x,y,0)
  inflec_ind <- round(inflec$first[5])
  
  
  i <- country %>% select(new_cases, date)
  i[is.na(i)] <- 0
  i <- i[complete.cases(i),]
  
  #create epidemi curve
  epid <- i %>% pull(new_cases)
  dates <- i %>% pull(date)
  names(epid) <- dates
  
  #set initial begin date and end date for sensitivity analysis
  begin_int <- min(which(data.frame(epid)[,'epid'] > 1))
  if (is.nan(inflec_ind)) {
    end_int <- min(which((data.frame(epid) %>% mutate(cmltv_i = cumsum(epid)))[,'cmltv_i'] > 1000))
  } else {
  end_int <- inflec_ind
  }
  
  # if (end_int-begin_int > 100) {
  #   end_int <- begin_int + inflec_ind
  # }
  
  begin_int <- 6
  if ((begin_int >= 11) & (length(epid)-end_int >=11)) {
    sense_mins <- 10
    sense_plus <- 10
  } else if (begin_int <= 11) {
    sense_mins <- begin_int-1
    sense_plus <- 20 - (begin_int-1)
  }
  
  #If insufficient cases, skip that country
  if (end_int == Inf) {
    next
  }
  
  #run sensitivity analysis to find data window with best fit of exponential growth
  t <- try(sensitivity.analysis(epid, GT, begin=(begin_int-sense_mins):(begin_int+sense_plus), end=(end_int-sense_mins):(end_int+sense_plus), est.method="EG", sa.type="time"))
  if ("try-error" %in% class(t)) {
    next
  } else {
    SI <- sensitivity.analysis(epid, GT, begin=(begin_int-sense_mins):(begin_int+sense_plus), end=(end_int-sense_mins):(end_int+sense_plus), est.method="EG", sa.type="time")
    max_R <- which.max(SI$df.clean$R)
    row_max_R <- SI$df.clean[max_R,]
    window_start <- which(names(epid) == row_max_R$Begin.dates)
    
    non_Z_vec <- which(epid>0, arr.ind=TRUE)
    window_start <- as.numeric(non_Z_vec[which(abs(non_Z_vec-window_start)==min(abs(non_Z_vec-window_start)))]) #get closest non-zero row to window_start
    
    window_end <- which(names(epid) == row_max_R$End.dates)
  }
  
  #Maximum Likelihood Method
  t <- try(est.R0.ML(epid, GT, begin=window_start*1, end=window_end*1, range=c(0.01,100), unknown.GT = FALSE),silent = T)
  if ("try-error" %in% class(t)) {
    ML.R0 = NA
  } else {
    ML.R0 <- est.R0.ML(epid, GT, begin=window_start*1, end=window_end*1, range=c(0.01,100), unknown.GT = FALSE)
  }
  
  #Exponential Growth Method
  t <- try(est.R0.EG(epid, GT, begin=window_start*1, end=window_end*1), silent = T)
  if ("try-error" %in% class(t)) {
    EG.R0 = NA
  } else { 
    EG.R0 <- est.R0.EG(epid, GT, begin=window_start*1, end=window_end*1)
  }
  
  #Assemble the Vectors
  country_vec <- append(country_vec,c)
  
  t <- try(country_MLr0_vec,ML.R0[1][[1]], silent = T)
  if ("try-error" %in% class(t)) {
    country_MLr0_vec <- append(country_MLr0_vec, 'NA')
    country_MLr0_vec_loCI <- append(country_MLr0_vec_loCI,'NA')
    country_MLr0_vec_hiCI <- append(country_MLr0_vec_hiCI,'NA')
  } else {  
    country_MLr0_vec <- append(country_MLr0_vec,ML.R0[1][[1]])
    country_MLr0_vec_loCI <- append(country_MLr0_vec_loCI,ML.R0[2][[1]][1])
    country_MLr0_vec_hiCI <- append(country_MLr0_vec_hiCI,ML.R0[2][[1]][2])
  }
  
  t <- try(country_EGr0_vec,EG.R0[1][[1]], silent = T)
  if ("try-error" %in% class(t)) {
    country_EGr0_vec <- append(country_EGr0_vec, 'NA')
    country_EGr0_vec_loCI <- append(country_EGr0_vec_loCI,'NA')
    country_EGr0_vec_hiCI <- append(country_EGr0_vec_hiCI,'NA')
  } else {  
    country_EGr0_vec <- append(country_EGr0_vec,EG.R0[1][[1]])
    country_EGr0_vec_loCI <- append(country_EGr0_vec_loCI,EG.R0[2][[1]][1])
    country_EGr0_vec_hiCI <- append(country_EGr0_vec_hiCI,EG.R0[2][[1]][2])
  }
  
}

R0.df <- data.frame(country_vec, 
                    
                    country_MLr0_vec, 
                    country_MLr0_vec_loCI, 
                    country_MLr0_vec_hiCI,
                    
                    country_EGr0_vec, 
                    country_EGr0_vec_loCI, 
                    country_EGr0_vec_hiCI
)
time_end <- Sys.time()
total_time <- time_end - time_start
print(total_time)
R0.df
write.csv(R0.df,"R0_dataframe_v5.csv", row.names = FALSE)
