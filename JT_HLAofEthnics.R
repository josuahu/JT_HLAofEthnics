#Author: Jozsef Toth
#Title: Predict ethnicity by patient's HLA set (Capstone)
#Date: July 4, 2024

#For this project I will use the IPD-IMGT/HLA database as source. 
#Database URL : https://www.ebi.ac.uk/ipd/imgt/hla/
#API help page : https://www.ebi.ac.uk/ipd/imgt/hla/about/help/api/
#License type: Creative Commons Attribution-NoDerivs Licence
#  https://www.ebi.ac.uk/ipd/imgt/hla/licence/
#  Barker DJ, Maccari G, Georgiou X, Cooper MA, Flicek P, Robinson J, Marsh SGE
#  The IPD-IMGT/HLA Database
#  Nucleic Acids Research (2023) 51:D1053-60

#Install required libraries:
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(rpart)) install.packages("rpart", repos = "http://cran.us.r-project.org")
if(!require(jsonlite)) install.packages("jsonlite", repos = "http://cran.us.r-project.org")
if(!require(xtable)) install.packages("xtable", repos = "http://cran.us.r-project.org")
if(!require(rpart)) install.packages("rpart", repos = "http://cran.us.r-project.org")
if(!require(rattle)) install.packages("rattle", repos = "http://cran.us.r-project.org")

###############################################################
#Download data from IPD-IMGT/HLA
#Download time is pretty long, about 10 hours,
#Please use the prepared imgt_cell_dat_raw.csv from Github
###############################################################

#Helper function to handle null values
ConvertNull <- function(obj) {
  ifelse(is.null(obj),NA,obj)
}

#Helper function to convert allele digits to 4 digit format
convert_4_digit_allele <- function(s) {
  result <- NA
  if(!is.na(s)) {
    s_parts <- unlist(strsplit(s, split = ":"))
    s1 <- str_trim(s_parts[1])
    s2 <- str_trim(s_parts[2])
    if(!is.na(s1) & !is.na(str_match(s1,"[0-9]{2,3}")) & s1==str_match(s1,"[0-9]{2,3}") & 
       !is.na(s2) & !is.na(str_match(s2,"[0-9]{2,3}")) &  s2==str_match(s2,"[0-9]{2,3}")) {
      result <- paste0(s1,":",s2)
    }
  }
  return(result)
}

#Target raw data filename for later use
rawdata_file <- "imgt_cell_dat_raw.csv"

#Only download if not exists, download can take hours...
if(!file.exists(rawdata_file)) {
  #Get Cell IDs :
  i<-1
  error_code <- 0
  baseurl<-"https://www.ebi.ac.uk/cgi-bin/ipd/api/cell"
  nexturl<-"?limit=1000"
  id_list <- c()
  while (error_code==0) {
    print(paste0("Downloading index list ",i))
    get_idlist_json <- fromJSON(paste0(baseurl,nexturl), flatten = TRUE)
    if(length(get_idlist_json$data)==2) {
      id_list <- c(id_list,get_idlist_json$data$cellid)
      if(!is.null(get_idlist_json$meta$`next`)) {
        nexturl <- get_idlist_json$meta$`next`
      }else {
        error_code <- -1
      }
    }else {
      error_code <- -2
    }
    i <- i+1
  }
  length(id_list)

  #Download cell data for each cell ID :
  rawdat <- map_df(1:length(id_list), function(i) {
    id <- id_list[i]
    print(paste0("Downloading cell data ",i,"/",length(id_list)))
    cell_url <- paste0("https://www.ebi.ac.uk/cgi-bin/ipd/api/cell/",id)
    cell_data <- fromJSON(cell_url, flatten = TRUE)
    i<-i+1
    data.frame(cellid = cell_data$cellid,
               ancestry = ConvertNull(cell_data$source$ancestry),
               HLA_A_1 = ConvertNull(cell_data$typing$A[1]),
               HLA_A_2 = ConvertNull(cell_data$typing$A[2]),
               HLA_B_1 = ConvertNull(cell_data$typing$B[1]),
               HLA_B_2 = ConvertNull(cell_data$typing$B[2]),
               HLA_C_1 = ConvertNull(cell_data$typing$C[1]),
               HLA_C_2 = ConvertNull(cell_data$typing$C[2]))
  })
  
  #Save result file for later use
  write.csv(rawdat,rawdata_file,row.names=FALSE)  
  
  #Remove variables
  rm(id_list)
  rm(rawdat)
}

#Reload the (pre) downloaded rawdata file
rawdat <- read.csv(rawdata_file,header=TRUE)

########### Data cleaning ###########

#Convert to 4 digit allele format
rawdat[,3:8] <- apply(rawdat[,3:8],c(1,2),convert_4_digit_allele)

#In the next step, we filter the cells to have the 6 alleles and convert ancestry texts to 
#descents. Only samples (cells) with valid descent are kept.

#Generate cleared data (all allele is fine, descent is valid):
data_cleared <- rawdat %>%
  filter(!is.na(HLA_A_1) & !is.na(HLA_A_2) &
           !is.na(HLA_B_1) & !is.na(HLA_B_2) &
           !is.na(HLA_C_1) & !is.na(HLA_C_2)) %>%
  mutate(Descent = NA,
         Descent = ifelse(str_starts(ancestry,"Admixed"),"Admixed",Descent),
         Descent = ifelse(str_starts(ancestry,"African"),"African",Descent),
         Descent = ifelse(str_starts(ancestry,"Asian"),"Asian",Descent),
         Descent = ifelse(str_starts(ancestry,"European"),"European",Descent),
         Descent = ifelse(str_starts(ancestry,"Greater Middle Eastern"),"Greater Middle Eastern",Descent),
         Descent = ifelse(str_starts(ancestry,"Hispanic or Latin American"),"Hispanic or Latin American",Descent),
         Descent = ifelse(str_starts(ancestry,"Native American"),"Native American",Descent),
         Descent = ifelse(str_starts(ancestry,"Oceanian"),"Oceanian",Descent),
         Descent = ifelse(str_starts(ancestry,"Aboriginal Australian"),"Aboriginal Australian",Descent),
         Descent = ifelse(str_starts(ancestry,"Undefined"),"Undefined",Descent)) %>%
  filter(!is.na(Descent)) %>%
  select(cellid,Descent,HLA_A_1,HLA_A_2,HLA_B_1,HLA_B_2,HLA_C_1,HLA_C_2)

#Add allele prefix before the HLA digits
data_cleared$HLA_A_1=paste0("HLA-A*",data_cleared$HLA_A_1)
data_cleared$HLA_A_2=paste0("HLA-A*",data_cleared$HLA_A_2)
data_cleared$HLA_B_1=paste0("HLA-B*",data_cleared$HLA_B_1)
data_cleared$HLA_B_2=paste0("HLA-B*",data_cleared$HLA_B_2)
data_cleared$HLA_C_1=paste0("HLA-C*",data_cleared$HLA_C_1)
data_cleared$HLA_C_2=paste0("HLA-C*",data_cleared$HLA_C_2)

########### Data exploration ###########

nrow(data_cleared)
head(data_cleared)

#The descent distribution is the following :
descents_distribution <- data_cleared %>%
  group_by(Descent) %>%
  summarize(N=n()) %>%
  arrange(desc(N))
head(descents_distribution)

#Unfortunately we have a lot of unknown descents. We will focus on the European and Asian 
#ones because they have around the similar counts.

#Filter data to European and Asian
dat <- data_cleared %>% 
  filter(Descent %in% c("Asian","European"))
nrow(dat)

#Before we continue, split the data into train and test set in the ratio of 90% vs 10%
set.seed(1)
test_index = createDataPartition(dat$Descent, times = 1, p = 0.1, list = FALSE)
train_set <- dat[-test_index,]
test_set <- dat[test_index,]

nrow(test_set)
nrow(train_set)

#Create long format of the train set:
train_set_long <- train_set %>% 
  pivot_longer(cols=-c("cellid","Descent"),names_to="HLAGroup",values_to="HLAValue")

#Check HLA frequency (use distinct because homozigotes)
hla_freq <- train_set_long %>%
  select(cellid,HLAValue) %>%
  distinct() %>%
  group_by(HLAValue) %>%
  summarise(N=n(),Percentage=100*N/nrow(train_set)) %>%
  arrange(desc(N))
head(hla_freq)

#We can see that there are common alleles. The top 100 HLA frequency looks like this :
head(hla_freq,100) %>%
  ggplot(aes(reorder(HLAValue, -Percentage),Percentage)) +
  geom_bar(stat="identity") + xlab("TOP100 HLA alleles") +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())

#Compare top HLA frequencies for the descents
asian_count <- sum(train_set$Descent=="Asian")
hla_freq_asian <- train_set_long %>%
  filter(Descent=="Asian") %>%
  select(cellid,HLAValue) %>%
  distinct() %>%
  group_by(HLAValue) %>%
  summarise(N=n(),Percentage=100*N/asian_count) %>%
  arrange(desc(N)) %>%
  head(5)
hla_freq_asian

european_count <- sum(train_set$Descent=="European")
hla_freq_european <- train_set_long %>%
  filter(Descent=="European") %>%
  select(cellid,HLAValue) %>%
  distinct() %>%
  group_by(HLAValue) %>%
  summarise(N=n(),Percentage=100*N/european_count) %>%
  arrange(desc(N)) %>%
  head(5)
hla_freq_european

########### Methods ###########

#Preparation with One-Hot Encoding
#In the first three model we will apply the classification tree, the KNN and the GLM methods. 
#The test set contains around 2000 HLAs what is a lot of HLA combinations for the subject's 
#6 MHCI allele. Converting the HLA names into numbers for the target methods method can lead 
#to problems, because the order of the alleles may interferes the result of the algorithm. 
#Having an HLA is basically a categorical data, so we can use One-Hot Encoding to deal with it. 
#For this, we have to convert our data set into a much wider format, where every possible HLA 
#has a column and the the value is 1 if that subject has the HLA and the value is 0 if the 
#subject does not have that HLA. For practical reasons we will not convert 2000 HLA to columns, 
#we will use only the top 100 HLA in the frequency list. This way the model is viable for much 
#larger populations as well, when we have much more HLAs.

#Select top 100 HLA by frequency
top100_hla_list <- head(hla_freq,100)$HLAValue

#One-Hot Encoding for the train set's top 100 HLA:
train_set_ohe <- train_set_long %>%
  filter(HLAValue %in% top100_hla_list) %>%
  select(cellid,Descent,HLAValue) %>% 
  distinct() %>%
  mutate(HaveHLA=1) %>%
  pivot_wider(names_from = HLAValue, values_from = HaveHLA) %>%
  select(-cellid) %>%
  mutate(Descent=factor(Descent,levels=c("European","Asian")))

#Convert NA values to zeroes
train_set_ohe[is.na(train_set_ohe)] <- 0

#Convert columns names for train compatible format
colnames(train_set_ohe) <- make.names(colnames(train_set_ohe))

#Our data looks like this after the One-Hot encoding (showing only the first rows and columns) : 
head(train_set_ohe[,1:6])

#Create long format of the train set:
test_set_long <- train_set %>% 
  pivot_longer(cols=-c("cellid","Descent"),names_to="HLAGroup",values_to="HLAValue")

#One-Hot Encoding for the top 100 HLA Test set
test_set_ohe <- test_set_long %>%
  filter(HLAValue %in% top100_hla_list) %>%
  select(cellid,Descent,HLAValue) %>% 
  distinct() %>%
  mutate(HaveHLA=1) %>%
  pivot_wider(names_from = HLAValue, values_from = HaveHLA) %>%
  select(-cellid) %>%
  mutate(Descent=factor(Descent,levels=c("European","Asian")))

#Convert NA values to zeroes
test_set_ohe[is.na(test_set_ohe)] <- 0

#Convert columns names for train compatible format
colnames(test_set_ohe) <- make.names(colnames(test_set_ohe))

#Now we are ready to apply some traditional models :

#Classification tree :
fit_class <- rpart(Descent ~ ., data = train_set_ohe)
rpart.plot::prp(fit_class,uniform = TRUE, compress = TRUE, branch = .2)

#Classification tree performance on the test set
y_hat_class <- predict(fit_class, test_set_ohe, type = "class")
mean(y_hat_class==test_set_ohe$Descent)

#The accuracy is not bad. We can see some alleles from the TOP5 lists because the common alleles
# play a role in determining descents. The top level allele is HLA-C*07:01, this means that if the
# subject has this allele we predict European. Let's check the the ratios for this allele for the
# test set :
onehla_freq <- test_set_long %>%
  filter(HLAValue=="HLA-C*07:01") %>%
  select(Descent,cellid,HLAValue) %>%
  distinct() %>%
  group_by(Descent) %>%
  summarise(N=n()) %>%
  mutate(Percentage=ifelse(Descent=="Asian",100*N/asian_count,100*N/european_count))

head(onehla_freq)

#There is a big difference between the 2 ratios. 
#This allele is really more possible for Europeans.

#For the KNN method we will try different k values from 3 to 29.

#Knn model: 
knn_fit <- train(Descent ~ ., method = "knn", data = train_set_ohe,
                 tuneGrid = data.frame(k = seq(3,29,by = 2))) 
knn_fit$bestTune

#Tuning parameters:
plot(knn_fit)

#Accuracy on the test set :
y_hat_knn <- predict(knn_fit, test_set_ohe, type = "raw")
mean(y_hat_knn==test_set_ohe$Descent)

#The KNN method has better accuracy. The optimal k parameter is surprisingly large, but the
#difference between the k values from 15 to 31 is not too high.

#GLM Model
glm_fit <- train(Descent ~ ., method = "glm", data = train_set_ohe)
y_hat_glm <- predict(glm_fit, test_set_ohe, type = "raw")
mean(y_hat_glm==test_set_ohe$Descent)

#The GLM method performs pretty well.

#Relative Frequency model
# In this part we create a new model. For this model we will calculate relative frequencies for
# all HLA. Relative frequency is defined to every allele independently from others as the 
# possibility to being Asian or being European descent.

#Calculate percentages of being EU or AS for each HLA
hla_freq_descent <- train_set_long %>% 
  distinct() %>%
  mutate(EU = ifelse(Descent=="European",1,0),
         AS = ifelse(Descent=="Asian",1,0)) %>%
  group_by(HLAValue) %>%
  summarize(N=n(),EU=sum(EU,na.rm = TRUE)/N,AS=sum(AS,na.rm = TRUE)/N) %>%
  arrange(desc(N))
head(hla_freq_descent)

# This means if we pick up HLA-A*02:01 that it is European descent in 72% and Asian in 28%. 
# If we pick up HLA-A*24:02, it is Asian descent in 64% and European in 36%. And so on for the 
# other alleles. By using this informations we can can build and algorith to summarize this 
# values for the subjects.

#Method performance on the test set:
test_set_long %>% 
  left_join(hla_freq_descent,by=c("HLAValue"))  %>%
  group_by(cellid,Descent) %>%
  summarize(EU=prod(EU,na.rm = TRUE),AS=prod(AS,na.rm = TRUE)) %>%
  mutate(predicted_descent = ifelse(EU>AS,"European","Asian"),
         ok=ifelse(predicted_descent==Descent,1,0)) %>%
  ungroup() %>%
  summarize(Method="Frequency model 1 (test set)",accuracy=sum(ok,na.rm = TRUE)/n()) %>%
  knitr::kable()

#Optimization for minimal ratio because zero and 100% values. Let's try to increase 0% values 
#and decrease 100% values with a constant.
result_optimization <- map_df(seq(0,0.05,by=0.001),function(min_ratio) {
  test_set_long %>% 
    left_join(hla_freq_descent,by=c("HLAValue"))  %>%
    mutate(EU=pmax(EU,min_ratio),AS=pmax(AS,min_ratio)) %>%
    mutate(EU=pmin(EU,1-min_ratio),AS=pmin(AS,1-min_ratio)) %>%
    group_by(cellid,Descent) %>%
    summarize(EU=prod(EU,na.rm = TRUE),AS=prod(AS,na.rm = TRUE)) %>%
    mutate(predicted_descent = ifelse(EU>AS,"European","Asian"),
           ok=ifelse(predicted_descent==Descent,1,0)) %>%
    ungroup() %>%
    summarize(N=n(),min_ratio=min_ratio,accuracy=sum(ok,na.rm = TRUE)/N)
})
result_optimization %>%
  ggplot(aes(min_ratio,accuracy)) +
  geom_point()

#As we see the the optimization does not really help to improve our previous accuracy. 
#So we stay with the previous version.

#This result is based on a relatively small test group, so better to do a 10-fold cross 
#validation to get the real accuracy. We do this only for our selected method.

#Manually create 10 fold for the cross validation:
set.seed(1)
split_index <- createFolds(dat$Descent, k = 10, list = FALSE, returnTrain = FALSE)
table(split_index)

cv_result <- map_df(seq(1,10,1), function(fold_index) {
  index_test <- which(split_index==fold_index)
  train_set <- dat[-index_test,]
  test_set <- dat[index_test,]
  #Create long format of the train set:
  #Determine HLA frequencies
  hla_freq_descent <- train_set %>% 
    pivot_longer(cols=-c("cellid","Descent"),names_to="HLAGroup",values_to="HLAValue") %>% 
    distinct() %>%
    mutate(EU = ifelse(Descent=="European",1,0),
           AS = ifelse(Descent=="Asian",1,0)) %>%
    group_by(HLAValue) %>%
    summarize(N=n(),EU=sum(EU,na.rm = TRUE)/N,AS=sum(AS,na.rm = TRUE)/N) 
  #Calculate accuracy
  test_set %>% 
    pivot_longer(cols=-c("cellid","Descent"),names_to="HLAGroup",values_to="HLAValue") %>% 
    left_join(hla_freq_descent,by=c("HLAValue"))  %>%
    group_by(cellid,Descent) %>%
    summarize(EU=prod(EU,na.rm = TRUE),AS=prod(AS,na.rm = TRUE)) %>%
    mutate(predicted_descent = ifelse(EU>AS,"European","Asian"),
           ok=ifelse(predicted_descent==Descent,1,0)) %>%
    ungroup() %>%
    summarize(Method="Frequency model",fold_index=fold_index,accuracy=sum(ok,na.rm = TRUE)/n()) 
})

real_accuracy <- mean(cv_result$accuracy)
real_accuracy

########### Results / Conclusion ###########
# The cross validated accuracy is good. This match our goal to have at least 90% for this two 
# descents. The developed Relative Frequency model is fast and scalable for subject size. 
# For this two descents we was able make a pretty accurate method but presumably the model 
# performance will drop if we try to add more descents. For example adding descents which are naturally lives close to each other will introduce more similarities in HLA sets and making harder to distinguish them. In any case, population-based or HLA-tailored vaccines may have a place in the vaccine planning process.


########### ##References ###########
# *The data is downloaded from : *
#     IPD-IMGT/HLA database \
#     https://www.ebi.ac.uk/ipd/imgt/hla/ \
#     https://www.ebi.ac.uk/ipd/imgt/hla/licence/
# *You can read more about MHC and HLA here: *
#     Major histocompatibility complex Wikipedia page \
#     https://en.wikipedia.org/wiki/Major_histocompatibility_complex 
