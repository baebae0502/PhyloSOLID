##Author: Zhirui Yang, Qing Yang
##Date: 2024-07-26
##Update: 2024/01/02
##Details: This is one script to calculate hFDR.


##### Time #####
library(stringr)
start_time <- Sys.time()
print(str_c("Start time: ", start_time))


##### library
library(optparse)
option_list = list(
    make_option(c("-i", "--inputfile"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/labeling/identifier.features.CI792_tumor.txt", help = "the name of inputfile - posterior raw data & reads information"),
    make_option(c("-o", "--outputfile"), type = "character", default = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_10xCI792/labeling/CI792.identifier.feature.sigFilter.txt", help = "The outputfile you can set")
)
parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)


inputfile <- as.character(opt$inputfile)
outputfile <- as.character(opt$outputfile)

# if (!file.exists(outputpath)) {
#   dir.create(outputpath)
#   print(paste("Path", outputpath, "created successfully"))
# } else {
#   print(paste("Path", outputpath, "already exists"))
# } 


#===========functions=================
library(pracma)
library(dplyr)
library(purrr)
library(ggplot2)


# function to get standard mutation signature matrix
get_default_labels <- function(choice) {
  if (!(choice %in% c("DNA", "RNA", "96", "192"))) {
    stop("Choice must be 'DNA'/'96' or 'RNA'/'192'")
  }
  if (choice == "DNA" | choice == "96") {
    mid_list <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    value <- 96
  } else {
    mid_list <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "A>C", "A>G", "A>T", "G>A", "G>C", "G>T")
    value <- 192
  }
  first <- c("A", "T", "C", "G")
  inner_bracket <- rep(rep(mid_list, each = 16), times = 1)
  outter_bracket <- expand.grid(first, first)
  result <- sapply(1:value, function(f) {
    paste0(outter_bracket[f %% 16 + 1, 1], "[", inner_bracket[f], "]", outter_bracket[f %% 16 + 1, 2])
  })
  return(result)
}


#=============handle data=============

# read features
df_candidate <- read.csv(inputfile,header = T,sep="\t")
if (!("identifier" %in% colnames(df_candidate))) {
  colnames(df_candidate)[colnames(df_candidate) == "X.identifier"] <- "identifier"
}

counts <- table(df_candidate$RNAMutationType)
sig <- as.data.frame(counts)
colnames(sig)<-c("MutationType","Count")

default_list=get_default_labels("RNA")
default_df=data.frame("MutationType" = default_list,"none" = rep(0,length((default_list))))
df_sigProfile<-left_join(default_df,sig,by="MutationType")

df_sigProfile$Count <- ifelse(is.na(df_sigProfile$Count), 0, df_sigProfile$Count)
df_sigProfile["none"]<-NULL
dim(df_sigProfile)
# [1] 192   2


##### test paired mutation type
# 定义碱基配对规则
complement_pairs <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

# 修正后的 get_pair 函数
get_pair <- function(mutation) {
  # 提取完整的 MutationType
  # 提取前缀、后缀和中间部分
  prefix <- substr(mutation, 1, 1) # 第一个碱基
  suffix <- substr(mutation, nchar(mutation), nchar(mutation)) # 最后一个碱基
  middle <- gsub(".*\\[|\\].*", "", mutation) # 中间部分（如 C>A）

  # 分解前缀、后缀和中间部分
  prefix_comp <- complement_pairs[prefix] # 配对前缀
  suffix_comp <- complement_pairs[suffix] # 配对后缀
  X <- substr(middle, 1, 1) # 中间部分的第一个碱基
  Y <- substr(middle, nchar(middle), nchar(middle)) # 中间部分的第二个碱基
  X_comp <- complement_pairs[X] # 配对中间的第一个碱基
  Y_comp <- complement_pairs[Y] # 配对中间的第二个碱基

  # 检查是否所有配对都存在
  if (!is.na(prefix_comp) && !is.na(suffix_comp) &&
      !is.na(X_comp) && !is.na(Y_comp)) {
    # 生成配对的 MutationType
    pair_mutation <- paste0(prefix_comp, "[", X_comp, ">", Y_comp, "]", suffix_comp)
    return(pair_mutation)
  } else {
    return(NA) # 如果任意配对失败，返回 NA
  }
}

# 对所有 MutationType 生成配对
df_sigProfile$Pair <- sapply(df_sigProfile$MutationType, get_pair)

# 创建统一的 pair 名字（两两配对无顺序差别）
get_unique_pair <- function(mutation, pair) {
  # 确保 MutationType 和 Pair 的顺序一致，取字典序最小的作为 key
  sorted_pair <- sort(c(mutation, pair))
  return(paste(sorted_pair, collapse = " + "))
}

# 计算唯一的配对
UniquePair <- mapply(get_unique_pair, df_sigProfile$MutationType, df_sigProfile$Pair)

# 为每个唯一配对分配一个 sig_pair 名字
unique_pairs <- unique(UniquePair) # 找到所有唯一配对
sig_pair_names <- paste0("pair", seq_along(unique_pairs)) # 为每个唯一配对命名
pair_map <- setNames(sig_pair_names, unique_pairs) # 创建配对映射

df_sigProfile$sig_pair <- pair_map[UniquePair]

# 使用 sapply 计算 binomial_test 相关的结果返回两列值
results <- t(sapply(seq_len(nrow(df_sigProfile)), function(i) {
  # print(i)
  # 当前行的 Count 和配对行的 Count
  a <- df_sigProfile$Count[i]
  pair_idx <- which(df_sigProfile$MutationType == df_sigProfile$Pair[i])
  
  # 如果没有配对数据，返回 NA 和 "no_pair"
  if (length(pair_idx) == 0) {
    return(c(greater_sig = "no_pair", p_value = NA))
  }
  
  b <- df_sigProfile$Count[pair_idx]
  
  # 确定 greater_sig 和较大较小值
  if (a > b) {
    greater_sig <- df_sigProfile$MutationType[i]
    larger <- a
    smaller <- b
  } else if (a < b) {
    greater_sig <- df_sigProfile$Pair[i]
    larger <- b
    smaller <- a
  } else {
    greater_sig <- "equal"
    larger <- a
    smaller <- b
  }
  
  # 计算 p 值
  if (a == 0 & b == 0) {
    return(c(greater_sig = greater_sig, p_value = 1))
  } else {
    test <- binom.test(smaller, larger, p = 0.5)
    return(c(greater_sig = greater_sig, p_value = test$p.value))
  }
}))

# 将结果保存到 df_sigProfile
df_sigProfile$greater_sig <- results[, "greater_sig"]
df_sigProfile$sig_pvalue <- as.numeric(results[, "p_value"])
# write.table(df_sigProfile, str_c(outputpath, "/muts_SigProfile.txt"), row.names=FALSE, col.names=TRUE,quote=FALSE, sep="\t")


##### Add information into site features dataframe
dim(df_candidate)
# [1] 176 182
# 1. 重命名 df_sigProfile 中的 MutationType 列为 RNAMutationType 以便匹配
colnames(df_sigProfile)[colnames(df_sigProfile) == "MutationType"] <- "RNAMutationType"
dim(df_sigProfile)
# [1] 192   7
# 2. 将 sig_pvalue 列加入 df_features，匹配 RNAMutationType 和 df_sigProfile 的 RNAMutationType
df_features <- merge(df_candidate, 
                     df_sigProfile[, c("RNAMutationType", "sig_pvalue")], 
                     by = "RNAMutationType", 
                     all.x = TRUE)
dim(df_features)
# [1] 176 183
# 3. 添加 signature_filter 列，根据 RNAMutationType 是否与 greater_sig 相同判断
df_removed <- df_sigProfile[
  df_sigProfile$RNAMutationType == df_sigProfile$greater_sig & df_sigProfile$sig_pvalue < 0.01, 
]
dim(df_removed)
# [1] 4 6
print(str_c("The number of mutation types would be removed is: ", as.character(nrow(df_removed))))
removed_muttype_by_sig <- df_removed$RNAMutationType
# [1] "G[T>G]A" "G[T>G]G" "T[G>A]A" "G[G>A]A"
df_features$signature_filter <- ifelse(
  df_features$RNAMutationType %in% removed_muttype_by_sig, 
  "fail", 
  "pass"
)
dim(df_features)
# [1] 176 184


##### geneate bed format columns
# 加载 tidyr 包
library(tidyr)
library(dplyr)

# 分割 identifier 列并保留原列
df_features <- df_features %>%
  mutate(original_identifier = identifier) %>%  # 复制一份 identifier 列
  separate(
    col = identifier, 
    into = c("chrom", "position", "ref", "alt"), 
    sep = "_", 
    remove = FALSE  # 保留原始 identifier 列
  ) %>%
  mutate(
    start = as.numeric(position) - 1,  # 第二个元素减 1 作为 start
    end = as.numeric(position)         # 第二个元素作为 end
  )

dim(df_features)


##### output and save
df_out <- df_features[, c('chrom', 'start', 'end', colnames(df_candidate), 'sig_pvalue', 'signature_filter')]; dim(df_out)
write.table(df_out, outputfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


##### check results
# # filter
# df_filter <- df_features[df_features$signature_filter=="pass", ]; dim(df_filter)
# # [1]  63 185

# df_true <- df_filter[df_filter$label=='mosaic',]; dim(df_true)
# # [1]  14 185

# dim(df_candidate[df_candidate$label=='mosaic',])
# # [1]  33 182


# # check
# df_check <- df_features[df_features$label == "mosaic" & df_features$signature_filter == "pass",]
# dim(df_check)
# # [1] 14  5
# head(df_check[, c("label", "RNAMutationType", "sig_pvalue", "signature_filter")])
#     label RNAMutationType   sig_pvalue signature_filter
# 5  mosaic         A[C>T]T 5.335212e-04             pass
# 9  mosaic         A[G>A]C 1.000000e+00             pass
# 15 mosaic         C[C>T]T 1.967493e-11             pass
# 20 mosaic         G[A>G]A 1.000000e+00             pass
# 40 mosaic         G[C>T]T 6.250000e-02             pass
# 41 mosaic         G[C>T]T 6.250000e-02             pass


##### Time #####
end_time <- Sys.time()
print(str_c("End time: ", end_time))
print(str_c("Program finished in ", as.character(round((end_time-start_time), 4)), " seconds"))

