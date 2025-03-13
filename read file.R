# 设置工作目录

setwd("/home/gary/下载/code/wam/Training Set")

extract_numbers_from_second_line <- function(filepath) {
  lines <- readLines(filepath)
  if (length(lines) >= 2) {
    second_line <- lines[2]
    matches <- regmatches(second_line, gregexpr("[0-9]+", second_line))[[1]]
    numbers <- as.numeric(matches)
    return(numbers)
  } else {
    return(numeric(0))  # 文件少于两行
  }
  
}
calculate_similarity <- function(predicted_data, actual_data) {
  # """
  # 计算两个 int 数组中相同数字的比例。
  #
  # Args:
  #   predicted_data (numeric): 预测的 int 数组。
  #   actual_data (numeric): 实际的 int 数组。
  #
  # Returns:
  #   numeric: 相同数字的比例。
  # """
  #
  # 找到两个数组中相同的数字
  common_numbers <- intersect(predicted_data, actual_data)
  
  # 计算相同数字的个数
  common_count <- length(common_numbers)
  
  # 计算比例（除以 actual_data 的长度）
  similarity <- common_count / length(actual_data)
  
  return(similarity)
}


# 列出所有txt文件
files <- list.files()

# 存储结果的列表
results <- list()

##########################3
# 遍历每个文件
results <- list() # 初始化结果列表

for (file in files) {
  # 读取文件内容
  content <- readLines(file)
  
  # 提取第二行（索引为2，因为R从1开始计数）
  cds_line <- content[2]
  
  # 使用正则表达式提取外显子位置
  exons_str <- gsub("[^0-9]", " ", cds_line)  # 去除非数字字符，用空格分隔
  exons_str <- trimws(exons_str) # 去除首尾空格
  exons_str <- strsplit(exons_str, "\\s+")[[1]] # 将字符串按空格分割成数字字符串向量
  exons <- as.numeric(exons_str) # 将数字字符串向量转换为数字向量
  
  # 按照起始和结束分组
  if (length(exons) > 0) {
    # 检查是否提取到外显子位置
    exon_positions <- matrix(exons, ncol = 2, byrow = TRUE) # 将数字向量转换为矩阵，每两列为一组起始和结束位置
    results[[file]] <- exon_positions # 存储结果
  } else {
    results[[file]] <- NULL # 如果没有提取到外显子位置，则存储NULL
  }
  # 提取核酸序列（从第三行开始）
  sequences <- content[3:length(content)]
  sequences <- gsub("\n", "", sequences)  # 去除换行符
  
  exon_sequences <- apply(exon_positions, 1, function(row) {
    start_start <- row[1] - 3
    start_end <- row[1] + 5
    end_start <- row[2] - 3
    end_end <- row[2] + 5
    
    # 将 sequences 连接成一个字符串
    genome_seq <- paste(sequences, collapse = "")
    
    # 检查起始和结束位置是否在基因组范围内
    if (start_start < 1) {
      start_start <- 1
    }
    if (start_end > nchar(genome_seq)) {
      start_end <- nchar(genome_seq)
    }
    if (end_start < 1) {
      end_start <- 1
    }
    if (end_end > nchar(genome_seq)) {
      end_end <- nchar(genome_seq)
    }
    
    start_seq <- substr(genome_seq, start_start, start_end)
    end_seq <- substr(genome_seq, end_start, end_end)
    
    return(c(start_seq, end_seq))
  })
  
  
  # 存储外显子信息和核酸序列
  results[[file]] <- list(exons = exons, res_sequences = exon_sequences)
}

############


# 查看结果（以第一个文件为例）
print(results[[1]])





















# 提取所有外显子序列
all_start_sequences <- unlist(lapply(results, function(x)
  x$res_sequences[1, ]))
all_end_sequences <- unlist(lapply(results, function(x)
  x$res_sequences[2, ]))

# 创建 WAM 模型 (分别针对起始和结束位置的序列)
library(Biostrings)
wam_start <- consensusMatrix(DNAStringSet(all_start_sequences))
wam_end <- consensusMatrix(DNAStringSet(all_end_sequences))

# 将频率矩阵转换为权重矩阵（可选）
wam_start_weights <- log2(wam_start / apply(wam_start, 2, sum))
wam_end_weights <- log2(wam_end / apply(wam_end, 2, sum))

# 定义打分函数
score_sequence <- function(sequence, wam) {
  sequence <- DNAString(sequence)
  score <- 0
  for (i in 1:length(sequence)) {
    base <- as.character(sequence[i])
    if (base %in% rownames(wam)) {
      score <- score + wam[base, i]
    }
  }
  return(score)
}

# 假设 other_files 是包含其他文件名的字符向量

setwd("/home/gary/下载/code/wam/Testing Set")


# 列出所有txt文件
other_files <- list.files()




# 初始化预测结果列表
all_predicted_start_sites <- c()
all_predicted_end_sites <- c()
all_predicted_scores_start <- c()
all_predicted_scores_end <- c()


# 初始化预测结果列表和得分列表
all_potential_start_sites <- c()
all_potential_end_sites <- c()
all_start_scores <- c()
all_end_scores <- c()

for (file in other_files) {
  # 读取文件内容
  other_content <- readLines(file)
  other_sequences <- paste(other_content[3:length(other_content)], collapse = "")
  numbers <- extract_numbers_from_second_line(file) - 3
  
  if (is.null(numbers)) {
    print("文件不存在或文件内容不正确。")
  } else if (length(numbers) > 0) {
    print(paste("提取到的数字:", paste(numbers, collapse = ", ")))
  } else {
    print("未找到数字。")
  }
  # 设置窗口大小（与 WAM 模型长度相同）
  window_size <- ncol(wam_start) # 假设起始和结束位置的序列长度相同
  
  # 设置得分阈值
  start_threshold <- -11.4 # 你需要根据实际情况调整阈值
  end_threshold <- -9.3
  
  # 滑动窗口打分
  potential_start_sites <- c()
  potential_end_sites <- c()
  start_scores <- c()
  end_scores <- c()
  for (i in 1:(nchar(other_sequences) - window_size + 1)) {
    window <- substr(other_sequences, i, i + window_size - 1)
    start_score <- score_sequence(window, wam_start_weights)
    end_score <- score_sequence(window, wam_end_weights)
    
    # 如果得分高于阈值，则将窗口位置添加到潜在位点列表
    if (start_score > start_threshold) {
      potential_start_sites <- c(potential_start_sites, i)
      start_scores <- c(start_scores, start_score) # 存储得分
      #print(start_score) # 注释掉打印语句
    }
    if (end_score > end_threshold) {
      potential_end_sites <- c(potential_end_sites, i)
      end_scores <- c(end_scores, end_score) # 存储得分
      #print(end_score) # 注释掉打印语句
    }
  }
  
  # 将当前文件的预测结果添加到总列表
  # all_potential_start_sites <- c(all_potential_start_sites, potential_start_sites)
  # all_potential_end_sites <- c(all_potential_end_sites, potential_end_sites)
  # all_start_scores <- c(all_start_scores, start_scores)
  # all_end_scores <- c(all_end_scores, end_scores)
  
  
  mergedArray <- c(potential_start_sites, potential_end_sites)
  
  # 打印潜在位点（可选，可以注释掉）
  if (length(potential_start_sites) > 0 ||
      length(potential_end_sites) > 0) {
    cat("File:", file, "\n")
    if (length(potential_start_sites) > 0) {
      cat("  Potential start sites:", potential_start_sites, "\n")
    }
    if (length(potential_end_sites) > 0) {
      cat("  Potential end sites:", potential_end_sites, "\n")
    }
  } else {
    cat("File:", file, "No potential sites found.\n")
  }
  similarity <- calculate_similarity(mergedArray, numbers)
  print(paste("Similarity:", similarity))
  
  
  
  
  
  
  
  
  
  
  
}

# 绘制ROC曲线部分，需提供实际位点信息actual_sites

# 准备数据，用于绘制 ROC 曲线
# labels_start <- ifelse(1:length(other_sequences) %in% actual_sites, 1, 0)
# labels_end <- ifelse(1:length(other_sequences) %in% actual_sites, 1, 0)
#
# scores_start <- rep(0, length(other_sequences))
# scores_start[all_potential_start_sites] <- all_start_scores
#
# scores_end <- rep(0, length(other_sequences))
# scores_end[all_potential_end_sites] <- all_end_scores
#
# data_start <- data.frame(labels = labels_start, scores = scores_start)
# data_end <- data.frame(labels = labels_end, scores = scores_end)

# 绘制 ROC 曲线
# library(pROC)
#
# roc_start <- roc(data_start$labels, data_start$scores)
# plot(roc_start, main = "ROC Curve (Start Sites)")
# auc(roc_start)
#
# roc_end <- roc(data_end$labels, data_end$scores)
# plot(roc_end, main = "ROC Curve (End Sites)")
# auc(roc_end)













# # for (file in other_files) {
# #   # 读取文件内容
# #   other_content <- readLines(file)
# #   other_sequences <- other_content[3:length(other_content)]
# #   other_sequences <- paste(other_sequences, collapse = "")
# #
# #   # 设置窗口大小（与 WAM 模型长度相同）
# #   window_size <- ncol(wam_start) # 假设起始和结束位置的序列长度相同
# #
# #   # 设置得分阈值
# #   start_threshold <- -11.4 # 你需要根据实际情况调整阈值
# #   end_threshold <- -9.2
# #
# #   # 滑动窗口打分
# #   potential_start_sites <- c()
# #   potential_end_sites <- c()
# #   for (i in 1:(nchar(other_sequences) - window_size + 1)) {
# #     window <- substr(other_sequences, i, i + window_size - 1)
# #     start_score <- score_sequence(window, wam_start_weights)
# #     end_score <- score_sequence(window, wam_end_weights)
# #
# #     # 如果得分高于阈值，则将窗口位置添加到潜在位点列表
# #     if (start_score > start_threshold) {
# #       potential_start_sites <- c(potential_start_sites, i)
# #       print(start_score)
# #     }
# #     if (end_score > end_threshold) {
# #       potential_end_sites <- c(potential_end_sites, i)
# #       print(end_score)
# #     }
# #   }
# #
# #   # 打印潜在位点
# #   if (length(potential_start_sites) > 0 || length(potential_end_sites) > 0) {
# #     cat("File:", file, "\n")
# #     if (length(potential_start_sites) > 0) {
# #       cat("  Potential start sites:", potential_start_sites, "\n")
# #     }
# #     if (length(potential_end_sites) > 0) {
# #       cat("  Potential end sites:", potential_end_sites, "\n")
# #     }
# #   } else {
# #     cat("File:", file, "No potential sites found.\n")
# #   }
# #  
# #   
# #   #for循环结束
# # }
# 
# 
# 
























