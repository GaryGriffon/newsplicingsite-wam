TP <- sum_right_point
FP <- sum_wrong_point
FN <- sum_point - sum_right_point

Recall <- TP / (TP + FN)
Precision <- TP / (TP + FP)
F1score <- 2*(Recall*Precision)/(Precision+Recall)

cat("  Sn/Recall:", Recall, "\n")
cat("  Precision:", Precision, "\n")
cat("  F1-score:", F1score, "\n")