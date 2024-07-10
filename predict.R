library("caret")

data=read.table("descriptors.csv",header=T)
print(colnames(data))
load("finmodel_9.rda")
write.table(file="results.csv",as.data.frame(predict(final_model2,data,type="prob")$X1),row.names=F)

load("finmodel_rfmix.rda")
write.table(file="results2.csv",as.data.frame(predict(final_model,data,type="prob")$X1),row.names=F)
print(as.data.frame(predict(final_model,data,type="prob")$X1))
