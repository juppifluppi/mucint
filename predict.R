library("caret")

data=read.table("descriptors.csv",header=T)
print(colnames(data))
load("finmodel_9.rda")
#write.csv("results.csv",predict(final_model2,data))
print(predict(final_model2,data,type="prob"))
