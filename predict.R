library("caret")

data=read.csv("descriptors.csv")
load("finmodel_9.rda")
write.csv("results.csv",predict(final_model2,data))
print(predict(final_model2,data))
