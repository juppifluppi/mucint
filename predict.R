library("caret")

data=read.csv("descriptors.csv")
load("finmodel_9.rda")
print(predict(final_model2,data))
