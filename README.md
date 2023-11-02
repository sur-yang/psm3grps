**Extending the Classical Propensity Score Matching Algorithm for Three-Group Datasets**

## Introduction

The psm3grps package innovatively extends the classical propensity score matching (PSM) algorithm to datasets featuring three distinct groups.

The key features include:

Goal: Implements 1-1-1 PSM among three treatment groups.
Scores: Uses Generalized Propensity Score (GPS) via multinomial logistic regression, producing a probability score vector for each individual's group belonging.
Data Processing: Applies the rectangular common support region method to refine the sample, followed by GPS recalculation.
Matching: Adopts Lopez's vector matching method:
a. Select a reference group, perform k-means clustering with the other groups.
b. Conduct matching with replacement within clusters, using the reference group’s score.
c. Iterate steps with different reference groups.
d. Include individuals matched in all datasets for final analysis.

## Basic usage

```R
rm(list=ls())
data("exampleData",package="psm3grps")

# check unbalanced dataset
compareGroups::descrTable(treat~.,method = 4,data = exampleData)

# run the code
df.psm<-psm3grps(dataset)

# check the result
compareGroups::descrTable(treat~.,method = 4,data = df.psm)
```

## Credits

This package was  developed by [Xiaobao Yang](https://github.com/sur-yang) in collaboration with [Ruizhi Chai](https://github.com/sur-berry).
