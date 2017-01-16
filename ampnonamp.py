# -*- coding: utf-8 -*-
from __future__ import division
"""
Created on Tue Dec 13 15:14:18 2016

@author: Ashiqul
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 16:49:34 2016

@author: Ashiqul
"""

import pandas as pd
import numpy as np
from sklearn import datasets 
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Perceptron
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.feature_extraction.text import CountVectorizer
from Bio import SeqIO
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier 
from sklearn.metrics import matthews_corrcoef
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import ShuffleSplit
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import average_precision_score
from sklearn.model_selection import cross_val_predict
corpus1=[]#make a blank list of corpus for the training set
tag=[]#make the list of outcomes
true = open("C:/Users/Ashiqul/Desktop/amp.fasta")
false = open("C:/Users/Ashiqul/Desktop/non_amp.fasta")
for line in SeqIO.parse(true, "fasta"):
    line = line.seq.tostring().lower().replace("x","")
    line = line.replace('-', "")
   # print line
    tag.append("1")
    fullstring = extract_motifs_pos(line, 1, 10, 1, 20, 2, 3)
    #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
    corpus1.append(fullstring) #apperd string from each protein to corpus
true.close()
for line in SeqIO.parse(false, "fasta"):
    line = line.seq.tostring().lower().replace("x","")
    line = line.replace('-', "")
    #print line
    tag.append("0")
    fullstring = extract_motifs_pos(line, 1, 10, 1, 20, 2, 3)
    #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
    corpus1.append(fullstring) #apperd string from each protein to corpus 
false.close()

corpus = np.array(corpus1) #convert corpus into numpy array
tag = np.array(tag)  # convert tag into numpy array   
#print corpus # print for debugging 
#print tag # print for debugging
for i in range(500, 3000, 250):
    count = CountVectorizer(max_features=15000000, vocabulary = None, max_df=0.3, min_df = 3, stop_words=[1,2])#giving the CountVectorizer function a short name
    #get the vocabulary of train set to use for the test set
    
    
    
    bag = count.fit_transform(corpus) #transform the corpus(bag of words into sparse martrix)
    #print (count.vocabulary_) #count the occurence of number of words
    ##get the vocabulary of train set to use for the test set. Next time put the "voc" in
    #the vocabulary parameter of count
    voc = count.get_feature_names() 
    #print len(voc)
    bag= bag.toarray() 
    np.place(bag, bag>0, [1])#convert the sparsematrix into an array
    #np.place(bag, bag>0, [1])
    #print bag
    
    forest = RandomForestClassifier(n_estimators = 500,
                                    random_state = 1,
                                    n_jobs =1)
    forest.fit(bag[:, 0:-1], tag)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
    # Print the feature ranking
    #print("Feature ranking:")
    important = list()
    for f in range(0,i):
        important.append(indices[f])
        
    bag=bag[:,important]
    voc = np.array(voc)[important]
    pca = PCA(n_components=30, svd_solver='full')
    pca.fit(bag)
    plt.plot(pca.explained_variance_, linewidth=2) 
    bag = pca.transform(bag)  
    bag=pd.DataFrame(bag)
    bag['tag']=tag
    
    x = bag.iloc[:, 0:-1]
    y = bag.iloc[:,-1]
    #lda = LinearDiscriminantAnalysis(n_components=200)
    #lda.fit(x,y)
    #X = lda.transform(x)
    #parameterize the Logistic Regression algorithm
    cv = ShuffleSplit(n_splits=10, test_size=0.2, random_state = 1)
    clf=  SVC(kernel = 'rbf', C=1.0, gamma = 0.05, random_state=0)
    print "When number of fearture is: "+ str(i) + "accuracy is: "+ str(np.mean(cross_val_score(clf, x, y, cv =cv))) 


    
"""
corpus1=[]#make a blank list of corpus for the training set
tag=[]#make the list of outcomes
true = open("C:/Users/Ashiqul/Desktop/testamp.fasta")
false = open("C:/Users/Ashiqul/Desktop/testnonamp.fasta")
for line in SeqIO.parse(true, "fasta"):
    line = line.seq.tostring().lower().replace("x","")
    line = line.replace('-', "")
   # print line
    tag.append("1")
    fullstring = extract_motifs_pos(line, 1, 10, 1, 20, 2, 3)
    #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
    corpus1.append(fullstring) #apperd string from each protein to corpus
true.close()
for line in SeqIO.parse(false, "fasta"):
    line = line.seq.tostring().lower().replace("x","")
    line = line.replace('-', "")
    #print line
    tag.append("0")
    fullstring = extract_motifs_pos(line, 1, 10, 1, 20, 2, 3)
    #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
    corpus1.append(fullstring) #apperd string from each protein to corpus 
false.close()

corpus = np.array(corpus1) #convert corpus into numpy array
tag = np.array(tag)  # convert tag into numpy array   
#print corpus # print for debugging 
#print tag # print for debugging

count = CountVectorizer(vocabulary = voc)
#get the vocabulary of train set to use for the test set



bag_test = count.fit_transform(corpus) #transform the corpus(bag of words into sparse martrix)
bag_test= bag_test.toarray()
np.place(bag_test, bag_test>0, [1])
bag_test = pca.transform(bag_test)
clf = SVC(kernel = 'rbf', C=1.0, gamma = 0.05, random_state=0)
    
x_train = x
y_train = y 
x_test = bag_test
y_test = tag
# fit the Logistic Regression model
clf.fit(x_train, y_train)
y_pred= clf.predict(x_test)


print 'Accuracy: ' + str(accuracy_score(y_test, y_pred))        
print matthews_corrcoef(y_test, y_pred)   

"""


