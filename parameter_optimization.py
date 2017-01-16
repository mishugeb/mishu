# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 22:19:22 2017

@author: Ashiqul
"""

def parameters(f1, f2, x1, x2, x3, x4, x5, x6, x7):

    corpus1=[]#make a blank list of corpus for the training set
    tag=[]#make the list of outcomes
    true = open(f1)
    false = open(f2)
    for line in SeqIO.parse(true, "fasta"):
        line = line.seq.tostring().lower().replace("x","")
        line = line.replace('-', "")
        line = line.upper()
       # print line
        tag.append("1")
        fullstring = extract_motifs_pos(line,x1, x2, x3, x4, x5, x6, x7)
        #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
        corpus1.append(fullstring) #apperd string from each protein to corpus
    true.close()
    for line in SeqIO.parse(false, "fasta"):
        line = line.seq.tostring().lower().replace("x","")
        line = line.replace('-', "")
        line = line.upper()
        #print line
        tag.append("0")
        fullstring = extract_motifs_pos(line, x1, x2, x3, x4, x5, x6, x7)
        #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
        corpus1.append(fullstring) #apperd string from each protein to corpus 
    false.close()
    """
    fh = open("C:/Users/Ashiqul/Desktop/Motif.txt")
    voc = []
    for i in fh:
        print i.split(" ")
        line = i.split(" ")[9]
        line= line.rstrip("\n")
        voc.append(line)
    fh.close()
    """
    try:    
        #voc=['RDGYP', 'GKCMN', 'GYCYW', 'KCYCY', 'WCVWD', 'CCEHL', 'HKWCK', 'CADWA', 'WKWCV', 'LPDNV',
             #'PMH', 'WYP', 'WWKCG', 'PCFTT', 'NYPLD', 'PMKCI', 'PWC', 'NDYCL', 'LWCRY', 'MWTCD']
        corpus = np.array(corpus1) #convert corpus into numpy array
        tag = np.array(tag)  # convert tag into numpy array   
        #print corpus # print for debugging 
        #print tag # print for debugging
        
        count = CountVectorizer(max_features=15000000, vocabulary = None, max_df=0.3, min_df = 3, stop_words=[1,2])#giving the CountVectorizer function a short name
        #get the vocabulary of train set to use for the test set
        
            
            
        
        
        bag = count.fit_transform(corpus) #transform the corpus(bag of words into sparse martrix)
        #print (count.vocabulary_) #count the occurence of number of words
        ##get the vocabulary of train set to use for the test set. Next time put the "voc" in
        #the vocabulary parameter of count
        voc = count.get_feature_names() 
        #print len(voc)
        bag= bag.toarray() #convert the sparsematrix into an array
        np.place(bag, bag>0, [1])
        #print bag
    except:
        return 0.0
    
    forest = RandomForestClassifier(n_estimators = 1000,
                                    random_state = 1,
                                    n_jobs =1)
    forest.fit(bag[:, 0:-1], tag)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
    
    # Print the feature ranking
    print("Feature ranking:")
    important = list()
    for f in range(0,500):
        important.append(indices[f])
        
    bag=bag[:,important]
    voc = np.array(voc)[important]
    pca = PCA(n_components=10, svd_solver='full')
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
    cv = ShuffleSplit(n_splits=10, test_size=0.2, random_state =1)
    clf= LogisticRegression(C=1.0,random_state=0)
    return np.mean(cross_val_score(clf, x, y, cv =cv)) 
    
def tune(f1, f2):
    
    x2=2;x4=1;x5=1; x6=0; x7=1;
    for i in range(0,5):
        print "parameter index: ", i
        x= 0
        p = 0.0
        
        
        p_best= -1
        
        while True:
            if p_best < p:
                count = 0
                p_best = p
            else:
                count+=1
            print count
            
            if i == 0:
                x4=  x
            elif i == 1:
                x6= x
            elif i== 2:
                x7= x
            elif i==3:
                x2=x
            elif i==4:
                x5= x
            p = parameters(f1, f2, 0, x2, 0, x4, x5, x6, x7)
            print x   
            x += 1 
            print "new", p
            if count == 10:
                break
    
    
        if i == 0:
            x4=  x-11
            print "x4: ", x4
        elif i == 1:
            x6= x-11
            print "x6: ", x6
        elif i== 2:
            x7= x-11
            print "x7: ", x7
        elif i==3:
            x2=x-11
            print "x2: ", x2
        elif i==4:
            x5= x-11
            print "x5: ", x5
            
    

    return (0, x2, 0, x4, x5, x6, x7)

print tune("C:/Users/Ashiqul/Desktop/tumor.fasta", "C:/Users/Ashiqul/Desktop/non_tumor.fasta")


def tune1(f1, f2):
    
    x2=1;x4=10;x5=1; x6=1; x7=1;
    for i in range(0,2):
        print "parameter index: ", i
        x= 1
        p = 0.0
        
        
        p_best= -1
        
        while True:
            if p_best < p:
                count = 0
                p_best = p
            else:
                count+=1
            print count
            
            if i == 0:
                x2=  x
            elif i == 1:
                x5= x
            
            p = parameters(f1, f2, 0, x2, 0, x4, x5, x6, x7)
            print x   
            x += 1 
            print "new", p
            if count == 10:
                break
    
    
        if i == 0:
            x2=  x-12
            print "x2: ", x2
        elif i == 1:
            x5= x-12
            print "x5: ", x5
        
            
    

    return (x2, x5)


def tune2(f1, f2):
    
    x2=10;x4=1;x5=1; x6=1; x7=1;
    for i in range(0,3):
        print "parameter index: ", i
        x= 1
        p = 0.0
        
        
        p_best= -1
        
        while True:
            if p_best < p:
                count = 0
                p_best = p
            else:
                count+=1
            print count
            
            if i == 0:
                x4=  x
            elif i == 1:
                x6= x
            elif i== 2:
                x7= x
            
            p = parameters(f1, f2, 0, x2, 0, x4, x5, x6, x7)
            print x   
            x += 1 
            print "new", p
            if count == 10:
                break
    
    
        if i == 0:
            x4=  x-12
            print "x4: ", x4
        elif i == 1:
            x6= x-12
            print "x6: ", x6
        elif i== 2:
            x7= x-12
            print "x7: ", x7
        
    

    return (x4, x6, x7)

def tune_new(f1, f2):
    x2, x5 = tune1(f1, f2)
    x4, x6, x7 = tune2(f1, f2)
    return 1, x2, 1, x4, x5, x6, x7
print tune_new("C:/Users/Ashiqul/Desktop/angio_test1.fasta", "C:/Users/Ashiqul/Desktop/nonangio_test1.fasta")

