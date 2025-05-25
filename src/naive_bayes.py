from sklearn.naive_bayes import MultinomialNB
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import LabelEncoder
import numpy as np

import pandas as pd
import os

ref = pd.read_csv(os.path.join("data", "database", "12S_Kankakee_db_1117.tsv"), sep="\t")

# data from db for training
ref_seqs = list(ref["Sequence"])
ref_taxa = list(ref["Taxon"])

def kmer_tokenizer(seq, k=7): # tokenize sequences into 7-mers
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

vectorizer = CountVectorizer(analyzer=lambda x: kmer_tokenizer(x, k=7))

le = LabelEncoder()
y_train = le.fit_transform(ref_taxa)

X_train = vectorizer.fit_transform(ref_seqs)
clf = MultinomialNB()
clf.fit(X_train, y_train)

asv_seqs = pd.read_csv(os.file.path())

X_test = vectorizer.transform(asv_seqs)
y_pred = clf.predict(X_test)
probs = clf.predict_proba(X_test)
