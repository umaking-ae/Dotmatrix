import streamlit as st

from Bio import SeqIO
#ファイルからの読み込み
record1 = next(SeqIO.parse("NC_045512_Wuhan.fa","fasta"))
record2 = next(SeqIO.parse("sars.fa","fasta"))
#配列の取り出し
seq1 = record1.seq
seq2 = record2.seq

import numpy as np
import matplotlib.pyplot as plt

win = 10 
len1 = len(seq1) - win + 1 #配列1の長さ
len2 = len(seq2) - win + 1 #配列2の長さ
width = 500 #実際に描く幅
height = 500 #実際に描く高さ
image = np.zeros((height,width)) #実際に描く幅・高さの行列

hash={}

for x in range(len1):
    subseq1 = seq1[x : x + win] 
    if subseq1 not in hash:
        hash[subseq1] = []
    hash[subseq1].append(x)

for y in range(len2):
    sub2 = seq2[y : y + win]
    py = int(y / len2 * height)  #yを画像位置pyに
    if sub2 in hash:
        for x in hash[sub2]:
            px = int(x / len1 * width)  #xを画像位置pxに
            image[py, px] = 1  #[py,px]にドット

plt.imshow(image, extent = (1, len1, len2, 1), cmap = "Grays") #行列表示

#plt.show()

st.pyplot(plt)