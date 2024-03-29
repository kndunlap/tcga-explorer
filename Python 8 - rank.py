#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:34:57 2024

@author: u1314571
"""

all = pd.read_csv('TCGA_all_tidy_with_metadata.csv')

all = all.drop(columns = ['Unnamed: 0'])

all = all[['Type', 'ASS1']]

all.groupby('Type')['ASS1'].median()




def rank (gene):
    all_data = (
       all
       .drop(columns=['Unnamed: 0'])
       .loc[:, ['Type', gene]]
       .groupby('Type')[gene]
       .median()
       .sort_values(ascending = False)
    )
    return all_data

rank('OTC')
