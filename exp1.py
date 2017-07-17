#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 10:50:08 2017

@author: frejj
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 08:45:21 2017

@author: frejj
"""

# import modules & set up logging
#import sys

#sys.path.append('/home/mrim/frejj/WORK/gensim-1.0.1/')


import logging
#import numpy as np
import xml.etree.ElementTree as ET


f = open('../data/collection/porter_stop_string_content', 'w')


logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
#path = '/home/mrim/frejj/Terrier4.2/terrier-core-4.2/share/vaswani_npl/corpus/'
#path = '/home/mrim/frejj/pyword2vec/'
print('before tree')
tree = ET.parse('../data/collection/PORTER_STOP_DOC_CHIC2012.xml')
root = tree.getroot()
content = []

print('after tree')

print(root[0][1].text)

for i in range(len(root)):
    #if i%100000 == 0:
    if root[i][1].text is not None:
        #content.append(root[i][1].text.split())
        f.write(root[i][1].text.encode('utf8') + '\n')     
        
    if root[i][1].text is None:
        #content.append(root[i][1].text.split())
        f.write('\n')  
        
      
"""

        

print('loading file')


start = time.time()

content = pickle.load( open( "content", "rb" ) )


end = time.time()
print('Time load : ' , end - start)

print('Size content : ',len(content))

print(content[0])

voc = []

start = time.time()


for i in range(len(content)):
    for j in content[i]:
        if not j in voc:
            voc.append(j)
            
end = time.time()
print('Time calc voc : ' , end - start)



pickle.dump( voc , open( "voc", "wb" ) )

voc2 = pickle.load( open( "voc", "rb" ) )




cf = [0] * len(voc)


for i in range(len(content)):
    for j in content[i]:
       cf[ voc.index(j) ] += 1
            
             
model = gensim.models.Word2Vec(content, min_count=1)




print(model.similarity('increase', 'decrease'))
print(model.similarity('increase', 'dissociated'))

model.save('mymodel')
new_model = gensim.models.Word2Vec.load('mymodel')

"""