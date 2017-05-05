#!/bin/env python

outtext = []
for i in range(0,25):
    start=i*10
    end=i*10+9
    file = open('predictions_%s_%s.csv' % (start, end), 'r')
    text = file.readlines()
    # Replace first element with count of entry we are working on
    for idx, line in enumerate(text[1:]):
        # Replace old number (0 thru 9) with new one
        text[idx+1] = str(start+idx)+line[1:] 

    # Strip header if not first case
    if i!=0:
        text=text[1:]

    outtext+=text

    file.close()

file=open('predictions.csv', 'w')
file.writelines(outtext)
