#!/usr/bin/env python3

#start = 0
#end   = 99
#divisor=7
#print("Printing out numbers from",start,"to",end, " not divisible by",divisor)
for i in range(100):
    if i == 0:
        print(i)
    if i % 7 == 0:
        continue
    print(i)
