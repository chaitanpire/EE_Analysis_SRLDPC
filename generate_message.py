import random
#generate a random number between 0 and 1
# p = random.random()
p = 0.5
#now for 5888 times with probability p, print "1" and with pobability 1-p, print "0"
for i in range(5888):
    if random.random() < p:
        print("1")
    else:
        print("0")

