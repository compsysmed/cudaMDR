import random
from itertools import combinations
NIND = 1000
NSNPS = 20000
NCOMBS = 320000
with open("geno6", "w") as gout:
	s =  ""
	for i in range(0,NSNPS):
		for j in range(0,NIND):
			if random.random() < 0.6:
				s += "0 "
			elif random.random() < 0.7:
				s += "1 "
			else:
				s += "2 "
		s.rstrip()
		s += "\n"
	gout.write(s)
print("geno done\n")
	
with open("pheno6", "w") as pout:
	s = ""
	for i in range(0,NIND):
		if random.random() < 0.85:
			s += "0" + "\n"
		else:
			s += "1" + "\n"
	pout.write(s)
print("pheno done\n")

with open("combinations6", "w") as cout:
	for i in range(0,NCOMBS):
		e0 = random.randint(0,NSNPS)
		e1 = random.randint(0,NSNPS) #!
		e2 = random.randint(0,NSNPS)
		cout.write(str(e0) + " " + str(e1) + " " + str(e2) + "\n")
print("combs done\n")
