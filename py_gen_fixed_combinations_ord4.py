NIND = 1000
NSNPS = 20000
NCOMBS = 96*(10**6)

def return_combs():
	c = 1
	gout = open("geno9_fixed_96m_order4", "w")
	for i in range(0,NSNPS):
		for j in range(i,NSNPS):
			for k in range(j,NSNPS):
				for l in range(k,NSNPS):
					gout.write(str(i) + " " + str(j) + " " + str(k) +  " " + str(l) + "\n")
					if c == NCOMBS:
						gout.close()
						return
					c += 1
	print("combs done\n")
	
return_combs()
