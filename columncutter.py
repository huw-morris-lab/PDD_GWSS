import csv

def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n): 
        yield l[i:i + n]
        
size = 10000
start =  6
        
f = open("${COHORT}_SNPsDF_recodeA.raw", "r")

line = f.readline()
bits = line.split()
FID = bits[0]
IID = bits[1]
PAT = bits[2]
MAT = bits[3]
SEX = bits[4]
PHENOTYPE  = bits[5]


partchunks = list(divide_chunks(bits[start:], size))
mylist = range(0,len(partchunks))
myliststr = map(str, mylist) 
pre_res = ["${COHORT}_file" + sub for sub in myliststr]
suf_res = [sub + ".raw" for sub in pre_res]

i=0
for filename in suf_res:
    with open(filename, 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow([FID,IID,PAT,MAT,SEX,PHENOTYPE] + partchunks[i])
        csvfile.close()
        i=i+1

while f:
    line  = f.readline()
    if line == "":
        break
    bits = line.split()
    FID = bits[0]
    IID = bits[1]
    PAT = bits[2]
    MAT = bits[3]
    SEX = bits[4]
    PHENOTYPE  = bits[5]
    partchunks = list(divide_chunks(bits[start:], size))
    mylist = range(0,len(partchunks))
    myliststr = map(str, mylist) 
    pre_res = ["${COHORT}_file" + sub for sub in myliststr]
    suf_res = [sub + ".raw" for sub in pre_res]

    i=0
    for filename in suf_res: 
        with open(filename, 'a') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow([FID,IID,PAT,MAT,SEX,PHENOTYPE] + partchunks[i])
            csvfile.close()
            i=i+1
