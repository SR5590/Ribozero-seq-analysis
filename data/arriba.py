import sys
from time import sleep
#aa="ccp_rna_module/data/8RPVZEsomatic.fusion.tsv"
aa=sys.argv[1]
out=sys.argv[3]
counter=1

outfile=open(sys.argv[2],"w")
with open(aa) as s:
    line1 = s.readline()
    outfile.write(line1)
    for line in s:
        print(line)
        line_original=line
        line=line.split("\t")
        #print(line.split("\t")[16])
        if line[16]=="high" or line[16]=="medium":
            sleep(1)
            print(int(line[11]))
            r1=int(line[11])
            r2=int(line[12])
            Dis_mates= int(line[13])
            cov1=int(line[14])
            cov2=int(line[15])

            vv=[x for x in [cov2,cov1,r1,r2,Dis_mates] if x != 0]
            sleep(1)

            if len(vv)>2:
                sleep(1)
                c=(float(r1)+float(r2)+float(Dis_mates))/(float(cov1)+float(cov2))
                sleep(1)
                print(float(c))

                if c >= 0.01:
                    with open(out+"Fusion"+str(counter), "w")as t:
                        t.write(line1)
                        t.write(line_original)
                        outfile.write(line_original)
                        counter=counter+1
                        sleep(1)

outfile.close()
with open(out+".done", "w") as t:
    t.write("found " + str(counter) + "fusions with high or medium confidence")
