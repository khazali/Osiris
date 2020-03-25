file2=open("Petsc_Functions.cpp", "rt")
file1=open("Do_Cycle2.cpp", "rt")

row1 = file1.readlines()
row2 = file2.readlines()

i=0
for line1 in row1:
    i=i+1
    flag=False
    for line2 in row2:
        if line1.strip()==line2.strip():
            flag=True
            break
    if (i < 12515) and (not flag) and (line1.find('ShapeFactor')==(-1)):
        print(i,' : ', line1)
            
    

file1.close()
file2.close()