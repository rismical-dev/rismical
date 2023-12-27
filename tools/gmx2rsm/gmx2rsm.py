#! /usr/bin/python3
#
#  Usage: gmx2rsm.py [gro] [top] [ffnonbon]
#
#  [gro] gro file contains coordinate and atom name
#  [top] topology file contains charge and atom type
#  [ffnonbon] Non-bonding force field file contains LJ params
#

import sys
import pandas as pd

args=sys.argv

#########################
# Read gro
f=open(args[1],'r')
grodata=f.readlines()
f.close()

# Define dataframe
cols=["charge","sigma","epsilon","X","Y","Z","type"]
df=pd.DataFrame(columns=cols)

# Number of atoms
natom=int(grodata[1])

#
for i in grodata[2:natom+2]:
    data1=i.split()
    df=df.append({"X":data1[3],"Y":data1[4],"Z":data1[5]},ignore_index=True)


#########################
# Read top
f=open(args[2],'r')
topdata=f.readlines()
f.close()

flag=False
count=0
for i in topdata:
    if i=="\n":
        continue
    data2=i.split()
    if data2[0].startswith(";"):
        continue
    if data2[0]=="[" and data2[1]=="atoms" and data2[2]=="]":
        flag=True
        continue
    if flag==True:
        df.at[int(data2[0])-1,"type"]=str(data2[1])
        df.at[int(data2[0])-1,"charge"]=float(data2[6])

        count+=1
        if count==natom:
            break

#########################
# Read ffnonbon file
f=open(args[3],'r')
ljdata=f.readlines()
f.close()

flag=False
count=0
ljdict_s={}
ljdict_e={}
for i in ljdata:
    if i=="\n":
        continue
    data3=i.split()

    # if data3[0]==";":
    if data3[0].startswith(";"):
        continue
    if data3[0]=="[" and data3[1]=="atomtypes" and data3[2]=="]":
        flag=True
        continue
    if flag==True:
        ljdict_s[data3[0]]=float(data3[5])
        ljdict_e[data3[0]]=float(data3[6])
        
#########################
# Set sigma and epsilon
for i in range(0,natom):
    atomtype=df.at[int(i),"type"]
    df.at[int(i),"sigma"]=ljdict_s[atomtype]
    df.at[int(i),"epsilon"]=ljdict_e[atomtype]

#########################
# Print parameters
print(natom)
for i in range(0,natom):
    print(f'{df.at[int(i),"type"]}  {float(df.at[int(i),"sigma"])*10:.5f}  {float(df.at[int(i),"epsilon"])*1000:.5f} {df.at[int(i),"charge"]:.5f}  {float(df.at[int(i),"X"])*10:.5f}  {float(df.at[int(i),"Y"])*10:.5f}  {float(df.at[int(i),"Z"])*10:.5f}')

