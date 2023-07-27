#CREATE INITIAL GRO FILE
import numpy as np
import os
import random
seq_to_resname = { \
	'A': 'ALA',
	'G': 'GLY',
	'S': 'SER',
	'I': 'ILE',
	'L': 'LEU',
	'V': 'VAL',
	'M': 'MET',
	'D': 'ASP',
	'Q': 'GLN',
	'N': 'ASN',
	'P': 'PRO',
	'T': 'THR',
	'K': 'LYS',
	'R': 'ARG',
	'E': 'GLU',
	'W': 'TRP',
	'F': 'PHE',
	'Y': 'TYR',
	'H': 'HIS',
	'Rsd': 'RSD'
}
seq_to_resno = { \
	'A': 3,
	'G': 3,
	'S': 6,
	'I': 4,
	'L': 4,
	'V': 4,
	'M': 4,
	'D': 5,
	'Q': 6,
	'N': 6,
	'P': 4,
	'T': 6,
	'K': 5,
	'R': 5,
	'E': 5,
	'W': 9,
	'F': 6,
	'Y': 8,
	'H': 7,
	'Rsd': 6
}
f = open("seq.txt", 'r')
A = f.readlines()
f.close()
string = A[0].strip()
Seq = [seq_to_resname[i] for i in string]
Seqno = [seq_to_resno[i] for i in string]
Dict = {'GLY':'BB BBp BBm','SER':'BB BBp BBm S1 S1p S1m', 'ASN':'BB BBp BBm S1 S1p S1m', 'LYS':'BB BBp BBm S1 S2', 'ALA':'BB BBp BBm', 'ILE':'BB BBp BBm S1', 'LEU':'BB BBp BBm S1', 'MET':'BB BBp BBm S1', 'PHE': 'BB BBp BBm S1 S2 S3', 'VAL': 'BB BBp BBm S1', 'GLU': 'BB BBp BBm S1 S2', 'ASP': 'BB BBp BBm S1', 'ARG': 'BB BBp BBm S1 S2', 'HIS': 'BB BBp BBm S1 S2 S2p S2m', 'TYR': 'BB BBp BBm S1 S2 S3 S3p S3m', 'GLN': 'BB BBp BBm S1 S1p S1m', 'THR': 'BB BBp BBm S1 S1p S1m', 'TRP': 'BB BBp BBm S1 S2 S2p S2m S3 S4', 'PRO': 'BB BBp BBm S1'}
f = open("temp.gro", "w+")
f.write( string +'\n')
f.write( str(sum(Seqno)) +'\n')
bb_dummy = 0.14
bb_bb = 0.385
bb_s1 = 0.25
s1_s2 = 0.28
s1_dummy = 0.14
s2_dummy = 0.14
basex = 0.001
Change_x = {'BB': 0.0, 'BBm': 0.0, 'BBp': 0.0, 'S1': 0.0,'S2': 0.0, 'S1p': 0.0, 'S1m': 0.0, 'S2p': 0.0, 'S2m': 0.0}
Change_y = {'BB': 0.0, 'BBm': +0.141, 'BBp': -0.142, 'S1': 0.0,'S2': 0.0, 'S1p': 0.141, 'S1m': -0.141, 'S2p': 0.141, 'S2m': -0.141}
Change_z = {'BB': 0.0, 'BBm': 0.0, 'BBp': 0.0, 'S1': 0.251,'S2': 0.501, 'S1p': 0.251, 'S1m': 0.251, 'S2p': 0.501, 'S2m': 0.501}
resid = 0
no = 0
b = 0.27
for i in Seq:
    toss = random.choice([-1,1])
    Change_y.update({n: toss * Change_y[n] for n in Change_y.keys()})
    print Change_y
    resid=resid+1
    basex = basex+0.385
    basey = 0.001
    basez = 0.001
    if i not in ['TYR', 'TRP', 'PHE']:
        for j in Dict[i].split():
            if (j == 'S1p' or j == 'S1m'):
                #print pos
                no=no+1
                pos_s1p = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                #print resid, i, j, no, pos_s1p
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos_s1p[0])) + '%8s' % (str(pos_s1p[1])) + '%8s' % (str(pos_s1p[2])) + '\n')
            elif (j == 'S2p' or j == 'S2m'):
                no=no+1
                pos_s1p = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                #print resid, i, j, no, pos_s1p
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos_s1p[0])) + '%8s' % (str(pos_s1p[1])) + '%8s' % (str(pos_s1p[2])) + '\n')
            else:
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
    if i == 'PHE':
        for j in Dict[i].split():
            if j == 'BB':
                pos = [basex, basey, basez]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'BBp':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'BBm':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S1':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
                s1 = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
            if j == 'S2':
                pos = [s1[0]-b/2.0, s1[1], s1[2]+(0.233)]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S3':
                pos = [s1[0]+b/2.0, s1[1], s1[2]+(0.233)]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
    if i == 'TYR':
        for j in Dict[i].split():
            if j == 'BB':
                pos = [basex, basey, basez]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'BBp':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'BBm':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S1':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
                s1 = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
            if j == 'S2':
                pos = [s1[0]-b/2.0, s1[1], s1[2]+(0.233)]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S3':
                pos = [s1[0]+b/2.0, s1[1], s1[2]+(0.233)]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S3p':
                newtoss = random.choice([-1,1])
                pos = [s1[0]+b/2.0, s1[1]+newtoss*0.141, s1[2]+(0.233)]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S3m':
                pos = [s1[0]+b/2.0, s1[1]-newtoss*0.141, s1[2]+(0.233)]
                no=no+1
                #print resid, i, j, no, pos
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
    if i == 'TRP':
        #print Dict[i].split()
        for j in Dict[i].split():
            if j == 'BB':
                pos = [basex, basey, basez]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'BBp':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'BBm':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S1':
                pos = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
                s1 = [basex+Change_x[j], basey+Change_y[j], basez+Change_z[j]]
            if j == 'S2':
                pos = [s1[0]-b/2.0, s1[1], s1[2]+0.233]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S2p':
                newtoss = random.choice([-1,1])
                pos = [s1[0]-b/2.0, s1[1]+newtoss*0.141, s1[2]+0.233]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S2m':
                pos = [s1[0]-b/2.0, s1[1]-newtoss*0.141, s1[2]+0.233]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S3':
                pos = [s1[0]+b/2.0, s1[1], s1[2]+0.233]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
            if j == 'S4':
                pos = [s1[0], s1[1], s1[2]+(0.381)]
                no=no+1
                f.write('%5s' % (str(resid)) + '%3s' % (str(i)) + '%7s' % (str(j)) + '%5s' % (str(no)) + '%8s' % (str(pos[0])) + '%8s' % (str(pos[1])) + '%8s' % (str(pos[2])) + '\n')
                  
f.write('\t'+str(basex+5)+'\t'+str(basex+5)+'\t'+str(basex+5))
f.close()

os.system("gmx genconf -f temp.gro -o protein.gro")
