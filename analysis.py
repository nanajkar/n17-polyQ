import numpy as np
import numpy as np
import MDAnalysis as mda
from MDAnalysis import *
import networkx as nx
from MDAnalysis.lib.distances import distance_array
from scipy.spatial import ConvexHull


def COM_dist(tpr,xtc):
    '''
    Function to get distance of PHE sidechain beads from agg center of mass
    '''

    u = Universe(tpr,xtc)
    phe = u.select_atoms('resname PHE and name S1 S2 S3')

    all_dat = []

    for ts in u.trajectory[0::1]:
        agg_com = u.select_atoms('name BB S1 S2 S3 S4').center_of_mass()
        phe = u.select_atoms('resname PHE and name S1 S2 S3').groupby('resids')
        dist = [np.abs(phe[i].center_of_mass() - agg_com) for i in phe]
        all_dat.append(np.mean(dist))

    return all_dat

def avg_GLN_edge(tpr,xtc):
    '''
    Calculates avg number edges per GLN residue. Excludes bonded interactions, but includes inter and intra-Q interaction
    Output in domain_int/phi
    '''
    u = Universe(tpr,xtc)
    gln = u.select_atoms('resname GLN and name BB S1')
    avg_edges = []
    for ts in u.trajectory:
        G = nx.Graph()
        D = np.where(distance_array(gln.positions,gln.positions,box=ts.dimensions)<7)
        aa_x =  D[0]//2
        aa_y =  D[1]//2
        for bead in range(len(aa_x)):
            if np.abs(aa_x[bead]-aa_y[bead])>1:
                ## Excluding bonded interactions, but considers both inter and intra Q interactions
                G.add_edge(aa_x[bead],aa_y[bead])

        avg_edges.append(len(G.edges)/len(G.nodes))

    return avg_edges

def avg_N17_w_cutoff(tpr,xtc):
    u = Universe(tpr,xtc)
    n17 = u.select_atoms('not resname GLN and name BB')
    avg_edges = []
    for ts in u.trajectory:
        G = nx.MultiGraph()
        J = nx.Graph()
        D = np.where(distance_array(n17.positions,n17.positions,box=ts.dimensions)<7)
        x = D[0]//17
        y = D[1]//17
        for bead in range(len(x)):
            #Create a multigraph G that adds all interacting N17 domains (excluding self interactions)
            if x[bead]!=y[bead]:
                G.add_edge(x[bead],y[bead])

        # Following chunk filters out node pairs that have less than 5 edges b/w each other.
        for i,j,k in G.edges:
            #i,j are the nodes and k is the number of edges at that point
            if len(G.get_edge_data(i,j))>4:
                #If there are 5 or more interactions between 2 N17 domains, consider as interacting
                # If interacting, add to new network J
                print(ts.time/1000)
                J.add_edge(i,j)
        if nx.is_empty(J):
            avg_edges.append(0)
        else:
            avg_edges.append(len(max(nx.connected_components(J), key=len)))
    return avg_edges

def dummy_per_pep(bbp,bbm,numpeps,Q_len):
	##To split the dummy array into a 2D array shapes as (numpeps x peplen)
	## called by old_interpep_Qcontacts, old_intrapep_Qcontacts, intrapep_Qcontacts, intrapep_Qcontacts
	bbp_pep=[]
	bbm_pep=[]
	count = 0
	for i in range(numpeps):
		bbp_pep.append(bbp.atoms[count:count+Q_len])	
		bbm_pep.append(bbm.atoms[count:count+Q_len])
		count = count+Q_len
	
	return bbp_pep, bbm_pep

def remove_helical(D):
	## function to remove indices of all helical contacts
	## called by new inter and intra peptide sheet contacts
	a = np.unique(D[0])
	b = np.unique(D[1])
	new_a = []
	new_b = []
	for i, j in zip(a,b):
		if np.abs(i-j)>4:
			new_a.append(i)
			new_b.append(j)
			#print("The following qualify as helical contacts: %d,%d" %(i,j))
	return new_a, new_b		
	

def intrapep_Qcontacts(tpr, xtc, Q_len,ti = 0,tf = -1, skip = 1):
	'''
	Get thr probability of forming an intra-peptide beta-sheet contact
	Normalise by the highest number of possible intra-peptide betasheet interactions.--> peptide forming a hairpin
	Thus, normalize by Qlen*6/2, which is half the number of Qs in the system
	'''
	numpeps = 6
	u = Universe(tpr,xtc)
	intra_contacts = []
	for ts in u.trajectory[ti:tf:skip]:
		#print ts
		bbp = u.select_atoms('resname GLN and name BBp')
		bbm = u.select_atoms('resname GLN and name BBm')
		bbp_pep,bbm_pep = dummy_per_pep(bbp,bbm,numpeps,Q_len)
		sheet_contact= 0	
		for pep in range(numpeps):
			D= distance_array(bbp_pep[pep].positions, bbm_pep[pep].positions,u.dimensions)
			D = np.where(D<3)
			#print('This is peptide number: '+str(pep))
			a,b = remove_helical(D)
			split_a = np.split(a, np.where(np.diff(a)!=1)[0]+1)
			split_b = np.split(b, np.where(np.diff(a)!=1)[0]+1)
			## Count as sheet contacts only if in a chain of 4 or more	
			for i in range(len(split_a)):
				if (len(split_a[i])>3):
					temp = np.split(split_b[i],np.where(np.diff(split_b[i])!=1)[0]+1)
					temp_lens = [np.shape(temp[i])[0] for i in range(len(temp))]
					sheet_contact+=np.max(temp_lens)
		intra_contacts.append(sheet_contact)		## Normalizing by a factor of Qlen/2*6

	intra_contacts = np.array(intra_contacts)
	intra_contacts = intra_contacts/(Q_len*3)

	return intra_contacts

			
def interpep_Qcontacts(tpr, xtc, Q_len,ti = 0,tf = -1, skip = 1):
	numpeps = 6
	u = Universe(tpr,xtc)
	inter_contacts = []
	for ts in u.trajectory[ti:tf:skip]:
		#print ts
		bbp = u.select_atoms('resname GLN and name BBp')
		bbm = u.select_atoms('resname GLN and name BBm')
		bbp_pep,bbm_pep = dummy_per_pep(bbp,bbm,numpeps,Q_len)
		sheet_contact = 0	##interpeptide contacts per frame
		for i in range(numpeps):	
			for j in range(numpeps):	
				if i!=j:
					#print(i,j)
					D = distance_array(bbp_pep[i].positions, bbm_pep[j].positions,u.dimensions)
					D = np.where(D<3)
					a = np.unique(D[0])
					b = np.unique(D[1])
					#print(a,b)					
					split_a = np.split(a, np.where(np.diff(a)!=1)[0]+1)
					split_b = np.split(b, np.where(np.diff(a)!=1)[0]+1)
					for pos in range(len(split_a)):
						if split_a:
							if (len(split_a[pos])>3):
								#print(np.shape(split_a[pos])[0])		
								temp = np.split(split_b[pos],np.where(np.diff(split_b[pos])!=1)[0]+1)
								temp_lens = [np.shape(temp[k])[0] for k in range(len(temp))]
								sheet_contact+=np.max(temp_lens)
								#print sheet_contact

		inter_contacts.append(sheet_contact)## number of interpeptide contacts per frame
	inter_contacts = np.array(inter_contacts)
	inter_contacts = inter_contacts/(Q_len*6)
	
	return inter_contacts


def get_interaction_vals(tpr,xtc):

	#definitions
	u = Universe(tpr,xtc)
	status = []
	G = nx.Graph()
	n17 = u.select_atoms('name BB S1 S2 S3 S4 and not resname GLN')
	gln = u.select_atoms('resname GLN and name BB S1')
	for ts in u.trajectory[:900]:
	
		## Do N17 - N17 interaction
		
		added_n17 = False	
		len_grp1 = len(n17)/6 # Because 6 peptides
		len_grp2 = len(n17)/6
		D = np.where(distance_array(n17,n17,ts.dimensions)<7)
		
		x = D[0]//len_grp1
		y = D[1]//len_grp2
		
		for i in range(len(x)):
			# Make sure interaction is between different peptides
			if x[i]!=y[i]:
				if G.has_edge(x[i],y[i]):
					# This interaction has been encountered before, updating the list with -1 placeholder
					continue
				else:
					# If this is a new interaction, add edge and update value in desired list
					G.add_edge(x[i],y[i])
					#count = count+1
					added_n17 = True
					
		## Do GLN - GLN interaction

		added_gln = False	
		len_grp1 = len(gln)/6 # Because 6 peptides
		len_grp2 = len(gln)/6
		D = np.where(distance_array(gln,gln,ts.dimensions)<7)
		
		x = D[0]//len_grp1
		y = D[1]//len_grp2
		
		for i in range(len(x)):
			# Make sure interaction is between different peptides
			if x[i]!=y[i]:
				if G.has_edge(x[i],y[i]):
					# This interaction has been encountered before, updating the list with -1 placeholder
					continue
				else:
					# If this is a new interaction, add edge and update value in desired list
					G.add_edge(x[i],y[i])
					#count = count+1
					added_gln = True


		## Do GLN - GLN interaction

		added_cross = False	
		len_grp1 = len(n17)/6 # Because 6 peptides
		len_grp2 = len(gln)/6
		D = np.where(distance_array(n17,gln,ts.dimensions)<7)
		
		x = D[0]//len_grp1
		y = D[1]//len_grp2
		
		for i in range(len(x)):
			# Make sure interaction is between different peptides
			if x[i]!=y[i]:
				if G.has_edge(x[i],y[i]):
					# This interaction has been encountered before, updating the list with -1 placeholder
					continue
				else:
					# If this is a new interaction, add edge and update value in desired list
					G.add_edge(x[i],y[i])
					#count = count+1
					added_cross = True


		print(added_n17, added_gln, added_cross)
		if added_n17 == False and added_gln == False and added_cross == False:
			status.append(0)
			continue
		elif added_n17 == True:
			status.append(1)
			continue
		elif added_gln == True:
			status.append(2)
			continue
		elif added_cross == True:
			status.append(3)
			continue
			
	return status

########################################################################################################
## Function call loop.
## Only the function name and the location/name of the outputted file need to be changed. 
## Outputted files then read by MBAR


Qlen = [7,15,35,40,45]
T = np.arange(380,520,20) ## temp range has to be from 0 to 6 

for pep in Qlen:
	a = 0
	print('Peptide segment : '+str(pep)+'\n')
	for temp in T:
		filename_out = './domain_int/'+str(pep)+'Q.'+str(a)+'int_val.chi'
		g = open(filename_out,'w')
		tpr = str(pep)+'Q/6p_'+str(pep)+'Q_'+str(temp)+'K.tpr'
		xtc = str(pep)+'Q/6p_'+str(pep)+'Q_'+str(temp)+'K.nj.xtc'
		
		dat = get_interaction_vals(tpr,xtc)
		for i in dat:
			g.write(str(i)+'\n')
		
		g.close()
		a = a+1

