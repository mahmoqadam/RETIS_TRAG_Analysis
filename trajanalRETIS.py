######## INPUT ############################################################
TRAJdir="/Users/Mehr/Desktop/REMOVAL/008"       #input directory for trajectory file
L=13.0           #Box length 
traj4analysis=range(1,1000)
OnlyAccepted=True
##selectframes=range(10400,11524)                      #selection of frames 
centerindex=[208,209] #places an atom or the mid-point of a group of atoms in the center. Leave empty if no centering is desired
Mirrors=True         #If true the molecules overlapping the PBC will be duplicated
###### END OF INPUT #####################################################


######## DEFAULT BONDING PARAMETERS ############################################################
#following values are obtained from limiting distances where jmol gives a bond:
#i.e. dX is half the threshold distance for the X-X bond
dH=0.91/2.
dO=1.81/2.
dSi=2.85/2.
dNa=2.39/2.
###### END OF DEFAULT BONDING PARAMETERS #####################################################

########## IMPORTS #######################################################
from numpy import *
from os import system
from copy import *
#########################################################################


#########################################################################
def CHANGESTRING(file,oldstring,newstring):
 """
 Change a certain string inside an existing file
 """
 st="sed 's/"+oldstring+"/"+newstring+"/g' "+\
     file+" > XXX.tmp;mv XXX.tmp "+file
 system(st)

#########################################################################
def DISTANCE(c1,c2,L=None):
  """
  Miniumum distance between coordinates c1 and c2 in periodic system
  or non-periodic system (if L is not provided)
  """
  vector=c1-c2
  if L is not None: vector-=L*around(vector/L) #apply periodic boundaries
  d=sqrt(sum(vector*vector))
  return d

#########################################################################
# This function returns a True is at1 and at2 are connected             #
#########################################################################
def CONNECTED(at1,at2,cutoffs,L=False):
  c1,el1=at1[0],at1[1]
  cutoff1=cutoffs[el1]
  c2,el2=at2[0],at2[1]
  cutoff2=cutoffs[el2]
  d=DISTANCE(c1,c2,L)
  return d<cutoff1+cutoff2
  

#########################################################################
def EXTRACTNEIGHBORSFROMLIST(atom,leftover,cutoffs,L):
  """
  Check neighbors of atom in list leftover. Remove them from leftover and return them
  as list extract.
  """
  indexleftover=0
  extract=[]
  while indexleftover<len(leftover):
    secatom=leftover[indexleftover]
    if CONNECTED(atom,secatom,cutoffs,L):
      extract+=[secatom]
      del leftover[indexleftover]
    else:
      indexleftover+=1
  return extract,leftover

#########################################################################
def MOLECLIST(atomlist,L,cutoffs):
 """
 Deliver back a list of molecules by checking with atoms are connected 
 """
 moleclist=[]
 leftover=deepcopy(atomlist)
 while len(leftover)>0: #until all atoms are assigned to molecules
   mol=[]               #create a new molecule
   mol+=[leftover[0]]   #first atom of new molecule is first atom in list leftover
   del leftover[0]      #remove it from leftover
   iat=0                #atom-index inside molecule
   while iat<len(mol): #stop if last atom does not generate new neighbors to be added.
     atom=mol[iat]
     neighbors,leftover=EXTRACTNEIGHBORSFROMLIST(atom,leftover,cutoffs,L)
     #the above subroutine call checks which atoms from list leftover are neighbors of atom
     #the list of neighbors is returned and taken out from list leftover
     mol+=neighbors #molecule is extended by all neighbors of atom
     iat+=1 #look at the next atom in the molecule that is being generated
   moleclist+=[mol] 
 return moleclist 

#########################################################################
def MIRRORCOORDINATES(mol,L,N,Mirrors):
  firstat=mol[0]   #the first atom defines the visual position of the molecule
  mirror=[firstat] #here we will store the same  coordinates as mol but shifted across PBC
  del mol[0]
  imirror=0
  listoftrans=[]
  while len(mol)>0: #stop if all atoms are transferred to mirror 
    atom=mirror[imirror]
    neighbors,mol=EXTRACTNEIGHBORSFROMLIST(atom,mol,cutoffs,L)
    #the above subroutine call checks which atoms from list mol are neighbors of atom
    #the list of neighbors is returned and taken out from list mol 
    for ni in neighbors:
      cni=ni[0]
      cat=atom[0]
      vector=cni-cat
      trans=around(vector/L)
      cni-=trans*L
      ni[0]=cni      
      mirror+=[ni] 
      if list(trans)!= [0,0,0] and list(trans) not in  [ list(item) for item in listoftrans ]: listoftrans+=[trans]
    imirror+=1 #look at the next atom in the mirror
  #add other periodic images of molecule 
  if Mirrors:
    mol=deepcopy(mirror)
    for at in mol:
      for trans in listoftrans:
        newat=deepcopy(at)
        newat[0]+=trans*L
        newat[2]+=N
        mirror+=[newat]
  return mirror 

#########################################################################
def WRITEMOLECLIST(g,moleclist,counter,commentline):
  """
  Write info about which atoms form molecules
  """
  g.write("********** counter="+str(counter)+" ****************\n")
  g.write(commentline)
  for mol in moleclist:
    g.write("nat="+str(len(mol))+" ")
    for at in mol:
      g.write(at[1])
    for at in mol:
      g.write(" "+str(at[2]))
    if len(mol)!=3: g.write("         !!!!!!!!!!!!") 
    if len(mol)>4:  g.write(4*"!!!!!!!!!!!!")
    g.write("\n")

#########################################################################
def WRITEMIRROR2MOV(mirrorlist,framecount,f):
  """
  reorder atoms and write to file
  """
  for m in mirrorlist:
    m[0]=list(m[0])
  listatindex=map(lambda x: x[2],mirrorlist)
  reorderedmirror=zip(*sorted(zip(listatindex,mirrorlist)))[1]
  f.write(str(len(mirrorlist))+"\n")
  f.write(str(framecount)+"\n")
  for at in reorderedmirror:
    f.write(at[1]+" "+str(at[0]).replace("[","").replace("]","").replace(",","")+" "+str(at[2])+"\n") 

#########################################################################
def WRITEFRAME(N,f,g,counter,L,strcoordinates,cutoffs,centerindex,Mirrors,commentline):
 """
 write the frame 
 """
 #determine center
 shift=array([0.,0.,0.])
 for i in centerindex:
   x,y,z=strcoordinates[i-1].split()[1:4]
   xyz=array([float(x),float(y),float(z)])
   shift+=xyz
 shift/=max(1,len(centerindex)) #if centerindex is empty, then shift remains 0,0,0
 shift-=L/2.
 atomlist=[] 
 for atindex in range(N):
    element,x,y,z=strcoordinates[atindex].split()[0:4]
    xyz=array([float(x),float(y),float(z)])
    xyz-=shift #the centerpoint from the atom-indices in centerindex, will be at L/2,L/2,L/2
    xyz-=L*floor(xyz/L) #ensure that xyz is within the cube [0,0,0]:[L,L,L]
    atom=[xyz,element,atindex+1]
    atomlist+=[atom]
 moleclist=MOLECLIST(atomlist,L,cutoffs) #determine which atoms form molecules
 WRITEMOLECLIST(g,moleclist,counter,commentline)     #write molecule info to moleclist file
 mirrorlist=[] #will contain all atomistic coordinates plus some mirrored coordinates
 for mol in moleclist:
   mol2=MIRRORCOORDINATES(mol,L,N,Mirrors) 
   mirrorlist+=mol2
 WRITEMIRROR2MOV(mirrorlist,framecount,f)

#########################################################################
def GETTRAJINFO(lines,lines2):
 """
 create information-file about the trajectories
 """ 
 pathinfo=[]
 f=open("pathlength_ALL.txt","w")
 g=open("pathlength_ACC.txt","w")
 indexes_pdat=map(lambda x:[int(x.split()[0]),int(x.split()[1]),int(x.split()[2])],lines2)
 accrej=map(lambda x:x.split()[7],lines2)

 #internal subroutine of GETTRAJINFO
 def giveoutput(pathinfo):
  if pindex in indexes_pdat:
    ar=accrej[indexes_pdat.index(pindex)]
    f.write(str(numtraj)+" "+str(L)+" "+ar+" "+str(pindex)+" "+str([lstart,linecounter-1])+"\n")
    if ar=="ACC": g.write(str(numtraj)+" "+str(L)+"\n")
    toadd=[L,ar,index,[lstart,linecounter-1]]
    pathinfo+=[toadd]

 numtraj=0
 linecounter=0
 L=0
 for l in lines:
   linecounter+=1
   if "index" in l:
     words=l.split()
     index,time=words[1:4],words[4]
     index=map(lambda x:int(x),index)
     if time=="1":  #new path 
       if numtraj>0:
         giveoutput(pathinfo)
       lstart=linecounter    
       numtraj+=1
       L=1
     else:
       L+=1
     pindex=index
 linecounter+=1;giveoutput(pathinfo) #final output that is not considered by the loop because of missing new trajectory
 f.close()
 g.close()
 return pathinfo

#########################################################################
def MAKEXYZ(lines,pathinfo,L,traj4analysis,OnlyAccepted,N):
  """
  make xyz file
  """
  f=open("movie.xyz","w")
  counter=0
  for p in pathinfo:
    counter+=1
    if counter in traj4analysis:
      if OnlyAccepted and p[1]!="ACC": continue 
      print counter,p
      first,last=p[-1]
      for l in lines[first-1:last]:
        if "index" in l:
          f.write(str(N)+"\n"+"#"+str(counter)+" "+l) 
        else:
          f.write(l[9:60]+"\n")
  f.close() 
  print "plain xyz-format movie written in movie.xyz"

######## MAIN PROGRAM ###########################################################
# This part is executed after the input definitions at the top of this file     #
#################################################################################
#read all lines of TRAJECTORY.dat file
ifile=TRAJdir+"/TRAJECTORY.dat";ifile2=TRAJdir+"/PATH.dat"
print "Analysis of ",ifile
f=open(ifile,"r");lines=f.readlines();f.close()
f=open(ifile2,"r");lines2=f.readlines();f.close()
pathinfo=GETTRAJINFO(lines,lines2)
numtraj=len(pathinfo)
numalltimeslices=sum(map(lambda x:x[0],pathinfo))
N=(len(lines)/numalltimeslices)-1
print "number of atoms:", N
print "number of trajectories", numtraj
L=array([L,L,L]) #Box dimensions 3D
print "periodic box", L
MAKEXYZ(lines,pathinfo,L,traj4analysis,OnlyAccepted,N)

#from here on I use simple copy past from a previous python script that needed an xyz-file as an input
#It could be programmed more efficient, but I took the approach which is the easiest for me
#
f=open("movie.xyz","r");lines=f.readlines();f.close()
print "\nNow, further analysis on movie.xyz\n"
nframes=len(lines)/(N+2) 
print "number of frames", nframes

cutoffs={"H":dH,"O":dO,"Si":dSi,"Na":dNa}
f=open("moviePBC.xyz","w");g=open("moleclist.txt","w")
lcount=0;framecount=0
for l in range(nframes):
  ##if framecount in selectframes:
  framecoordinates=lines[lcount+2:lcount+2+N]
  commentline=lines[lcount+1]
  WRITEFRAME(N,f,g,framecount,L,framecoordinates,cutoffs,centerindex,Mirrors,commentline)
  framecount+=1 
  lcount=lcount+2+N
f.close();g.close()
print "molecular information in 'moleclist.txt'"
print "new movie file in moviePBC.xyz" 

st="cp jmol.base jmol.run";system(st)
#Adjust jmol.run
newstring=str(L).replace("[","").replace("]","").replace(",","")
CHANGESTRING("jmol.run","boxdimensions",newstring)
newstring=""
for aa in cutoffs.keys():
  for bb in cutoffs.keys():
    newstring+="connect "+str(cutoffs[aa]+cutoffs[bb])+" (_"+aa+") (_"+bb+")\\\n "
CHANGESTRING("jmol.run","connections",newstring)
print ""
print "Some visualization defaults (like showing hydrogen bonds) are provided as a script for jmol in jmol.run"
print "If you have jmol installed as online-command, just type\n" 
print "jmol jmol.run\n"
print "If jmol is installed but not as command-line: execute jmol and load jmol.run from pull-down menu: File-> Open -> select file 'jmol.run'\n"
print "movie.xyz can be used for other molecular visualization programs as well, but some might not work with the"
print "option Mirrors since it generates a xyz-trajectory file with fluctuating number of atoms."
print "In that case use Mirrors=False as input that can be provided in the top lines of makemovfromxyz.py"
