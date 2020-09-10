import numpy as np

class SegmentTreeFast2D:
    #only applies if qo(uo(a,v),b)=uo(qo(a,b),v)
    #def uor(v,0)=ui, uor(v,k)=uo(uor(v,k-1),v)
    class Node:
        def __init__(self,A,AZ):
            self.A=A
            self.AZ=AZ
    def __init__(self,A,uo,ui,qo,qi,uor):
        self.N=len(A)
        self.M=len(A[0])
        self.uo=uo #update operation
        self.ui=ui #update identity (can be created artificially)
        self.qo=qo #query operation
        self.qi=qi #query identity
        self.uor=uor
        self.F=(lambda a,v,sz:self.uo(a,self.uor(v,sz)))
        self.nodes=[None for i in range(4*self.N)]
        self.build(0,0,self.N-1,A)
    def build(self,nid,nx0,nx1,A):
        a1D=[]
        if nx0==nx1:
            a1D=A[nx0].copy()
        else:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            al=self.build(nl,nx0,m,A)
            ar=self.build(nr,m+1,nx1,A)
            a1D=[self.qo(al[j],ar[j]) for j in range(M)]
        self.nodes[nid]=self.Node(SegmentTree1D(a1D,uo,ui,qo,qi,F),
                                  SegmentTree1D([ui for j in range(M)],uo,ui,uo,ui,F))
        #print("[{},{}]={}".format(nx0,nx1,a1D))
        return a1D
    def U(self,x0,x1,y0,y1,v):
        self.UN(0,0,self.N-1,x0,x1,y0,y1,v)
    def UN(self,nid,nx0,nx1,x0,x1,y0,y1,v):
        if x0<=nx0 and nx1<=x1:
            self.nodes[nid].AZ.U(y0,y1,v)
        elif x0<=nx1 and nx0<=x1:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            self.UN(nl,nx0,m,x0,x1,y0,y1,v)
            self.UN(nr,m+1,nx1,x0,x1,y0,y1,v)
            self.nodes[nid].A.U(y0,y1,self.uor(v,min(nx1,x1)-max(nx0,x0)+1))
    def Q(self,x0,x1,y0,y1):
        return self.QN(0,0,self.N-1,x0,x1,y0,y1)
    def QN(self,nid,nx0,nx1,x0,x1,y0,y1):
        if x0<=nx0 and nx1<=x1:
            return self.uo(self.nodes[nid].A.Q(y0,y1),self.uor(self.nodes[nid].AZ.Q(y0,y1),nx1-nx0+1))
        elif x0<=nx1 and nx0<=x1:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            return self.uo(self.qo(self.QN(nl,nx0,m,x0,x1,y0,y1),self.QN(nr,m+1,nx1,x0,x1,y0,y1))
                           ,self.uor(self.nodes[nid].AZ.Q(y0,y1),min(nx1,x1)-max(nx0,x0)+1))
        else:
            return qi

class SegmentTree1D:
    class Node:
        def __init__(self,V,Z):
            self.V=V
            self.Z=Z
    def __init__(self,A,uo,ui,qo,qi,F):
        self.N=len(A)
        self.uo=uo #update operation
        self.ui=ui #update identity (can be created artificially)
        self.qo=qo #query operation
        self.qi=qi #query identity
        self.F=F #special function for updating every element in a set
        self.nodes=[None for i in range(4*self.N)]
        self.build(0,0,self.N-1,A)
    def build(self,nid,nx0,nx1,A):
        if nx0==nx1:
            self.nodes[nid]=self.Node(A[nx0],ui)
        else:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            self.build(nl,nx0,m,A)
            self.build(nr,m+1,nx1,A)
            self.nodes[nid]=self.Node(self.qo(self.nodes[nl].V,self.nodes[nr].V),ui)
    def U(self,x0,x1,v):
        self.UN(0,0,self.N-1,x0,x1,v)
    def UN(self,nid,nx0,nx1,x0,x1,v):
        #(nid,nx0,nx1) is the node with id number ==nid, covering [nx0,nx1]
        if x0<=nx0 and nx1<=x1:
            self.nodes[nid].Z=self.uo(self.nodes[nid].Z,v)
        elif x0<=nx1 and nx0<=x1:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            self.UN(nl,nx0,m,x0,x1,v)
            self.UN(nr,m+1,nx1,x0,x1,v)
            self.nodes[nid].V=qo(self.F(self.nodes[nl].V,self.nodes[nl].Z,m+1-nx0),
                                    self.F(self.nodes[nr].V,self.nodes[nr].Z,nx1-m))
    def Q(self,x0,x1):
        return self.QN(0,0,self.N-1,x0,x1)
    def QN(self,nid,nx0,nx1,x0,x1):
        if x0<=nx0 and nx1<=x1:
            return self.F(self.nodes[nid].V,
                          self.nodes[nid].Z,nx1+1-nx0)
        elif x0<=nx1 and nx0<=x1:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            return self.F(self.qo(self.QN(nl,nx0,m,x0,x1),self.QN(nr,m+1,nx1,x0,x1)),
                          self.nodes[nid].Z,min(nx1,x1)-max(nx0,x0)+1)
        else:
            return qi
    def arr(self):
        A=[None for i in range(self.N)]
        self.arrHelp(A,0,0,self.N-1,ui)
        return A
    def arrHelp(self,A,nid,nx0,nx1,lazy):
        nlazy=self.uo(lazy,self.nodes[nid].Z)
        if nx0==nx1:
            A[nx0]=self.uo(self.nodes[nid].V,nlazy)
        else:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            self.arrHelp(A,nl,nx0,m,nlazy)
            self.arrHelp(A,nr,m+1,nx1,nlazy)

class BruteForce2D:
    def __init__(self,A,uo,ui,qo,qi):
        self.N=len(A)
        self.M=len(A[0])
        self.uo=uo #update operation
        self.ui=ui #update identity (can be created artificially)
        self.qo=qo #query operation
        self.qi=qi #query identity
        self.A=[[v for v in row] for row in A]
    def U(self,x0,x1,y0,y1,v):
        for i in range(x0,x1+1):
            for j in range(y0,y1+1):
                self.A[i][j]=self.uo(self.A[i][j],v)
    def Q(self,x0,x1,y0,y1):
        out=self.qi
        for i in range(x0,x1+1):
            for j in range(y0,y1+1):
                out=qo(out,self.A[i][j])
        return out

#speedtest
uo=lambda a,b:a+b
ui=0
qo=uo
qi=ui
uor=lambda v,k:v*k
N=1000
M=1000
UQ=10000
MV=1_000_000
print("N={},M={},UQ={},MV={}".format(N,M,UQ,MV))
testBF=True
rnd.seed(1)
A=[[rnd.randint(-MV,MV) for j in range(M)] for i in range(N)]
'''for r in A:
    print(r)
print('')'''
actions=[]
for i in range(UQ):
    t=rnd.randrange(2)
    i0=rnd.randrange(N)
    i1=rnd.randrange(N)
    x0=min(i0,i1)
    x1=max(i0,i1)
    i0=rnd.randrange(M)
    i1=rnd.randrange(M)
    y0=min(i0,i1)
    y1=max(i0,i1)
    act=[t,x0,x1,y0,y1]
    if t==0:
        act.append(rnd.randint(-MV,MV))
    actions.append(act)

def results2D(actions,ds):
    out=[]
    for a in actions:
        if a[0]==0:
            ds.U(a[1],a[2],a[3],a[4],a[5])
        else:
            out.append(ds.Q(a[1],a[2],a[3],a[4]))
    return out

if testBF:
    stt=time.time()
    bf=BruteForce2D(A,uo,ui,qo,qi)
    print('bf init time='+str(time.time()-stt))
    stt=time.time()
    bfRet=results2D(actions,bf)
    print('bf run time='+str(time.time()-stt))

stt=time.time()
st=SegmentTreeFast2D(A,uo,ui,qo,qi,uor)
print('st init time='+str(time.time()-stt))
stt=time.time()
stRet=results2D(actions,st)
print('st run time='+str(time.time()-stt))

if testBF:
    print('st '+('good' if bfRet==stRet else 'mismatch'))
