#specialized version of class SegmentTree2DFast
#(+,+) update-query
class SegmentTree2DSum:
    #(uo,ui,qo,qi,uor)=(+,0,+,0,*)
    class Node:
        def __init__(self,A,AZ):
            self.A=A
            self.AZ=AZ
    class SegmentTree1DSum:
        class Node:
            def __init__(self,V,Z):
                self.V=V
                self.Z=Z
        def __init__(self,A):
            self.N=len(A)
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
                self.nodes[nid]=self.Node(self.nodes[nl].V+self.nodes[nr].V,0)
        def U(self,x0,x1,v):
            self.UN(0,0,self.N-1,x0,x1,v)
        def UN(self,nid,nx0,nx1,x0,x1,v):
            #(nid,nx0,nx1) is the node with id number ==nid, covering [nx0,nx1]
            if x0<=nx0 and nx1<=x1:
                self.nodes[nid].Z+=v
            elif x0<=nx1 and nx0<=x1:
                m=(nx0+nx1)//2
                nl=2*nid+1
                nr=2*nid+2
                self.UN(nl,nx0,m,x0,x1,v)
                self.UN(nr,m+1,nx1,x0,x1,v)
                self.nodes[nid].V=(self.nodes[nl].V+self.nodes[nl].Z*(m+1-nx0))+(self.nodes[nr].V+self.nodes[nr].Z*(nx1-m))
        def Q(self,x0,x1):
            return self.QN(0,0,self.N-1,x0,x1)
        def QN(self,nid,nx0,nx1,x0,x1):
            if x0<=nx0 and nx1<=x1:
                return self.nodes[nid].V+self.nodes[nid].Z*(nx1+1-nx0)
            elif x0<=nx1 and nx0<=x1:
                m=(nx0+nx1)//2
                nl=2*nid+1
                nr=2*nid+2
                return (self.QN(nl,nx0,m,x0,x1)+self.QN(nr,m+1,nx1,x0,x1))+self.nodes[nid].Z*(min(nx1,x1)-max(nx0,x0)+1)
            else:
                return qi
    def __init__(self,A):
        self.N=len(A)
        self.M=len(A[0])
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
            a1D=[al[j]+ar[j] for j in range(M)]
        self.nodes[nid]=self.Node(self.SegmentTree1DSum(a1D),
                                  self.SegmentTree1DSum([0 for j in range(M)]))
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
            self.nodes[nid].A.U(y0,y1,v*(min(nx1,x1)-max(nx0,x0)+1))
    def Q(self,x0,x1,y0,y1):
        return self.QN(0,0,self.N-1,x0,x1,y0,y1)
    def QN(self,nid,nx0,nx1,x0,x1,y0,y1):
        if x0<=nx0 and nx1<=x1:
            return self.nodes[nid].A.Q(y0,y1)+self.nodes[nid].AZ.Q(y0,y1)*(nx1-nx0+1)
        elif x0<=nx1 and nx0<=x1:
            m=(nx0+nx1)//2
            nl=2*nid+1
            nr=2*nid+2
            return (self.QN(nl,nx0,m,x0,x1,y0,y1)+self.QN(nr,m+1,nx1,x0,x1,y0,y1))+self.nodes[nid].AZ.Q(y0,y1)*(min(nx1,x1)-max(nx0,x0)+1)
        else:
            return qi
